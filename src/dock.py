"""Docking into the four CYP cavities: the (ligand, cavity) PAIR quantity, and the last open lever.

Why this exists and what it is not. Item 168 separates two things: descriptors OF THE ENZYME are
closed by arithmetic (a one-hot is a sufficient statistic for four enzymes all seen in training),
while the INTERACTION form -- a quantity of the (ligand, cavity) pair -- is the one class that is
not a re-encoding of the SMILES. Item 189 collected six ligand-only blocks that all returned zero.
A docked pose depends on the cavity, so it is outside that class by construction.

THE PRIOR IS LOW AND WAS LOWERED TWICE, which is why this file records the cost honestly rather
than selling the idea. Item 179: hand-built active-site blocks bought nothing, and CYP2D6's own
salt-bridge angles made it WORSE (-0.0094). Item 296: the cheap proxy for this very campaign -- the
co-crystal overlay contrast -- came in at +0.0017 on CYP1A2 against a floor of 0.0061, a quarter of
the bar, and on CYP2D6 it did not beat its own wrong-isoform control at all. So a pose has to supply
what neither a guessed cavity nor an overlay onto one reference ligand could.

THE FEATURE IS THE CONTRAST, not the four scores. The common part of the four affinities is how
large and how greasy the molecule is, and the ligand block already carries that in 2295 columns.
Centring each molecule across the four cavities cancels it and leaves complementarity to a
particular pocket. The contrast is also the wrong-isoform control built into the feature itself:
if the own-cavity column carried nothing but bulk, centring would annihilate it -- which is exactly
what happened to CYP2C9 and CYP3A4 in item 296 (+0.077 -> +0.002 and +0.082 -> +0.002).

MEASURED COST, because "docking is expensive" is not a number. smina (master:dc3dfab, Vina 1.1.2),
exhaustiveness 8, one CPU per process, on a six-molecule set spanning our size distribution
(p10..p99, 21-34 heavy atoms): 38.9 s/molecule on 2HI4 and 26.2 s/molecule on 3NXU. At ~32 s over
5652 ligands x 4 cavities = 22608 runs that is about 204 CPU-hours, or roughly a day of wall clock
on this machine. An earlier probe on a 12-heavy-atom molecule gave 5.5 s and was unrepresentative
by a factor of six; the size distribution here is narrow and high (p25=23, p50=24, p95=30).

RESUMABILITY IS NOT OPTIONAL at that length. Work is split into chunks; each (receptor, chunk) writes
its own output and is skipped if already complete. A crash at hour twenty costs one chunk, not the
run. Progress is appended to results/logs/dock_progress.jsonl so the state is inspectable while it
runs.

RECEPTOR PREPARATION, and the one place a silent catastrophe was avoided. 4WNV carries four protein
copies and 3NXU two; our co-crystal ligand SDFs had to be matched to the right one or the autobox
would have landed in empty space. Measured: all four SDFs coincide with chain A at 0.00 A, the other
copies sitting 50-89 A away. So each receptor is chain A plus its own heme, with waters, glycerol,
DMSO and ions stripped. The heme stays: it is part of the site, and without it there is a hole where
the iron should be.

Reads data/lig3d_{train,test}.sdf (built by the 3D pass, names carry the rows.csv index),
data/receptors/*_rec.pdb, data/cocrystal/*.sdf. Writes data/dock.npz with `train` (4905x4),
`test` (750x4), `names`, and `ok_train`/`ok_test` masks -- NaN where a ligand had no 3D structure
(three of 4905) or a docking run failed, so the matrix stays rectangular and aligned to rows.csv.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES

import argparse, json, os, subprocess, time
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit import RDLogger

RDLogger.DisableLog("rdApp.*")

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
# enzyme -> (receptor pdb id, co-crystal ligand sdf). Chain A in every case, verified at 0.00 A.
CAVITY = {
    "CYP1A2": ("2HI4", "CYP1A2_2hi4_BHF.sdf"),
    "CYP2C9": ("1R9O", "CYP2C9_1r9o_FLP.sdf"),
    "CYP2D6": ("4WNV", "CYP2D6_4wnv_QI9.sdf"),
    "CYP3A4": ("3NXU", "CYP3A4_3nxu_RIT.sdf"),
}
SMINA = os.path.expanduser("~/smina")


def split_chunks(sdf_path, outdir, tag, size):
    """Split a multi-model SDF into chunk files, preserving the _Name row index."""
    _pl.Path(outdir).mkdir(parents=True, exist_ok=True)
    chunks, buf, k = [], [], 0
    supp = Chem.ForwardSDMolSupplier(sdf_path, removeHs=False)
    for m in supp:
        if m is None:
            continue
        buf.append(m)
        if len(buf) >= size:
            p = f"{outdir}/{tag}_{k:04d}.sdf"
            if not os.path.exists(p):
                w = Chem.SDWriter(p)
                for x in buf:
                    w.write(x)
                w.close()
            chunks.append(p); buf = []; k += 1
    if buf:
        p = f"{outdir}/{tag}_{k:04d}.sdf"
        if not os.path.exists(p):
            w = Chem.SDWriter(p)
            for x in buf:
                w.write(x)
            w.close()
        chunks.append(p)
    return chunks


def chunk_done(out_sdf, n_expected):
    """A chunk counts as done only if its output holds as many poses as the input had ligands."""
    if not os.path.exists(out_sdf):
        return False
    try:
        with open(out_sdf) as f:
            return sum(1 for line in f if line.startswith("$$$$")) >= n_expected
    except OSError:
        return False


def parse_scores(out_sdf):
    """name -> best affinity. smina writes minimizedAffinity on each pose."""
    out = {}
    for m in Chem.ForwardSDMolSupplier(out_sdf, removeHs=False):
        if m is None:
            continue
        nm = m.GetProp("_Name") if m.HasProp("_Name") else None
        if nm is None:
            continue
        for key in ("minimizedAffinity", "CNNaffinity", "affinity"):
            if m.HasProp(key):
                v = float(m.GetProp(key))
                # keep the best (most negative) if several poses share a name
                if nm not in out or v < out[nm]:
                    out[nm] = v
                break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=7, help="одновременных процессов smina, по одному ядру")
    ap.add_argument("--chunk", type=int, default=60, help="лигандов в чанке; меньше = дешевле потеря при падении")
    ap.add_argument("--exhaustiveness", type=int, default=8)
    ap.add_argument("--scoring", default="vina", help="мерено на vina; смена scoring --- методологическое изменение")
    ap.add_argument("--enzymes", default=",".join(CYPS))
    ap.add_argument("--assemble-only", action="store_true")
    a = ap.parse_args()
    enzymes = [e for e in a.enzymes.split(",") if e.strip()]

    work = RES + "dock/"
    _pl.Path(work).mkdir(parents=True, exist_ok=True)
    prog = RES + "logs/dock_progress.jsonl"
    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)

    chunks = {}
    for tag, sdf in (("tr", D + "lig3d_train.sdf"), ("te", D + "lig3d_test.sdf")):
        chunks[tag] = split_chunks(sdf, work + "chunks", tag, a.chunk)
        print(f"{tag}: {len(chunks[tag])} чанков по <= {a.chunk}", flush=True)

    jobs = []
    for enz in enzymes:
        pdb, lig = CAVITY[enz]
        for tag in ("tr", "te"):
            for c in chunks[tag]:
                base = os.path.basename(c).replace(".sdf", "")
                out = f"{work}{pdb}_{base}.sdf"
                jobs.append((enz, pdb, lig, c, out))
    print(f"всего заданий: {len(jobs)}  (пропущу уже готовые)", flush=True)

    if not a.assemble_only:
        # сколько лигандов в каждом чанке --- для проверки полноты вывода
        n_in = {}
        for tag in chunks:
            for c in chunks[tag]:
                with open(c) as f:
                    n_in[c] = sum(1 for line in f if line.startswith("$$$$"))

        todo = [j for j in jobs if not chunk_done(j[4], n_in[j[3]])]
        print(f"к запуску: {len(todo)}, уже готово: {len(jobs) - len(todo)}", flush=True)

        running, t_start, done = [], time.time(), 0
        while todo or running:
            while todo and len(running) < a.jobs:
                enz, pdb, lig, cin, out = todo.pop(0)
                cmd = [SMINA, "-r", f"{D}receptors/{pdb}_rec.pdb", "-l", cin,
                       "--autobox_ligand", f"{D}cocrystal/{lig}", "--autobox_add", "4",
                       "--scoring", a.scoring, "--exhaustiveness", str(a.exhaustiveness),
                       "--num_modes", "1", "--cpu", "1", "--seed", "42", "-o", out]
                p = subprocess.Popen(cmd, stdout=subprocess.DEVNULL, stderr=subprocess.STDOUT)
                running.append((p, enz, pdb, cin, out, time.time()))
            time.sleep(2)
            for item in list(running):
                p, enz, pdb, cin, out, t0 = item
                if p.poll() is None:
                    continue
                running.remove(item)
                done += 1
                rec = {"enz": enz, "pdb": pdb, "chunk": os.path.basename(cin),
                       "rc": p.returncode, "sec": round(time.time() - t0, 1),
                       "done": done, "left": len(todo) + len(running)}
                with open(prog, "a") as f:
                    f.write(json.dumps(rec, ensure_ascii=False) + "\n")
                el = time.time() - t_start
                rate = done / el if el > 0 else 0
                eta = (len(todo) + len(running)) / rate / 3600 if rate > 0 else float("nan")
                print(f"  [{done}/{len(jobs)}] {pdb} {os.path.basename(cin)} rc={p.returncode} "
                      f"{rec['sec']:.0f}s  осталось {rec['left']}  ETA {eta:.1f} ч", flush=True)

    # ---- сборка ------------------------------------------------------------------------------
    rows = pd.read_csv(D + "rows.csv")
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    A = {"tr": np.full((len(rows), len(CYPS)), np.nan),
         "te": np.full((len(te), len(CYPS)), np.nan)}
    for enz in enzymes:
        e = CYPS.index(enz)
        pdb, _ = CAVITY[enz]
        for tag in ("tr", "te"):
            for c in chunks[tag]:
                out = f"{work}{pdb}_{os.path.basename(c).replace('.sdf','')}.sdf"
                if not os.path.exists(out):
                    continue
                for nm, v in parse_scores(out).items():
                    if "_" not in nm:
                        continue
                    t, i = nm.split("_", 1)
                    if t == tag:
                        A[tag][int(i), e] = v
    cov = {c: (int(np.isfinite(A['tr'][:, i]).sum()), int(np.isfinite(A['te'][:, i]).sum()))
           for i, c in enumerate(CYPS)}
    print("\nпокрытие (train, test) по ферментам:", cov, flush=True)
    ok_tr = np.isfinite(A["tr"]).all(1)
    ok_te = np.isfinite(A["te"]).all(1)
    print(f"полных строк: train {int(ok_tr.sum())}/{len(rows)}, test {int(ok_te.sum())}/{len(te)}",
          flush=True)
    np.savez_compressed(D + "dock.npz", train=A["tr"], test=A["te"],
                        names=np.array(CYPS), ok_train=ok_tr, ok_test=ok_te)
    print(f"сохранено: {D}dock.npz", flush=True)
    print("""
Как читать. Это СЫРЫЕ аффинности по четырём полостям, не признак. Признак --- КОНТРАСТ:
вычесть из каждой строки её среднее по четырём полостям. Общая часть (размер, липофильность) у нас
уже есть, и она же несёт объёмный конфаунд CYP3A4 (пункт 199: корреляция с числом атомов +0.283
даже после центрирования, потому что полость огромна, а эталон --- ритонавир на 98 атомов).
Оценивать блок только предрегистрированной абляцией по рангу после аффинной пары, на свежих сидах,
и помнить, что контраст между полостями и есть встроенный контроль неверной изоформой.""")


if __name__ == "__main__":
    main()
