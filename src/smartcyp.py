"""SMARTCyp turnover scores as molecular features, via the CDK port of the Java tool.

What this is and why it is not redundant with what we have. SMARTCyp predicts **sites of
metabolism** from precomputed DFT activation energies plus accessibility -- that is *turnover*,
how fast the enzyme processes the molecule. Our labels are *inhibition*. The two are related and
distinct, and the distinction is the reason to try: item 156 measured that the model has already
learned the binding pharmacophore from the fingerprint on its own -- pyrimidine, N-aryl azole,
pyridine nitrogen, all sp2 nitrogen with a lone pair for the haem iron -- so writing binding
chemistry by hand adds nothing. Turnover is a different phenomenon and arrives from outside as a
finished number rather than being estimated from inside our own data, where item 118 found it
underdetermined.

Provenance, stated because it is not clean. The official distribution sites are gone: smartcyp.sdu.dk
does not resolve and the Copenhagen page returns 503. The only accessible source is the Chemistry
Development Kit organisation's port, `github.com/cdk/smartcyp`, LGPL-3.0, whose own README says
**"It has not been tested yet."** It is version 2.5.0-SNAPSHOT based on SMARTCyp 2.4.2, not 3.0.
Built here with Maven against Java 8 bytecode in thirty seconds; the fat jar is 10.8 MB.

What the tool emits, per atom:

    Score, Ranking, Energy       общая модель
    2D6score, 2D6ranking, N+Dist   модель CYP2D6 и расстояние до протонированного азота
    2Cscore, 2Cranking, COODist    модель CYP2C9 и расстояние до карбоксилата
    Relative Span, Span2End, 2DSASA

`N+Dist` and `COODist` are the same geometry our own mechanistic block computes topologically for
CYP2D6 (`topo_bN_to_arom_min`) and the acid block for CYP2C9 -- derived independently, by a
different method. That makes the ablation sharper than a plain addition: SMARTCyp's CYP2D6 columns
go on top of a block already worth +0.0454 on that enzyme, so if they overlap they will show it,
while CYP1A2 and CYP3A4 have nothing for them to overlap with.

A cheap precondition was checked before writing this: on 67 per cent of a thirty-molecule sample
the CYP2D6 model ranks a **different** atom first than the general model, so the per-enzyme models
are not a monotone rescaling of the general one and there is something to add.

Aggregation. The tool is per atom and the model needs per molecule, so each score is reduced by
its minimum -- the most labile site, which is what "will this be metabolised" means -- plus the
mean over atoms and the count of sites below a threshold, which distinguishes one soft spot from
many.

**The molecule index is the trap.** SMARTCyp numbers its output 1..N in input order, so a molecule
CDK cannot parse would shift every subsequent row silently -- the same failure `data/rows.csv`
exists to prevent. Every chunk therefore asserts that its highest index equals its input length,
and a chunk that fails is retried one molecule at a time so the failure is isolated to the rows it
belongs to rather than shifting the ones after it.

**Two defences against CDK, both needed.** The first full-set run died with
`DeduceBondSystemTool ... Timed out after 100 seconds` -- CDK's legacy aromatic bond-order
deduction is exponential on some ring systems, and one molecule took the whole batch down. So the
input is written as **Kekule** SMILES, with explicit alternating bond orders, which is what that
routine exists to recover; given it, there is nothing to deduce. And the run is chunked, so a
molecule that defeats it anyway costs its chunk rather than the set. Molecules that fail even
alone are recorded as NaN rows and counted in the output.

Writes data/smartcyp.npz with the train and test blocks and the column names. ~5 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D

import argparse
import glob
import os
import subprocess
import tempfile
import time

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

JAVA = "/Library/Java/JavaVirtualMachines/microsoft-11.jdk/Contents/Home/bin/java"
JAR = ("/private/tmp/claude-501/-Users-nikitapolomosnov-PycharmProjects-CYPChallenge/"
       "5e67ea47-cb5e-4c94-b84b-5c8193beeec5/scratchpad/smartcyp/target/smartcyp.jar")
SCORES = ["Score", "2D6score", "2Cscore"]
GEOM = ["N+Dist", "COODist", "2DSASA", "Relative Span", "Span2End"]
SOFT = 40.0        # порог «мягкого места»: ниже него сайт считается лабильным


def kekule(smi):
    """Кекулевская форма: у CDK нечего дедуцировать, и экспоненциальная ветвь не запускается."""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    try:
        m = Chem.Mol(m)
        Chem.Kekulize(m, clearAromaticFlags=True)
        return Chem.MolToSmiles(m, kekuleSmiles=True)
    except Exception:
        return None


def _one_call(smiles, java, jar, timeout):
    """Один вызов jar на список. Возвращает таблицу или None, если упал/разъехался."""
    d = tempfile.mkdtemp(prefix="scyp_")
    inp = os.path.join(d, "in.smi")
    open(inp, "w").write("\n".join(smiles) + "\n")
    try:
        subprocess.run([java, "-jar", jar, "-nohtml", "-outputdir", d, inp],
                       capture_output=True, text=True, timeout=timeout)
    except subprocess.TimeoutExpired:
        return None
    got = sorted(glob.glob(os.path.join(d, "SMARTCyp_Results*.csv")))
    if not got:
        return None
    df = pd.read_csv(got[-1])
    if df.empty or int(df.Molecule.max()) != len(smiles) \
            or df.Molecule.nunique() != len(smiles):
        return None
    return df


def run(smiles, java, jar, tag, chunk=250):
    """Пакетно кусками, с одиночным повтором упавшего куска. Индексы --- глобальные."""
    t0 = time.time()
    kek = [kekule(s) for s in smiles]
    n_bad_smi = sum(v is None for v in kek)
    parts, dropped = [], []
    for a0 in range(0, len(smiles), chunk):
        b0 = min(a0 + chunk, len(smiles))
        idx = [i for i in range(a0, b0) if kek[i] is not None]
        df = _one_call([kek[i] for i in idx], java, jar, 1800) if idx else None
        if df is not None:
            loc = {j + 1: idx[j] for j in range(len(idx))}
            df = df.assign(Molecule=df.Molecule.map(loc) + 1)
            parts.append(df)
            continue
        # Кусок упал --- изолируем виновника, а не теряем двести пятьдесят строк.
        print(f"    кусок {a0}-{b0} упал, разбираю по одной", flush=True)
        for i in idx:
            d1 = _one_call([kek[i]], java, jar, 120)
            if d1 is None:
                dropped.append(i)
            else:
                parts.append(d1.assign(Molecule=i + 1))
    df = pd.concat(parts, ignore_index=True)
    print(f"  {tag}: {len(smiles)} на входе, RDKit не кекулизовал {n_bad_smi}, "
          f"CDK не осилил {len(dropped)}, получено {df.Molecule.nunique()} молекул, "
          f"{len(df)} строк  ({time.time()-t0:.0f} с)")
    if dropped:
        print(f"    выброшены индексы: {dropped[:12]}{' ...' if len(dropped) > 12 else ''}")
    return df


def aggregate(df, n):
    """По атомам -> по молекуле. Минимум --- самый лабильный сайт."""
    cols, names = [], []
    g = df.groupby("Molecule")
    for c in SCORES:
        for how, fn in (("min", "min"), ("mean", "mean")):
            v = g[c].agg(fn).reindex(range(1, n + 1)).to_numpy(float)
            cols.append(v); names.append(f"{c}_{how}")
        v = df.assign(f=df[c] < SOFT).groupby("Molecule").f.sum() \
              .reindex(range(1, n + 1)).fillna(0).to_numpy(float)
        cols.append(v); names.append(f"{c}_n_soft")
    for c in GEOM:
        v = g[c].min().reindex(range(1, n + 1)).to_numpy(float)
        cols.append(v); names.append(f"{c}_min")
    # Расходятся ли поферментные модели с общей: 1 если первый атом разный.
    top = df[df.Ranking == 1].groupby("Molecule").Atom.first().reindex(range(1, n + 1))
    for c, r in (("2D6", "2D6ranking"), ("2C", "2Cranking")):
        t = df[df[r] == 1].groupby("Molecule").Atom.first().reindex(range(1, n + 1))
        cols.append((top != t).to_numpy(float)); names.append(f"{c}_top_differs")
    return np.column_stack(cols).astype(np.float32), names


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--java", default=JAVA)
    ap.add_argument("--jar", default=JAR)
    ap.add_argument("--out", default=D + "smartcyp.npz")
    a = ap.parse_args()
    if not os.path.exists(a.jar):
        raise SystemExit(f"нет jar: {a.jar}\nсобрать: git clone https://github.com/cdk/smartcyp "
                         "&& cd smartcyp && mvn -B clean package")

    rows = pd.read_csv(D + "rows.csv")
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    print("SMARTCyp 2.5.0-SNAPSHOT (порт CDK от SMARTCyp 2.4.2, LGPL-3.0)")
    print("README порта: «It has not been tested yet» --- проверяется здесь, а не принимается.\n")

    out = {}
    for tag, smi in (("train", list(rows.SMILES)), ("test", list(te.SMILES))):
        df = run(smi, a.java, a.jar, tag)
        X, names = aggregate(df, len(smi))
        out[tag] = X
        out["names"] = np.array(names, dtype=object)
        bad = int(np.isnan(X).sum())
        print(f"    {X.shape[1]} колонок, NaN {bad} "
              f"({100*bad/X.size:.2f} %) --- их даёт молекула без подходящих сайтов")

    np.savez_compressed(a.out, **out)
    print(f"\nсохранено: {a.out}")
    nm = list(out["names"])
    Xtr = out["train"]
    print(f"\n{'колонка':22s} {'медиана':>10s} {'sd':>10s} {'доля NaN':>10s}")
    for i, n in enumerate(nm):
        v = Xtr[:, i]
        print(f"{n:22s} {np.nanmedian(v):10.2f} {np.nanstd(v):10.2f} "
              f"{100*np.isnan(v).mean():9.1f} %")
    print("""
Что дальше. Абляция ставится ПОФЕРМЕНТНО и читается против поферментных полов пункта 165
(1A2 0.0061, 2C9 0.0071, 2D6 0.0049, 3A4 0.0033), а не против макро --- пункт 167 показал,
что макро топит поферментные эффекты.

Предрегистрация, потому что у двух ферментов уже есть чем перекрыться, а у двух нет.
Колонки 2D6 ложатся поверх блока, который на этом ферменте стоит +0.0454, а COODist поверх
кислотного блока, который на 2C9 намерен в -0.0023. Если SMARTCyp там ничего не добавит ---
это совпадение с уже имеющимся, а не провал инструмента. Решают CYP1A2 и CYP3A4, где
перекрываться не с чем и где, по пункту 167, не работает вообще ничто.""")


if __name__ == "__main__":
    main()
