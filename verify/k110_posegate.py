"""Pre-upload validity gate for the structure track: what PoseBusters would zero.

WHY THIS IS THE FIRST THING BUILT FOR THAT TRACK, before any pose of our own exists. The
organisers' scorer does not merely penalise an invalid pose, it ERASES it: per
`evaluation.evaluate_predictions.score_single_structure`, a model copy failing more than
`POSEBUSTERS_MAX_FAILURES` checks has LDDT-PLI set to 0, LDDT-LP set to 0 AND BiSyRMSD set to
`BISYRMSD_NAN_PENALTY`, and that override happens BEFORE the best-pair sort, so a bad copy cannot
be rescued by its OST scores. `coverage` is computed from LDDT-PLI.notna() BEFORE the fill, so a
zeroed structure still counts as covered: **Coverage 1.0 on the leaderboard does not mean twenty
structures scored, only that twenty parsed.** The board is scored by a plain mean over the twenty
compounds, so one erased pose that would have scored 0.42 costs 0.021 of the mean — and on
1 October the whole top of that board, rank 1 at 0.4252 down to rank 5 at 0.3921, spans 0.033.
About one and a half erased poses.

THE THRESHOLD CANNOT BE READ FROM THE PUBLIC CODE, and this file is built around that rather than
around a guess. `evaluation/evaluate_predictions.py` imports POSEBUSTERS_MAX_FAILURES,
BISYRMSD_NAN_PENALTY and STRUCTURE_METRICS from `.config`, and the shipped `evaluation/config.py`
defines none of the three -- proven by execution, not by grep: the import raises
`ImportError: cannot import name 'BISYRMSD_NAN_PENALTY'`. Upstream commit 832ae6d is HEAD, so this
is the current published state. Therefore this gate reports the DISTRIBUTION of failed-check counts
per pose, never a verdict against an assumed tolerance. A pose failing many checks is erased under
any plausible threshold; a pose failing one is a question nobody outside the organisers can answer.

THE CONTROL IS NOT OPTIONAL. "All twenty clean" is the shape of a query that could not have failed
-- this project's most frequent defect, and the reason item 304 lost 2632 of 4902 poses to an
unchecked RDKit None. So the run plants a 21st pose: a real ligand translated onto the heme iron.
If the harness does not FLAG it, the run says so and refuses its own result. Each pose is also
checked before busting -- a residue named exactly LIG must exist and carry the heavy-atom count the
SMILES implies -- so a silently empty selection cannot masquerade as a clean pass.

INTERPRETER. PoseBusters is not in the project environment and must not be added to it: the
project caps scikit-learn below 1.9 and `pyproject.toml` is load-bearing for reproducing
`results/preds/oof.json` bit for bit. Build a scratch environment instead and run this file with
it. The organisers' own bond-order helper is used rather than RDKit's
`AllChem.AssignBondOrdersFromTemplate`, because theirs tries every topological match -- their
docstring says the plain one can raise a spurious valence error or silently mis-assign a locally
symmetric ring.

    uv venv /tmp/pbenv
    uv pip install --python /tmp/pbenv/bin/python posebusters loguru
    /tmp/pbenv/bin/python verify/k110_posegate.py --cifdir <dir of *.cif> --eval /tmp/ost_check

Exit: 0 the gate ran AND the planted pose was flagged, 1 the planted pose was missed (result void),
2 inputs missing. Writes results/logs/k110_posegate.json.
"""

import argparse
import json
import pathlib
import sys
import urllib.request

LIGCSV = ("https://huggingface.co/datasets/openadmet/cyp-challenge-train-test/"
          "resolve/main/cyp-challenge-TEST-BLINDED_structures.csv")


def parse_cif(path):
    """Atoms of a Boltz-style mmCIF as (comp_id, element, x, y, z, atom_name, seq_id)."""
    lines = pathlib.Path(path).read_text().splitlines()
    start = None
    for k, line in enumerate(lines):
        if line.strip() == "loop_" and k + 1 < len(lines) \
                and lines[k + 1].strip().startswith("_atom_site."):
            start = k
            break
    if start is None:
        raise SystemExit(f"{path}: нет цикла _atom_site")
    cols, j = [], start + 1
    while lines[j].strip().startswith("_atom_site."):
        cols.append(lines[j].strip())
        j += 1
    idx = {c: n for n, c in enumerate(cols)}
    need = ["label_comp_id", "type_symbol", "Cartn_x", "Cartn_y", "Cartn_z",
            "label_atom_id", "auth_seq_id"]
    for n in need:
        if "_atom_site." + n not in idx:
            raise SystemExit(f"{path}: нет колонки {n}")
    out = []
    while j < len(lines) and lines[j].strip() and not lines[j].startswith("#") \
            and not lines[j].strip().startswith("loop_"):
        p = lines[j].split()
        if len(p) == len(cols):
            g = lambda n: p[idx["_atom_site." + n]]
            out.append((g("label_comp_id"), g("type_symbol").upper(),
                        float(g("Cartn_x")), float(g("Cartn_y")), float(g("Cartn_z")),
                        g("label_atom_id"), g("auth_seq_id")))
        j += 1
    return out


def pdb_block(atoms, hetero_for):
    """Fixed-column PDB. Serial stays inside five characters by construction; item 302 is the
    reason that is stated rather than assumed -- a six-digit serial silently destroys the columns
    every downstream fixed-format reader depends on."""
    if len(atoms) > 99999:
        raise SystemExit(f"атомов {len(atoms)}, серийный номер не влезает в пять знаков")
    rows = []
    for i, (comp, el, x, y, z, name, seq) in enumerate(atoms, 1):
        rec = "HETATM" if comp in hetero_for else "ATOM  "
        nm = name if len(name) >= 4 else " " + name.ljust(3)
        rows.append(f"{rec}{i:5d} {nm[:4]}{'':1}{comp:>3} A{int(seq) % 10000:4d}    "
                    f"{x:8.3f}{y:8.3f}{z:8.3f}{1.0:6.2f}{0.0:6.2f}          {el:>2}")
    rows.append("END")
    return "\n".join(rows) + "\n"


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--cifdir", required=True, help="каталог с *.cif предсказанных комплексов")
    ap.add_argument("--eval", default="/tmp/ost_check",
                    help="каталог, содержащий пакет evaluation/ организаторов")
    ap.add_argument("--out", default="results/logs/k110_posegate.json")
    ap.add_argument("--work", default="/tmp/k110_work")
    a = ap.parse_args(argv)

    sys.path.insert(0, a.eval)
    try:
        from posebusters import PoseBusters
        from rdkit import Chem, RDLogger
        from evaluation.posebusters_utils import assign_bond_orders_from_template
    except ImportError as e:
        print(f"нет зависимостей шлюза: {e}\n"
              f"см. докстринг: отдельное окружение, posebusters + loguru, пакет evaluation/")
        return 2
    RDLogger.DisableLog("rdApp.*")

    import pandas as pd
    work = pathlib.Path(a.work)
    work.mkdir(parents=True, exist_ok=True)
    lig_csv = work / "ligands.csv"
    if not lig_csv.exists():
        urllib.request.urlretrieve(LIGCSV, lig_csv)
    smi = {str(r.Molecule_Name).upper(): str(r.SMILES)
           for r in pd.read_csv(lig_csv).itertuples()}
    print(f"лигандов в списке трека: {len(smi)}")

    cifs = sorted(pathlib.Path(a.cifdir).glob("*.cif"))
    if not cifs:
        print(f"в {a.cifdir} нет ни одного .cif")
        return 2
    print(f"комплексов на входе: {len(cifs)}\n")

    pb = PoseBusters(config="dock")
    rows, planted = [], None

    for cif in cifs:
        ident = None
        for key in smi:
            if key.split("-")[-1] in cif.name:
                ident = key
                break
        if ident is None:
            print(f"  {cif.name}: идентификатор не сопоставлен со списком трека, пропуск")
            continue
        atoms = parse_cif(cif)
        lig = [t for t in atoms if t[0].startswith("LIG")]
        rest = [t for t in atoms if not t[0].startswith("LIG")]
        if not lig:
            print(f"  {ident}: в файле НЕТ лиганда")
            rows.append({"id": ident, "ok": False, "why": "лиганд не найден"})
            continue

        # The ligand residue must be renamed LIG, which the submission format requires anyway.
        lig = [("LIG",) + t[1:] for t in lig]
        for tag, pose in (("", lig), ("ПОДЛОЖЕННАЯ", None)):
            if pose is None:
                # CONTROL: translate this ligand onto the heme iron. Must be flagged.
                fe = [t for t in rest if t[0] == "HEM" and t[1] == "FE"]
                if not fe:
                    continue
                cx, cy, cz = fe[0][2:5]
                lx = sum(t[2] for t in lig) / len(lig)
                ly = sum(t[3] for t in lig) / len(lig)
                lz = sum(t[4] for t in lig) / len(lig)
                pose = [(c, e, x - lx + cx, y - ly + cy, z - lz + cz, n, s)
                        for (c, e, x, y, z, n, s) in lig]
                if planted is not None:
                    continue

            prot_pdb = work / f"{ident}_prot.pdb"
            prot_pdb.write_text(pdb_block(rest, {"HEM"}))
            lig_pdb = work / f"{ident}{'_planted' if tag else ''}_lig.pdb"
            lig_pdb.write_text(pdb_block(pose, {"LIG"}))

            mol = Chem.MolFromPDBFile(str(lig_pdb), removeHs=False, sanitize=False)
            if mol is None:
                rows.append({"id": ident, "ok": False, "why": "RDKit не прочитал лиганд"})
                continue
            ref = Chem.MolFromSmiles(smi[ident])
            heavy_ref = ref.GetNumAtoms()
            heavy_got = sum(1 for at in mol.GetAtoms() if at.GetSymbol() != "H")
            if heavy_got != heavy_ref:
                rows.append({"id": ident, "ok": False,
                             "why": f"тяжёлых атомов {heavy_got}, SMILES требует {heavy_ref}"})
                continue
            try:
                mol = assign_bond_orders_from_template(ref, mol)
            except Exception as e:
                rows.append({"id": ident, "ok": False, "why": f"порядки связей: {e}"})
                continue

            df = pb.bust([mol], None, str(prot_pdb))
            cols = [c for c in df.columns if df[c].dtype == bool]
            failed = [c for c in cols if not bool(df.iloc[0][c])]
            rec = {"id": ident, "ok": True, "checks": len(cols),
                   "failed_n": len(failed), "failed": failed}
            if tag:
                planted = rec
                print(f"  {'ПОДЛОЖЕННАЯ ' + ident:<34} провалено {len(failed)}/{len(cols)}"
                      f"  {failed[:4]}")
            else:
                rows.append(rec)
                print(f"  {ident:<34} провалено {len(failed)}/{len(cols)}"
                      f"  {failed[:4] if failed else ''}")

    good = [r for r in rows if r.get("ok")]
    print(f"\nразобрано поз: {len(good)} из {len(cifs)}")
    if good:
        dist = {}
        for r in good:
            dist[r["failed_n"]] = dist.get(r["failed_n"], 0) + 1
        print("распределение числа проваленных проверок:")
        for k in sorted(dist):
            print(f"    провалено {k}: {dist[k]} поз")
        worst = sorted(good, key=lambda r: -r["failed_n"])[:3]
        for r in worst:
            if r["failed_n"]:
                print(f"    худшая: {r['id']} -> {r['failed']}")

    print("\nКОНТРОЛЬ --- подложенная поза на железе гема")
    if planted is None:
        print("    НЕ ПОСТРОЕНА (нет FE в гемe?) --- результат шлюза НЕДЕЙСТВИТЕЛЕН")
        code = 1
    elif planted["failed_n"] == 0:
        print("    НЕ ПОЙМАНА: 0 проваленных проверок --- шлюз слеп, результат НЕДЕЙСТВИТЕЛЕН")
        code = 1
    else:
        print(f"    поймана: провалено {planted['failed_n']} проверок -> шлюз видит плохую позу")
        code = 0

    pathlib.Path("results/logs").mkdir(parents=True, exist_ok=True)
    json.dump({"poses": rows, "planted": planted,
               "note": "порог POSEBUSTERS_MAX_FAILURES организаторами не опубликован; "
                       "здесь распределение, а не вердикт"},
              open(a.out, "w"), ensure_ascii=False, indent=1)
    print(f"\nсохранено: {a.out}")
    return code


if __name__ == "__main__":
    raise SystemExit(main())
