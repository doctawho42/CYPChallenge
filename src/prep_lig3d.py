"""Build data/lig3d_{train,test}.sdf -- the 3D conformers the structural scripts consume.

The second half of a reproducibility hole. `src/dock.py`, `src/shape3d.py`, `src/overlay.py`
and `src/quantum.py` all read these two SDFs, `data/` keeps them out of git, and no committed
script wrote them: `git grep lig3d -- '*.py'` finds exactly one hit, and it is a read. So the
inputs four scripts depend on existed on one machine with no recipe. `src/prep_receptors.py`
closed the receptor half; this closes the ligand half.

The parameters are not invented here. They are already pinned twice in the repository --
`src/shape3d.py` and `src/quantum.py` both use ETKDGv3 with `randomSeed = 0xC0FFEE` and
`MMFFOptimizeMolecule(maxIters=400)` on a hydrogen-added molecule, one conformer -- and the
files on disk match that (tr_0 carries 21 atoms of which 9 are hydrogens, one conformer).
This script only assembles that recipe into the step that writes the files.

**The naming is the trap.** `_Name` is the POSITIONAL ROW INDEX, not a running counter:
`tr_<i>` where i indexes `data/rows.csv`, and `te_<i>` where i indexes the blinded test CSV.
Three training molecules fail to embed, so the train file holds 4902 records whose names run
to tr_4904 with three gaps. Renumbering them would silently desynchronise everything that
maps a pose back to a label -- the same class of failure `data/rows.csv` exists to prevent.
`data/lig3d_index.npz` records which indices succeeded (`ok_train` 4902, `bad_train` 3,
`ok_test` 750, `bad_test` 0).

Default-safe: `--dry-run` rebuilds and COMPARES against what is on disk without writing, and
`--limit N` checks the first N molecules so the comparison costs seconds rather than an hour.
That order -- verify, then write -- is not a style preference. Rewriting the receptors before
running the check that was supposed to justify it destroyed the originals earlier today, in a
gitignored directory with no backup, and the retraction had to be argued from line counts.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D

import argparse
import os

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")

SEED = 0xC0FFEE
MMFF_ITERS = 400


def embed(smiles):
    """One ETKDGv3 conformer, hydrogens kept, short MMFF cleanup. None if it fails."""
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    mh = Chem.AddHs(m)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = SEED
    if AllChem.EmbedMolecule(mh, ps) != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mh, maxIters=MMFF_ITERS)
    except Exception:
        pass                      # MMFF has no parameters for some atoms; keep the ETKDG pose
    return mh


def build(smiles_list, tag, limit=None):
    """Returns (molecules, ok_indices, bad_indices). Names carry the positional index."""
    mols, ok, bad = [], [], []
    n = len(smiles_list) if limit is None else min(limit, len(smiles_list))
    for i in range(n):
        mh = embed(smiles_list[i])
        if mh is None:
            bad.append(i)
            continue
        mh.SetProp("_Name", "%s_%d" % (tag, i))
        mols.append(mh)
        ok.append(i)
    return mols, ok, bad


def read_existing(path, limit=None):
    """Name -> atom count for the molecules already on disk."""
    if not os.path.exists(path):
        return None
    out = {}
    supp = Chem.ForwardSDMolSupplier(path, removeHs=False)
    for i, m in enumerate(supp):
        if limit is not None and i >= limit:
            break
        if m is None:
            continue
        out[m.GetProp("_Name") if m.HasProp("_Name") else "?%d" % i] = m.GetNumAtoms()
    return out


def compare(built, existing, tag):
    """Report whether the rebuild reproduces the file, per name and per atom count."""
    if existing is None:
        print("    %-6s на диске файла нет -- сравнивать не с чем" % tag)
        return None
    got = {m.GetProp("_Name"): m.GetNumAtoms() for m in built}
    shared = set(got) & set(existing)
    only_new = sorted(set(got) - set(existing))[:5]
    only_old = sorted(set(existing) - set(got))[:5]
    diff = [k for k in shared if got[k] != existing[k]]
    print("    %-6s имён совпало %d, только в сборке %s, только на диске %s"
          % (tag, len(shared), only_new or "нет", only_old or "нет"))
    print("           расхождений по числу атомов: %d%s"
          % (len(diff), (" -- " + ", ".join(diff[:5])) if diff else ""))
    return len(shared) > 0 and not diff and not only_new


def write_sdf(path, mols, dry):
    if dry:
        return "не записан (сухой прогон)"
    tmp = path + ".tmp"
    w = Chem.SDWriter(tmp)
    for m in mols:
        w.write(m)
    w.close()
    os.replace(tmp, path)
    return "записан (%d записей)" % len(mols)


def main():
    ap = argparse.ArgumentParser(description="Сборка 3D-конформеров для структурных скриптов")
    ap.add_argument("--dry-run", action="store_true",
                    help="пересобрать и СРАВНИТЬ с тем, что на диске, ничего не записывая")
    ap.add_argument("--limit", type=int, default=None,
                    help="только первые N молекул каждой половины (для быстрой сверки)")
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    print("Сборка data/lig3d_{train,test}.sdf")
    print("ETKDGv3, сид 0x%X, MMFF %d итераций, водороды сохраняются -- параметры те же,"
          % (SEED, MMFF_ITERS))
    print("что в src/shape3d.py и src/quantum.py, а не выбранные здесь заново.")
    if a.dry_run:
        print("РЕЖИМ: сухой прогон, файлы не записываются%s"
              % ("" if a.limit is None else ", только первые %d молекул" % a.limit))
    print()

    verdicts = {}
    for tag, smiles, path in (("tr", list(rows.SMILES), D + "lig3d_train.sdf"),
                              ("te", list(te.SMILES), D + "lig3d_test.sdf")):
        mols, ok, bad = build(smiles, tag, a.limit)
        print("    %-6s встроено %d, не удалось %d%s"
              % (tag, len(ok), len(bad), (" (индексы %s)" % bad) if bad else ""))
        verdicts[tag] = compare(mols, read_existing(path, a.limit), tag)
        print("           %s" % write_sdf(path, mols, a.dry_run))
        if not a.dry_run and a.limit is None:
            np.savez(D + "lig3d_index.npz",
                     **{"ok_" + ("train" if tag == "tr" else "test"): np.array(ok),
                        "bad_" + ("train" if tag == "tr" else "test"): np.array(bad)})

    print()
    if all(v is True for v in verdicts.values()):
        print("    ВЕРДИКТ: пересборка воспроизводит файлы на диске в проверенной части.")
    elif any(v is False for v in verdicts.values()):
        print("    ВЕРДИКТ: РАСХОЖДЕНИЕ. Ничего не записывать, пока не объяснено --")
        print("    эти SDF читают dock.py, shape3d.py, overlay.py и quantum.py.")
    else:
        print("    ВЕРДИКТ: сравнивать было не с чем.")


if __name__ == "__main__":
    main()
