"""Build a TabICL input bundle from any set of frozen embeddings, with the order control each one
needs. Supersedes `src/tabicl_in.py` and `src/tabicl_test_in.py`, which were one-offs for the first
two bundles and are kept because their logs cite them.

    uv run python src/tabicl_bundle.py --out data/tabicl_in3.npz \
        --reps C=emb_chemprop_chemeleon.npz Mi=emb_chemprop_minimol.npz \
               U=emb_chemprop_medium.npz+emb_chemprop_chemeleon.npz+emb_chemprop_surrogate_adme.npz

A `+` in a spec concatenates the files' train arrays column-wise, which is the arm item 117 did NOT
measure: that item concatenated ONE embedding to our own FP+DESC+MECH matrix under the ablation's
head and lost 0.0075 of rank. Concatenating several EMBEDDINGS and capping the result by PCA for an
in-context learner is a different object, and it is place 41's published recipe.

THE ORDER CONTROL, per file, and both branches must succeed:

  * a file WITH `molecule_name` is checked element by element against `data/rows.csv`, and the same
    comparison on the names rolled by one position must FAIL -- a comparison that cannot fail is
    not a control.
  * a file WITHOUT names gets `k121.align_control`: a RidgeCV probe on the true order must beat the
    same probe on the matrix rolled by one row by more than 0.05 of Spearman. Rolling leaves the
    matrix, its scale and its rank untouched and breaks only the correspondence to the labels.

A concatenation is controlled twice: every part is checked by its own branch first, and then
`align_control` runs on the assembled matrix, because a column-stacking bug would survive the parts
passing individually.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "verify"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import os

import numpy as np
import pandas as pd

from cypsplit import butina_folds, fold_digest
from submit import CYPS
from k121_chemeleon import align_control

GOLDEN = "2d93c19815e14261"


def load_one(fn, rows, Y, M, fold0):
    z = np.load(D + fn, allow_pickle=True)
    E = z["train"].astype(np.float64)
    if E.shape[0] != len(rows):
        raise SystemExit("%s: %d строк против %d в rows.csv" % (fn, E.shape[0], len(rows)))
    if not np.isfinite(E).all():
        raise SystemExit("%s: есть не-конечные значения" % fn)
    if "molecule_name" in z.files:
        nm = np.asarray(z["molecule_name"], object)
        r = rows.Molecule_Name.to_numpy(object)
        ok, rolled = bool((nm == r).all()), bool((np.roll(nm, 1) == r).all())
        if not ok or rolled:
            raise SystemExit("%s: порядок не совпал (совпал=%s, сдвинутый=%s)" % (fn, ok, rolled))
        print("    %-44s %s  имена сверены, сдвинутый контроль False" % (fn, E.shape), flush=True)
    else:
        if not align_control(E, Y, M, fold0, e=0):
            raise SystemExit("%s: порядок строк не подтверждён метками" % fn)
        print("    %-44s %s  имён нет, контроль метками пройден" % (fn, E.shape), flush=True)
    return E


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--reps", nargs="+", required=True,
                    help="KEY=file.npz или KEY=a.npz+b.npz для конкатенации")
    ap.add_argument("--seeds", default="0,1,2,3")
    a = ap.parse_args(argv)

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr["%s_pIC50_direct_inhibition" % c].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)
    fold0, _ = butina_folds(list(rows.SMILES))
    if fold_digest(fold0) != GOLDEN:
        raise SystemExit("разбиение сдвинулось: %s" % fold_digest(fold0))

    out = {"Y": Y, "cyps": np.array(CYPS)}
    for spec in a.reps:
        key, _, files = spec.partition("=")
        parts = files.split("+")
        print("  %s <- %s" % (key, " + ".join(parts)), flush=True)
        mats = [load_one(f, rows, Y, M, fold0) for f in parts]
        E = mats[0] if len(mats) == 1 else np.hstack(mats)
        if len(mats) > 1:
            assert E.shape[1] == sum(m.shape[1] for m in mats), E.shape
            print("    собрано %s из %s" % (E.shape, "+".join(str(m.shape[1]) for m in mats)),
                  flush=True)
            if not align_control(E, Y, M, fold0, e=0):
                raise SystemExit("%s: собранная матрица не прошла контроль метками" % key)
            print("    КОНТРОЛЬ сборки: пройден", flush=True)
        out[key] = E.astype(np.float32)

    digs = []
    for s in [int(x) for x in a.seeds.split(",") if x.strip()]:
        f, _ = butina_folds(list(rows.SMILES), seed=s)
        out["fold%d" % s] = f
        digs.append(fold_digest(f))
    out["digests"] = np.array(digs)
    np.savez_compressed(a.out, **out)
    print("\n  записано %s, %.1f МБ; дайджесты %s"
          % (a.out, os.path.getsize(a.out) / 1e6, ",".join(d[:8] for d in digs)))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
