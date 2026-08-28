"""Ablation: does an explicit CYP2D6 mechanistic block beat generic descriptors?

Same Butina-cluster folds, same learner, only the feature matrix changes.
Saves out-of-fold predictions so the paired bootstrap can be run afterwards.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd, sys, time, json
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from sklearn.ensemble import HistGradientBoostingRegressor
from cypsplit import butina_folds
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]

z = np.load(D + "feats.npz")
FP, DESC, MECH = z["FP"], z["DESC"], z["MECH"]
rows = pd.read_csv(D + "rows.csv")
tr = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
tr = tr.set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index()
assert len(tr) == len(FP)

fold, n_clusters = butina_folds(list(rows.SMILES))
print("clusters:", n_clusters, "| fold sizes:", np.bincount(fold), flush=True)

SETS = {
    "FP":            FP,
    "FP+DESC":       np.hstack([FP, DESC]),
    "FP+DESC+MECH":  np.hstack([FP, DESC, MECH]),
    "DESC+MECH":     np.hstack([DESC, MECH]),
    "MECH":          MECH,
    "DESC":          DESC,
}
out = {}
for name, X in SETS.items():
    t0 = time.time()
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        msk = tr[col].notna().to_numpy()
        y = tr.loc[msk, col].to_numpy(); Xi = X[msk]; fi = fold[msk]
        pred = np.zeros_like(y)
        for f in range(5):
            a, b = fi != f, fi == f
            if b.sum() == 0: continue
            m = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.06,
                                              max_leaf_nodes=31, l2_regularization=1.0,
                                              random_state=0)
            m.fit(Xi[a], y[a]); pred[b] = m.predict(Xi[b])
        out[f"{name}|{c}"] = pred.tolist()
    print(f"{name:16s} done in {time.time()-t0:.0f}s", flush=True)
json.dump(out, open(RES + "preds/oof.json", "w"))

# Provenance, so that four people on four machines can tell whose oof.json this is and
# whether it is still the one the document quotes. oof.json itself is a single 750 KB
# line and cannot be merged or eyeballed; this file can.
import subprocess, sklearn, platform
from cypsplit import fold_digest
def _git(*a):
    try: return subprocess.run(("git",)+a, capture_output=True, text=True).stdout.strip()
    except Exception: return "?"
# Dirtiness of the code that could have changed these predictions. results/preds/ is
# excluded because this script has just rewritten oof.json, so a plain `git status`
# would report dirty on every single run; docs/ is excluded because prose cannot move
# a number.
_dirty = [ln for ln in _git("status", "--porcelain").splitlines()
          if "results/preds/" not in ln and "docs/" not in ln]
json.dump({
    "split_digest": fold_digest(fold),
    "n_clusters": int(n_clusters),
    "sklearn": sklearn.__version__,
    "numpy": np.__version__,
    "python": platform.python_version(),
    "git_commit": _git("rev-parse", "--short", "HEAD"),
    "code_dirty": bool(_dirty),
}, open(RES + "preds/oof.meta.json", "w"), indent=1)
print("saved")
