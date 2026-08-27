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
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
from sklearn.ensemble import HistGradientBoostingRegressor
RDLogger.DisableLog('rdApp.*')
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]

z = np.load(D + "feats.npz")
FP, DESC, MECH = z["FP"], z["DESC"], z["MECH"]
rows = pd.read_csv(D + "rows.csv")
tr = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
tr = tr.set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index()
assert len(tr) == len(FP)

gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
bits = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]
dists = []
for i in range(1, len(bits)):
    dists.extend([1 - x for x in DataStructs.BulkTanimotoSimilarity(bits[i], bits[:i])])
cl = Butina.ClusterData(dists, len(bits), 0.35, isDistData=True)
cid = np.zeros(len(bits), int)
for k, c in enumerate(cl):
    for i in c: cid[i] = k
fold = np.random.default_rng(0).integers(0, 5, len(cl))[cid]
print("clusters:", len(cl), "| fold sizes:", np.bincount(fold), flush=True)

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
print("saved")
