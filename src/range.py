"""Is 2D6 hard because the model can't learn it, or because the DRC labels are range-restricted?

Same features, same folds. Compare how well the model predicts the FULL-RANGE readout
(single-concentration log2fc, covers actives and inactives alike) versus the RANGE-RESTRICTED
one (DRC pIC50, only for compounds that were promoted because they were hits).
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from scipy.stats import spearmanr
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
RDLogger.DisableLog('rdApp.*')
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
z = np.load(D + "feats.npz"); X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
rows = pd.read_csv(D + "rows.csv")
tr = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv").set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index()
sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
piv = sc.pivot_table(index="Molecule_Name", columns="enzyme", values="log2fc_estimate")
tr = tr.join(piv, on="Molecule_Name")
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
bits = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]
dists = []
for i in range(1, len(bits)): dists.extend([1 - x for x in DataStructs.BulkTanimotoSimilarity(bits[i], bits[:i])])
cl = Butina.ClusterData(dists, len(bits), 0.35, isDistData=True)
cid = np.zeros(len(bits), int)
for k, c in enumerate(cl):
    for i in c: cid[i] = k
fold = np.random.default_rng(0).integers(0, 5, len(cl))[cid]
res = []
for c in CYPS:
    r = {"cyp": c}
    for tag, col in [("pIC50_DRC", f"{c}_pIC50_direct_inhibition"), ("log2fc_full", c)]:
        msk = tr[col].notna().to_numpy(); y = tr.loc[msk, col].to_numpy()
        Xi, fi = X[msk], fold[msk]; p = np.zeros_like(y)
        for f in range(5):
            a, b = fi != f, fi == f
            if b.sum() == 0: continue
            p[b] = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.06, random_state=0).fit(Xi[a], y[a]).predict(Xi[b])
        r[f"n_{tag}"] = int(msk.sum()); r[f"rho_{tag}"] = round(spearmanr(p, y).statistic, 3)
        r[f"R2_{tag}"] = round(1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum(), 3)
        r[f"iqr_{tag}"] = round(float(np.percentile(y, 75) - np.percentile(y, 25)), 2)
    res.append(r)
print(pd.DataFrame(res).to_string(index=False))
