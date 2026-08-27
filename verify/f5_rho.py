"""Тот же парный бутстрап, но для Спирмена -- документы утверждают, что по ранжированию
механистический блок «заметен», особенно на 2D6."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import numpy as np, pandas as pd, json, sys
from scipy.stats import spearmanr
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
rows=pd.read_csv(D+"rows.csv")
tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv").set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index()
oof=json.load(open(RES+"preds/oof.json"))
rng=np.random.default_rng(5); B=3000
print(f"{'фермент':8s} {'FP+DESC':>8s} {'+MECH':>8s} {'Δrho':>8s}  {'95% ДИ':>18s}  {'P(Δ<=0)':>8s}")
macro_d=np.zeros(B)
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; m=tr[col].notna().to_numpy()
    y=tr.loc[m,col].to_numpy()
    pa=np.asarray(oof[f"FP+DESC|{c}"]); pb=np.asarray(oof[f"FP+DESC+MECH|{c}"])
    n=len(y); d=np.empty(B)
    for b in range(B):
        i=rng.integers(0,n,n)
        d[b]=spearmanr(pb[i],y[i]).statistic-spearmanr(pa[i],y[i]).statistic
    r0=spearmanr(pa,y).statistic; r1=spearmanr(pb,y).statistic
    print(f"{c:8s} {r0:8.3f} {r1:8.3f} {r1-r0:+8.3f}  [{np.percentile(d,2.5):+7.3f},{np.percentile(d,97.5):+7.3f}]  {np.mean(d<=0):8.3f}")
