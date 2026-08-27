import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import pandas as pd, numpy as np

tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv"); sc=pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
print("=== per-enzyme log2fc distribution ===")
print(sc.groupby("enzyme").log2fc_estimate.describe(percentiles=[.01,.05,.25,.5,.75,.95,.99]).round(2).to_string())
piv=sc.pivot_table(index="Molecule_Name",columns="enzyme",values="log2fc_estimate")
m=tr.set_index("Molecule_Name").join(piv,how="inner")
print("\n=== Hill check: median log2fc in pIC50 bins (C = 49.5 uM -> pC = 4.31) ===")
bins=[0,3,3.5,4,4.31,4.6,5,5.5,6,6.5,10]
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"
    d=m[[col,c]].dropna().copy(); d["b"]=pd.cut(d[col],bins)
    g=d.groupby("b",observed=True)[c].agg(['median','count']).round(2)
    print(f"\n{c}:"); print(g.T.to_string())
