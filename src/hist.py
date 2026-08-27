import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import numpy as np, pandas as pd

tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
sc=pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
piv=sc.pivot_table(index="Molecule_Name",columns="enzyme",values="log2fc_estimate")
tr=tr.join(piv,on="Molecule_Name")
def bars(v,lo,hi,nb):
    h,e=np.histogram(v.dropna(),bins=nb,range=(lo,hi)); return h/h.max(), e
for name,v,lo,hi in [("3A4_pIC50",tr.CYP3A4_pIC50_direct_inhibition,2,7.5),
                     ("2D6_pIC50",tr.CYP2D6_pIC50_direct_inhibition,2,7.5)]:
    h,e=bars(v,lo,hi,22)
    print(name, "n=",int(v.notna().sum()), "iqr=",round(float(v.quantile(.75)-v.quantile(.25)),2))
    print("  ", " ".join(f"{x:.2f}" for x in h))
print("edges", " ".join(f"{x:.2f}" for x in e))
print()
for name,v,lo,hi in [("2D6_log2fc",tr.CYP2D6,-3.5,1.0)]:
    h,e=bars(v,lo,hi,22); print(name,"n=",int(v.notna().sum()),"iqr=",round(float(v.quantile(.75)-v.quantile(.25)),2))
    print("  "," ".join(f"{x:.2f}" for x in h)); print("edges"," ".join(f"{x:.2f}" for x in e))
# quartiles for shading
for c in ["CYP2D6","CYP3A4"]:
    v=tr[f"{c}_pIC50_direct_inhibition"]
    print(c,"q25",round(v.quantile(.25),2),"q75",round(v.quantile(.75),2),"min",round(v.min(),2),"max",round(v.max(),2))
