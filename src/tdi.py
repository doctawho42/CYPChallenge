import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import pandas as pd, numpy as np
from sklearn.metrics import matthews_corrcoef
t=pd.read_csv(D+"cyp-challenge-TRAIN_TDI.csv")
for c in ["CYP3A4","CYP2D6"]:
    d=t[[f"{c}_is_TDI",f"{c}_pIC50_direct_inhibition",f"{c}_pIC50_TDI_condition"]].dropna()
    di=d[f"{c}_pIC50_direct_inhibition"]; td=d[f"{c}_pIC50_TDI_condition"]; y=d[f"{c}_is_TDI"].astype(bool)
    rule=np.where(di>4, (td-di)>np.log10(2), td>4.301)
    print(f"{c}: rule reproduces label in {np.mean(rule==y)*100:.2f}% of {len(d)} rows; MCC={matthews_corrcoef(y,rule):.3f}")
    mis=d[rule!=y]
    print("   mismatches:",len(mis))
    if len(mis): print(mis.head(5).round(3).to_string())
    # how much would a perfect direct-inhibition model + perfect TDI-arm model give?
    print(f"   base rate={y.mean():.3f};  label vs (direct pIC50>4.6): MCC={matthews_corrcoef(y,di>4.6):.3f}")
    print(f"   label vs (TDI-arm pIC50>4.8): MCC={matthews_corrcoef(y,td>4.8):.3f}")
    best=max(((matthews_corrcoef(y,td>th),th) for th in np.arange(4.0,6.0,0.05)))
    print(f"   best single-threshold on TDI-arm pIC50: MCC={best[0]:.3f} at {best[1]:.2f}")
    best=max(((matthews_corrcoef(y,(td-di)>th),th) for th in np.arange(0.0,1.2,0.02)))
    print(f"   best single-threshold on shift: MCC={best[0]:.3f} at {best[1]:.2f}")
