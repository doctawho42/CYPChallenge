"""Есть ли у решающего слоя потолок вообще? Оракульный вариант: истинные границы полосы."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd, json, sys
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
rows=pd.read_csv(D+"rows.csv")
tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv").set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index()
oof=json.load(open(RES+"preds/oof.json"))
print("=== асимметрия полосы: a = y-lo (вниз), b = hi-y (вверх) ===")
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; d=tr[tr[col].notna()]
    y=d[col].to_numpy(); a=y-d[col+"_conf_low"].to_numpy(); b=d[col+"_conf_high"].to_numpy()-y
    for lab,msk in [("активные y>=5",y>=5),("средние 4-5",(y>=4)&(y<5)),("неактивные y<4",y<4)]:
        if msk.sum()<10: continue
        print(f"  {c} {lab:16s} n={msk.sum():4d}  a={np.median(a[msk]):.2f}  b={np.median(b[msk]):.2f}  b-a={np.median(b[msk]-a[msk]):+.2f}")
print("\n=== оракул: сдвиг с ИСТИННЫМИ границами полосы ===")
rng=np.random.default_rng(0); S=600; res=[]
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; m=tr[col].notna().to_numpy()
    y=tr.loc[m,col].to_numpy(); lo=tr.loc[m,col+"_conf_low"].to_numpy(); hi=tr.loc[m,col+"_conf_high"].to_numpy()
    yhat=np.array(oof[f"FP+DESC+MECH|{c}"])
    a=y-lo; b=hi-y; sig=np.std(y-yhat)
    eps=rng.standard_normal((S,len(y)))
    ys=yhat[None,:]+sig*eps
    pooled=np.concatenate([ys-a[None,:],ys+b[None,:]],axis=0)
    y_or=np.median(pooled,axis=0)
    # и «идеальная проекция»: если прогноз вне полосы, притянуть к ближайшей границе
    y_proj=np.clip(yhat,lo,hi)
    f_=lambda p: round(strae(y,p,y_true_upper=hi,y_true_lower=lo),4)
    res.append(dict(cyp=c,raw=f_(yhat),oracle_shift=f_(y_or),oracle_clip=f_(y_proj),
                    mean_move=round(float(np.mean(y_or-yhat)),3)))
r=pd.DataFrame(res); r.loc[len(r)]=dict(cyp="MACRO",**{k:round(r[k].mean(),4) for k in r.columns if k!="cyp"})
print(r.to_string(index=False))
print("\noracle_clip -- недостижимая верхняя граница (знает истинную полосу целиком),")
print("нужна только чтобы понять, сколько ошибки вообще лежит вне полос.")
