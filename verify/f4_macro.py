"""Парный бутстрап на МАКРО-уровне -- именно эта величина сортирует лидерборд.
Ресэмплим СОЕДИНЕНИЯ (строки общей таблицы), потом для каждого фермента берём те из них,
у кого есть метка. Так сохраняется корреляция между ферментами."""
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
N=len(tr)
# для каждого фермента: индексы строк в общей таблице и позиция внутри вектора предсказаний
IDX={}; TRU={}
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; m=tr[col].notna().to_numpy()
    idx=np.where(m)[0]
    pos=np.full(N,-1); pos[idx]=np.arange(len(idx))
    IDX[c]=(idx,pos)
    TRU[c]=(tr.loc[m,col].to_numpy(), tr.loc[m,col+"_conf_low"].to_numpy(), tr.loc[m,col+"_conf_high"].to_numpy())

def macro(setname, take):
    """take -- массив индексов строк общей таблицы (с повторами)."""
    v=[]
    for c in CYPS:
        idx,pos=IDX[c]; y,lo,hi=TRU[c]; p=np.asarray(oof[f"{setname}|{c}"])
        sel=pos[take]; sel=sel[sel>=0]
        if len(sel)<20: return np.nan
        v.append(strae(y[sel],p[sel],y_true_upper=hi[sel],y_true_lower=lo[sel]))
    return float(np.mean(v))

SETS=["FP","DESC","MECH","DESC+MECH","FP+DESC","FP+DESC+MECH"]
print("точечные макро-оценки:")
for s in SETS: print(f"  {s:14s} {macro(s,np.arange(N)):.4f}")

rng=np.random.default_rng(11); B=3000
base="FP+DESC"
boot={s:np.empty(B) for s in SETS}
for b in range(B):
    take=rng.integers(0,N,N)
    for s in SETS: boot[s][b]=macro(s,take)
print(f"\nПарный бутстрап макро-ST-RAE относительно {base}, B={B}")
print(f"{'набор':16s} {'Δмакро':>8s}  {'95% ДИ':>20s}  {'P(хуже базы)':>12s}")
for s in SETS:
    if s==base: continue
    d=boot[s]-boot[base]
    print(f"{s:16s} {d.mean():+8.4f}  [{np.percentile(d,2.5):+7.4f},{np.percentile(d,97.5):+7.4f}]  {np.mean(d>0):12.3f}")
# и разброс самой макро-оценки
print(f"\nразброс макро-оценки FP+DESC+MECH по бутстрапу: sd={boot['FP+DESC+MECH'].std():.4f}, "
      f"95% [{np.percentile(boot['FP+DESC+MECH'],2.5):.4f},{np.percentile(boot['FP+DESC+MECH'],97.5):.4f}]")
np.save(RES + "macro_boot.npy", np.array([boot[s] for s in SETS]))
