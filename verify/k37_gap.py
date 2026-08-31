"""Разрыв между тем, что берёт наша модель, и тем, что берёт одна скрининговая колонка.

Обе величины считаются здесь заново и одинаково: пятикратная кросс-валидация по тем же
кластерным фолдам, ST-RAE после аффинной пары. Скрининг подаётся как ОДНО число через
изотонику, подогнанную на обучающих фолдах, --- то есть без химии вообще.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / 'src'))
import json, numpy as np, pandas as pd
from cyppaths import D, RES, tutorial
tutorial()
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from sklearn.isotonic import IsotonicRegression
from cypsplit import butina_folds
from shrinkchoice import fit_apply
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
rows=pd.read_csv(D+"rows.csv")
tr=(pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
      .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
sc=pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
O=json.load(open(RES+"preds/oof.json"))
fold,_=butina_folds(list(rows.SMILES))

print(f"{'фермент':8s} {'модель+пара':>12s} {'скрининг один':>14s} {'разрыв':>8s} {'n со скрином':>13s}")
gaps={}
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; m=tr[col].notna().to_numpy()
    y=tr.loc[m,col].to_numpy(); lo=tr.loc[m,col+"_conf_low"].to_numpy(); hi=tr.loc[m,col+"_conf_high"].to_numpy()
    fi=fold[m]; u=np.ones(len(y))/len(y)
    p=np.asarray(O[f"FP+DESC+MECH|{c}"],float)
    mod=float(strae(y,fit_apply(p,lo,hi,fi,u),y_true_upper=hi,y_true_lower=lo))
    sub=sc[sc.enzyme==c].set_index("Molecule_Name")
    l2=rows.set_index("Molecule_Name").join(sub["log2fc_estimate"])["log2fc_estimate"].to_numpy(float)[m]
    ok=~np.isnan(l2)
    q=np.full(len(y),np.nan)
    for f in range(5):
        te=(fi==f)&ok; trn=(fi!=f)&ok
        if te.sum()==0 or trn.sum()<50: continue
        q[te]=IsotonicRegression(increasing=False,out_of_bounds="clip").fit(l2[trn],y[trn]).predict(l2[te])
    g=~np.isnan(q)
    scr=float(strae(y[g],q[g],y_true_upper=hi[g],y_true_lower=lo[g]))
    gaps[c]=(mod,scr)
    print(f"{c:8s} {mod:12.4f} {scr:14.4f} {mod-scr:8.3f} {int(g.sum()):13d}")
print(f"{'МАКРО':8s} {np.mean([v[0] for v in gaps.values()]):12.4f} "
      f"{np.mean([v[1] for v in gaps.values()]):14.4f} "
      f"{np.mean([v[0]-v[1] for v in gaps.values()]):8.3f}")
