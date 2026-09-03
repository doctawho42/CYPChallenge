"""Срез по порогу прибора: доказывает или ломает прочтение пункта 134.

Прочтение: разрыв на 1A2 несут соединения, чья метка --- экстраполяция ниже порога
разрешения, где скрининг ИЗМЕРЯЕТ неактивность, а модель обязана предсказать число для
неразрешённой кривой. Тогда ограничение обеих колонок на y > pC0 обязано схлопнуть разрыв
на 1A2 и почти не тронуть 2D6.

Тест не тавтологичен, и вот почему: у 3A4 экстраполированных меток вдвое больше (51.4 %),
а разрыв наименьший (0.269). Значит сама доля неразрешённых меток разрыв не объясняет.
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
PC0=-np.log10(49.5e-6)
rows=pd.read_csv(D+"rows.csv")
tr=(pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
      .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
sc=pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
O=json.load(open(RES+"preds/oof.json")); fold,_=butina_folds(list(rows.SMILES))
print(f"порог прибора pC0 = {PC0:.3f}\n")
print(f"{'фермент':8s} {'доля ниже':>10s} | {'разрыв весь':>12s} {'разрыв y>pC0':>13s} {'схлопнулся на':>14s}")
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; m=tr[col].notna().to_numpy()
    y=tr.loc[m,col].to_numpy(); lo=tr.loc[m,col+"_conf_low"].to_numpy(); hi=tr.loc[m,col+"_conf_high"].to_numpy()
    fi=fold[m]; u=np.ones(len(y))/len(y)
    p=fit_apply(np.asarray(O[f"FP+DESC+MECH|{c}"],float),lo,hi,fi,u)
    sub=sc[sc.enzyme==c].set_index("Molecule_Name")
    l2=rows.set_index("Molecule_Name").join(sub["log2fc_estimate"])["log2fc_estimate"].to_numpy(float)[m]
    ok=~np.isnan(l2); q=np.full(len(y),np.nan)
    for f in range(5):
        te=(fi==f)&ok; trn=(fi!=f)&ok
        if te.sum()==0 or trn.sum()<50: continue
        q[te]=IsotonicRegression(increasing=False,out_of_bounds="clip").fit(l2[trn],y[trn]).predict(l2[te])
    g=~np.isnan(q)
    def gap(sel):
        a=float(strae(y[sel],p[sel],y_true_upper=hi[sel],y_true_lower=lo[sel]))
        b=float(strae(y[sel],q[sel],y_true_upper=hi[sel],y_true_lower=lo[sel]))
        return a-b
    full=gap(g); above=gap(g&(y>PC0))
    print(f"{c:8s} {100*(y<PC0).mean():9.1f} % | {full:12.3f} {above:13.3f} {100*(1-above/full):13.0f} %")
