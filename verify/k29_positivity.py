"""Правильное предусловие слоя отбора: ПЕРЕКРЫТИЕ, а не дисперсия весов.

Перевзвешивание представляет целевую популяцию только там, где у неё есть носитель в
отобранной. Если P(отбор | x) обращается в ноль на целой области, эту область не вытянет
никакой вес --- это нарушение позитивности, и ESS его не видит: он меряет разброс весов
СРЕДИ отобранных, а не покрытие неотобранных.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / 'src'))
import numpy as np, pandas as pd
from cyppaths import D
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
rows=pd.read_csv(D+"rows.csv")
tr=(pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
      .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
sc=pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
piv=sc.pivot_table(index="Molecule_Name",columns="enzyme",values="log2fc_estimate")
S=rows.set_index("Molecule_Name").join(piv)[CYPS].to_numpy(float)

print(f"{'фермент':8s} {'отобрано':>9s} {'не отобр.':>10s} {'мин P среди отобр.':>19s} "
      f"{'доля неотобр. ниже него':>24s}")
for e,c in enumerate(CYPS):
    lab=tr[f"{c}_pIC50_direct_inhibition"].notna().to_numpy()
    s=S[:,e]; ok=~np.isnan(s)
    x=s[ok].reshape(-1,1); yb=lab[ok].astype(int)
    pr=cross_val_predict(LogisticRegression(max_iter=1000),x,yb,cv=5,method="predict_proba")[:,1]
    ps, pu = pr[yb==1], pr[yb==0]
    q=np.quantile(ps,0.01)          # 1-й процентиль склонности среди отобранных
    print(f"{c:8s} {len(ps):9d} {len(pu):10d} {q:19.4f} {100*(pu<q).mean():23.1f} %")

print(f"\n{'фермент':8s} {'скрин у отобранных':>20s} {'скрин у неотобранных':>22s} {'перекрытие':>12s}")
for e,c in enumerate(CYPS):
    lab=tr[f"{c}_pIC50_direct_inhibition"].notna().to_numpy()
    s=S[:,e]; ok=~np.isnan(s)
    a,b=s[ok&lab], s[ok&~lab]
    lo,hi=np.quantile(a,[0.01,0.99])
    print(f"{c:8s} {np.median(a):9.3f} [{np.quantile(a,0.05):.2f},{np.quantile(a,0.95):.2f}]"
          f" {np.median(b):11.3f} [{np.quantile(b,0.05):.2f},{np.quantile(b,0.95):.2f}]"
          f" {100*((b>=lo)&(b<=hi)).mean():11.1f} %")
print("""
Как читать. Второй столбец справа в первой таблице --- доля целевой популяции, у которой
склонность ниже, чем у 99% отобранных. Это доля, которую перевзвешивание не достаёт в
принципе. Ниже 10% --- слой строится; выше половины --- позитивность нарушена, и никакой
вес не спасёт.""")
