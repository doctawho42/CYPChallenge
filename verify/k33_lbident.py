"""Идентифицируется ли сдвиг по счёту лидерборда, и что этому мешает.

ST-RAE --- отношение. Числитель зависит от сдвига И от качества предсказаний; знаменатель ---
от разброса ТЕСТОВЫХ меток, которого у нас нет. Два способа ошибиться, и они меряются здесь
по отдельности: сначала во сколько сдвиг путается с качеством, потом --- со знаменателем.

Итог в пункте 130: качество не мешает (загрязнение 2-12 % нынешней неопределённости по δ),
мешает знаменатель, и на CYP2D6 он съедает всё. Пинится он через tau.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / 'src'))
import json, numpy as np, pandas as pd
from cyppaths import D, RES, tutorial
tutorial()
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
rows=pd.read_csv(D+"rows.csv")
tr=(pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
      .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
P=json.load(open(RES+"preds/oof_pool_all.json"))["preds"]
G=json.load(open(RES+"preds/oof_gp.json"))["preds"]
W=json.load(open(RES+"preds/oof_weak.json"))["preds"]
rng=np.random.default_rng(0)
print(f"{'фермент':8s} {'dScore/dСдвиг':>14s} {'сдвиг, эквивалентный':>22s} {'он же в единицах':>18s}")
print(f"{'':8s} {'при сдвиге 0':>14s} {'ухудшению ранга на 0.01':>22s} {'нашей неопр. по δ':>18s}")
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; m=tr[col].notna().to_numpy()
    y=tr.loc[m,col].to_numpy(); lo=tr.loc[m,col+"_conf_low"].to_numpy(); hi=tr.loc[m,col+"_conf_high"].to_numpy()
    p=np.mean([np.asarray(P[f"0|независимо|{c}"],float),np.asarray(P[f"0|пул|{c}"],float),
               np.asarray(G[f"0|GP|{c}"],float),np.asarray(W[f"0|гребневая|{c}"],float)],axis=0)
    f=lambda s,q=p: float(strae(y,q+s,y_true_upper=hi,y_true_lower=lo))
    h=0.05; slope=(f(h)-f(-h))/(2*h)
    # Порча качества: подмешать шум так, чтобы ранг упал на 0.01
    from scipy.stats import spearmanr
    r0=spearmanr(y,p).statistic; sd=p.std()
    lo_a,hi_a=0.0,2.0
    for _ in range(30):
        a=(lo_a+hi_a)/2
        pn=p+a*sd*rng.standard_normal(len(p))
        if spearmanr(y,pn).statistic > r0-0.01: lo_a=a
        else: hi_a=a
    pn=p+((lo_a+hi_a)/2)*sd*rng.standard_normal(len(p))
    dscore=f(0.0,pn)-f(0.0)
    equiv=dscore/slope if abs(slope)>1e-9 else float('nan')
    print(f"{c:8s} {slope:14.4f} {equiv:22.3f} {abs(equiv)/0.5:17.2f}x")
print("""
Как читать. Второй столбец --- на сколько меняется счёт при сдвиге на единицу. Третий ---
какой СДВИГ дал бы то же изменение счёта, что падение ранга всего на 0.01, то есть на полтора
шумовых пола. Если это доли pIC50, сдвиг и качество путаются, и одно число лидерборда их не
разделяет; наша нынешняя неопределённость по δ порядка 0.5, и четвёртый столбец говорит,
во сколько раз эквивалент больше или меньше неё.""")
import sys; sys.path.insert(0,'.'); sys.path.insert(0,'src')
import json, numpy as np, pandas as pd
from cyppaths import D, RES, tutorial
tutorial()
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
rows=pd.read_csv(D+"rows.csv")
tr=(pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
      .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
P=json.load(open(RES+"preds/oof_pool_all.json"))["preds"]
G=json.load(open(RES+"preds/oof_gp.json"))["preds"]
W=json.load(open(RES+"preds/oof_weak.json"))["preds"]
print(f"{'фермент':8s} {'счёт':>8s} {'наклон':>8s} {'сдвиг при ошибке':>18s} {'то же при 20%':>15s}")
print(f"{'':8s} {'':8s} {'':8s} {'знаменателя 10%':>18s} {'':15s}")
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; m=tr[col].notna().to_numpy()
    y=tr.loc[m,col].to_numpy(); lo=tr.loc[m,col+"_conf_low"].to_numpy(); hi=tr.loc[m,col+"_conf_high"].to_numpy()
    p=np.mean([np.asarray(P[f"0|независимо|{c}"],float),np.asarray(P[f"0|пул|{c}"],float),
               np.asarray(G[f"0|GP|{c}"],float),np.asarray(W[f"0|гребневая|{c}"],float)],axis=0)
    f=lambda s: float(strae(y,p+s,y_true_upper=hi,y_true_lower=lo))
    s0=f(0.0); h=0.05; slope=(f(h)-f(-h))/(2*h)
    print(f"{c:8s} {s0:8.4f} {slope:8.4f} {abs(0.10*s0/slope):18.3f} {abs(0.20*s0/slope):15.3f}")
# Насколько вообще может разъехаться знаменатель: разброс |y-mean| по фолдам
print(f"\nразброс знаменателя между нашими фолдами (косвенная оценка того, насколько")
print(f"он может отличаться на тесте того же размера):")
from cypsplit import butina_folds
fold,_=butina_folds(list(rows.SMILES))
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; m=tr[col].notna().to_numpy()
    y=tr.loc[m,col].to_numpy(); fi=fold[m]
    v=[np.mean(np.abs(y[fi==f]-y[fi==f].mean())) for f in range(5) if (fi==f).sum()>30]
    print(f"  {c}: {np.mean(v):.3f} +- {np.std(v):.3f}  ({100*np.std(v)/np.mean(v):.1f} %)")
print("""
Как читать. Столбцы справа --- на сколько сместится восстановленный сдвиг, если знаменатель
угадан с ошибкой 10 и 20 процентов. Сравнивать с нашей нынешней неопределённостью по δ около
0.5 и с загрязнением от качества 0.01-0.06.""")
