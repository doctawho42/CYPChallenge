"""Гипотеза: plug-in порог проигрывает потому, что вероятности не откалиброваны.
Проверяем: калибруем изотонической регрессией во вложенной схеме, потом снова plug-in.

ДЕФЕКТ, НАЙДЕННЫЙ И ИСПРАВЛЕННЫЙ (пункт 235). Калибровочные группы здесь рисовались как
`rng.integers(0, 5, n)` --- случайно, а не по Бутине. Близкие аналоги оцениваемых строк
попадали в подгоночную половину, и измеренная выгода калибровки оказалась завышена: §10
цитировал +0.028 MCC на CYP3A4, а на каноническом сплите с полным вложением та же величина
равна +0.0235, причём каждый посидовый интервал накрывает ноль.

Здесь группы теперь берутся из `butina_folds` --- те же фолды, на которых получены сами
вневыборочные вероятности (`src/tdiprob.py` кластеризует подмножество TDI, `rows.SMILES[keep]`,
и это правильно: другой набор молекул --- другие фолды).

Числа этого файла после исправления НЕ совпадают с тем, что печаталось раньше. Развёрнутое
четырёхсидовое измерение живёт в `verify/k68_tdicalib.py`; здесь остался быстрый однопроходный
вариант, и цитировать следует k68."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import numpy as np, json
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import matthews_corrcoef, brier_score_loss
from cypsplit import butina_folds
P=json.load(open(RES+"preds/tdi_probs.json"))
_rows=__import__("pandas").read_csv(D+"rows.csv")
_tdi=__import__("pandas").read_csv(D+"cyp-challenge-TRAIN_TDI.csv").set_index("Molecule_Name")
_keep=_rows.Molecule_Name.isin(_tdi.index).to_numpy()
_T=_tdi.loc[_rows.Molecule_Name[_keep]].reset_index()
_FOLD,_=butina_folds(list(_rows.SMILES[_keep]))
def mcc_c(tp,tn,fp,fn):
    d=np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)); return 0.0 if d<=0 else (tp*tn-fp*fn)/d
def plugin_thr(p, ts):
    v=[mcc_c(p[p>t].sum(),(1-p[p<=t]).sum(),(1-p[p>t]).sum(),p[p<=t].sum()) for t in ts]
    return ts[int(np.argmax(v))]
ts=np.unique(np.round(np.linspace(0.02,0.95,200),4))
for key in ["CYP3A4","CYP2D6"]:
    p=np.asarray(P[key]["p"],float); y=np.asarray(P[key]["y"],int); n=len(y)
    # Фолды Бутины, а не случайные группы: иначе близкие аналоги оцениваемых строк
    # оказываются в подгоночной половине и выгода калибровки завышается (пункт 235).
    g=_FOLD[_T[f"{key}_is_TDI"].notna().to_numpy()]
    assert len(g)==n, f"{key}: фолдов {len(g)}, вероятностей {n} -- маска разошлась"
    pc_iso=np.zeros(n); pc_pla=np.zeros(n)
    for f in range(5):
        a,b=g!=f,g==f
        pc_iso[b]=IsotonicRegression(out_of_bounds="clip",y_min=0,y_max=1).fit(p[a],y[a]).predict(p[b])
        lr=LogisticRegression(C=1e6).fit(np.log(np.clip(p[a],1e-6,1-1e-6)/(1-np.clip(p[a],1e-6,1-1e-6))).reshape(-1,1),y[a])
        pc_pla[b]=lr.predict_proba(np.log(np.clip(p[b],1e-6,1-1e-6)/(1-np.clip(p[b],1e-6,1-1e-6))).reshape(-1,1))[:,1]
    print(f"\n=== {key}: n={n}, доля положительных {y.mean():.3f} ===")
    for tag,pp in [("сырые",p),("изотоническая",pc_iso),("Платт",pc_pla)]:
        emp=np.array([matthews_corrcoef(y,(pp>t).astype(int)) for t in ts])
        tp_=plugin_thr(pp,ts); mcc_plug=emp[np.argmin(abs(ts-tp_))]
        print(f"  {tag:14s} E[p]={pp.mean():.3f} Брайер={brier_score_loss(y,pp):.4f} | "
              f"plug-in порог {tp_:.3f} -> MCC {mcc_plug:.4f} | оптимум {emp.max():.4f} | "
              f"потеря {emp.max()-mcc_plug:+.4f}")
