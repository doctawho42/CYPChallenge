"""Гипотеза: plug-in порог проигрывает потому, что вероятности не откалиброваны.
Проверяем: калибруем изотонической регрессией во вложенной схеме, потом снова plug-in."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, json
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import matthews_corrcoef, brier_score_loss
P=json.load(open(RES+"preds/tdi_probs.json"))
def mcc_c(tp,tn,fp,fn):
    d=np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)); return 0.0 if d<=0 else (tp*tn-fp*fn)/d
def plugin_thr(p, ts):
    v=[mcc_c(p[p>t].sum(),(1-p[p<=t]).sum(),(1-p[p>t]).sum(),p[p<=t].sum()) for t in ts]
    return ts[int(np.argmax(v))]
ts=np.unique(np.round(np.linspace(0.02,0.95,200),4))
rng=np.random.default_rng(3)
for key in ["CYP3A4","CYP2D6"]:
    p=np.asarray(P[key]["p"],float); y=np.asarray(P[key]["y"],int); n=len(y)
    g=rng.integers(0,5,n)                       # внешние группы для честной калибровки
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
