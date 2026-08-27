"""Сколько на самом деле стоит выбор порога plug-in методом вместо настоящего оптимума.
Считаем на СОХРАНЁННЫХ вероятностях TDI-классификатора, не на симуляции."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, json
from sklearn.metrics import matthews_corrcoef

P=json.load(open(RES+"preds/tdi_probs.json"))
print("ключи:", list(P.keys())[:8])
def mcc_from_counts(tp,tn,fp,fn):
    d=np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return 0.0 if d<=0 else (tp*tn-fp*fn)/d
for key in P:
    d=P[key]
    if not isinstance(d,dict) or "p" not in d: print("  пропуск",key,type(d)); continue
    p=np.asarray(d["p"],float); y=np.asarray(d["y"],int)
    ts=np.unique(np.round(np.linspace(0.02,0.95,200),4))
    emp=np.array([matthews_corrcoef(y,(p>t).astype(int)) for t in ts])
    plug=np.array([mcc_from_counts(p[p>t].sum(),(1-p[p<=t]).sum(),(1-p[p>t]).sum(),p[p<=t].sum()) for t in ts])
    t_emp=ts[emp.argmax()]; t_plug=ts[plug.argmax()]
    print(f"\n{key}: n={len(y)}, доля положительных {y.mean():.3f}, калибровка E[p]={p.mean():.3f}")
    print(f"  порог по plug-in E[MCC]: {t_plug:.3f} -> реальный MCC {emp[np.argmin(abs(ts-t_plug))]:.4f}")
    print(f"  оптимальный порог:       {t_emp:.3f} -> реальный MCC {emp.max():.4f}")
    print(f"  порог 0.5:                     -> реальный MCC {emp[np.argmin(abs(ts-0.5))]:.4f}")
    print(f"  ЦЕНА plug-in вместо оптимума: {emp.max()-emp[np.argmin(abs(ts-t_plug))]:+.4f}")
    print(f"  прирост plug-in против 0.5:   {emp[np.argmin(abs(ts-t_plug))]-emp[np.argmin(abs(ts-0.5))]:+.4f}")
