"""Центральное теоретическое утверждение документов, до сих пор не проверявшееся:
отбор в подробное измерение зависит от НАБЛЮДАЕМОГО скрининга, поэтому при совместном
правдоподобии множитель отбора не зависит от theta и выпадает из градиента.
Следствие: подгонка только по отобранным метках смещена, совместная -- нет.

Проверяем симуляцией, где истина известна."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np
rng=np.random.default_rng(2026)
P=6; N=6000; REP=40
beta=np.array([1.2,-0.8,0.5,0.0,0.3,-0.4])
s_lat, s_scr, s_drc = 0.7, 0.45, 0.25

def run(seed, thr_q=0.72):
    r=np.random.default_rng(seed)
    Xd=r.normal(0,1,(N,P))
    a = Xd@beta + r.normal(0,s_lat,N)             # латентная активность
    yscr = a + r.normal(0,s_scr,N)                # скрининг, есть у всех
    thr = np.quantile(yscr, thr_q)
    R = yscr > thr                                # отбор ТОЛЬКО по наблюдаемому скринингу
    ydrc = a + r.normal(0,s_drc,N)
    # A: обычная регрессия по отобранным
    XA=np.hstack([np.ones((R.sum(),1)),Xd[R]])
    bA=np.linalg.lstsq(XA,ydrc[R],rcond=None)[0][1:]
    # B: совместное правдоподобие -- скрининг по всем + DRC где есть, общий линейный предиктор.
    #    Взвешенный МНК: вес = 1/дисперсия соответствующего канала.
    wS=1/(s_lat**2+s_scr**2); wD=1/(s_lat**2+s_drc**2)
    Xs=np.hstack([np.ones((N,1)),Xd]); Xdd=np.hstack([np.ones((R.sum(),1)),Xd[R]])
    A_=wS*Xs.T@Xs + wD*Xdd.T@Xdd
    b_=wS*Xs.T@yscr + wD*Xdd.T@ydrc[R]
    bB=np.linalg.solve(A_,b_)[1:]
    # C: контроль -- регрессия только по скринингу, все строки
    bC=np.linalg.lstsq(Xs,yscr,rcond=None)[0][1:]
    # D: отбор, зависящий от ЛАТЕНТНОЙ величины напрямую (настоящий MNAR) -- должно ломаться и у B
    R2 = a > np.quantile(a,thr_q)
    Xd2=np.hstack([np.ones((R2.sum(),1)),Xd[R2]])
    A2=wS*Xs.T@Xs + wD*Xd2.T@Xd2
    b2=wS*Xs.T@np.where(R2,yscr,yscr) + wD*Xd2.T@ydrc[R2]
    # при MNAR скрининг тоже надо наблюдать у всех; ломаем иначе: скрининг виден только у отобранных
    Xs2=np.hstack([np.ones((R2.sum(),1)),Xd[R2]])
    A3=wS*Xs2.T@Xs2 + wD*Xd2.T@Xd2
    b3=wS*Xs2.T@yscr[R2] + wD*Xd2.T@ydrc[R2]
    bD=np.linalg.solve(A3,b3)[1:]
    return bA,bB,bC,bD,R.mean()

res={k:[] for k in "ABCD"}; frac=[]
for s in range(REP):
    bA,bB,bC,bD,f=run(1000+s); frac.append(f)
    res["A"].append(bA); res["B"].append(bB); res["C"].append(bC); res["D"].append(bD)
print(f"доля соединений, попавших в подробное измерение: {np.mean(frac):.3f}")
print(f"истинные коэффициенты: {beta}")
print()
print(f"{'метод':52s} {'смещение по норме':>18s} {'наклон калибровки':>18s}")
names={"A":"A. только отобранные DRC-метки (как сейчас)",
       "B":"B. совместное правдоподобие, отбор по скринингу",
       "C":"C. только скрининг, все строки (контроль сверху)",
       "D":"D. отбор по ЛАТЕНТНОЙ величине (настоящий MNAR)"}
for k in "ABCD":
    Bm=np.array(res[k]); bias=Bm.mean(0)-beta
    slope=np.sum(Bm.mean(0)*beta)/np.sum(beta*beta)
    print(f"{names[k]:52s} {np.linalg.norm(bias):18.4f} {slope:18.4f}")
print()
Bm=np.array(res["A"]); print("A покоэффициентно, смещение:", np.round(Bm.mean(0)-beta,3))
Bm=np.array(res["B"]); print("B покоэффициентно, смещение:", np.round(Bm.mean(0)-beta,3))
Bm=np.array(res["D"]); print("D покоэффициентно, смещение:", np.round(Bm.mean(0)-beta,3))
print("""
Читается так: если отбор смотрит только на скрининг, который мы видим у всех,
совместное правдоподобие возвращает истинные коэффициенты (наклон 1.00), а подгонка
по одним отобранным меткам сжимает их. Если отбор смотрит на скрытую истину, а
скрининг у неотобранных не наблюдаем, ломается и совместная схема -- то есть
аргумент держится ровно на том, что первая ступень отбора наблюдаема.""")
