"""Раскладывается ли совместное распределение на ТРИ множителя, как написано в спеке?
Проверяем два утверждения по отдельности:
  (A) множитель отбора p(R|y_scr) не зависит от theta и выпадает -- это про MAR;
  (B) p(y_drc, y_scr | theta) = p(y_drc|theta) * p(y_scr|theta) -- это про независимость каналов.
(B) неверно, когда у соединения есть собственный латентный остаток. Смотрим, что именно ломается."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np
rng=np.random.default_rng(7)
P=5; N=4000; REP=300
beta=np.array([1.0,-0.7,0.45,0.25,-0.35])
s_scr, s_drc = 0.40, 0.22

def one(seed, s_lat, thr_q=0.70):
    r=np.random.default_rng(seed)
    X=r.normal(0,1,(N,P)); Xs=np.hstack([np.ones((N,1)),X])
    a = X@beta + r.normal(0,s_lat,N)          # латентная активность соединения
    yscr = a + r.normal(0,s_scr,N)
    ydrc = a + r.normal(0,s_drc,N)
    R = yscr > np.quantile(yscr,thr_q)
    Xd=Xs[R]; yd=ydrc[R]
    # НАИВНАЯ (произведение двух маргиналов): веса = 1/маргинальную дисперсию канала
    wS=1/(s_lat**2+s_scr**2); wD=1/(s_lat**2+s_drc**2)
    A=wS*Xs.T@Xs + wD*Xd.T@Xd
    b=wS*Xs.T@yscr + wD*Xd.T@yd
    Ainv=np.linalg.inv(A); bn=Ainv@b
    se_model=np.sqrt(np.diag(Ainv))           # «модельная» ошибка, как если бы каналы были независимы
    # ЧЕСТНАЯ: у отобранных два коррелированных наблюдения одного латентного, обобщённый МНК 2x2
    Sig=np.array([[s_lat**2+s_scr**2, s_lat**2],[s_lat**2, s_lat**2+s_drc**2]])
    Si=np.linalg.inv(Sig); w1=Si.sum()        # вес для суммы, если оба канала есть
    # для отобранных: эффективное наблюдение = (Si@[yscr,ydrc]).sum()/w1 с весом w1
    yy=np.vstack([yscr[R],yd]); comb=(Si@yy).sum(0)/w1
    Aj=w1*Xd.T@Xd + wS*Xs[~R].T@Xs[~R]
    bj=w1*Xd.T@comb + wS*Xs[~R].T@yscr[~R]
    bj_=np.linalg.solve(Aj,bj)
    return bn[1:], bj_[1:], se_model[1:]

for s_lat in [0.0, 0.35, 0.70]:
    BN=[];BJ=[];SE=[]
    for s in range(REP):
        a,b,c=one(1000+s,s_lat); BN.append(a);BJ.append(b);SE.append(c)
    BN=np.array(BN);BJ=np.array(BJ);SE=np.array(SE)
    bias_n=np.abs(BN.mean(0)-beta).max(); bias_j=np.abs(BJ.mean(0)-beta).max()
    sd_true=BN.std(0,ddof=1).mean(); se_claim=SE.mean(0).mean()
    print(f"латентный остаток соединения sigma={s_lat:.2f}")
    print(f"   наивная (три множителя): макс. смещение {bias_n:.4f} | честная: {bias_j:.4f}")
    print(f"   истинное ско оценки {sd_true:.4f} | заявленная модельная ошибка {se_claim:.4f} "
          f"| занижение в {sd_true/se_claim:.2f} раза")
    print(f"   эффективность: ско честной / ско наивной = {BJ.std(0,ddof=1).mean()/sd_true:.3f}")
    print()
print("""Читается так. Смещения у наивной схемы нет ни при каком латентном остатке -- каждый
множитель по отдельности задан верно, а произведение верно заданных маргиналов
(составное правдоподобие) даёт состоятельную оценку. Ломается не точка, а её точность:
чем больше собственный латентный разброс соединения, тем сильнее наивная схема
считает два коррелированных наблюдения за два независимых и тем сильнее занижает
собственную погрешность.""")
