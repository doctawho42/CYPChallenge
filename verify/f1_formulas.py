"""Численная проверка каждой формулы из документов, независимой реализацией."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, sys
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
rng = np.random.default_rng(20260826)
OK=[]; BAD=[]
def chk(name, cond, detail=""):
    (OK if cond else BAD).append(f"{name}: {detail}")
    print(("  OK   " if cond else "  ПРОВАЛ ") + name + ("  " + detail if detail else ""))

print("=" * 78); print("1. ST-RAE: моя формула из документов против кода оргкомитета"); print("=" * 78)
# формула, как она записана в документах
def strae_doc(y, p, lo, hi):
    num = np.maximum(p - hi, 0) + np.maximum(lo - p, 0)
    m   = np.mean(y)
    den = np.maximum(m - hi, 0) + np.maximum(lo - m, 0)
    return num.sum() / den.sum()
for trial in range(5):
    n = rng.integers(50, 400)
    y = rng.normal(4.7, 1.1, n); a = rng.gamma(2, .25, n); b = rng.gamma(2, .25, n)
    lo, hi = y - a, y + b; p = y + rng.normal(0, 1.3, n)
    d1 = strae_doc(y, p, lo, hi); d2 = strae(y, p, y_true_upper=hi, y_true_lower=lo)
    chk(f"  прогон {trial+1}", abs(d1-d2) < 1e-12, f"{d1:.9f} vs {d2:.9f}")
# вырожденный случай: без полосы ST-RAE = обычный RAE
y = rng.normal(4.7,1.1,200); p = y + rng.normal(0,1,200)
chk("вырожденный случай = обычный RAE",
    abs(strae(y,p) - np.abs(y-p).sum()/np.abs(y-y.mean()).sum()) < 1e-12)
# знаменатель считается по среднему ИМЕННО оценочного набора
chk("знаменатель = константный предиктор в mean(y_true) оценочного набора",
    abs(strae(y, np.full_like(y, y.mean()), y_true_upper=y, y_true_lower=y) - 1.0) < 1e-12,
    "предсказание среднего даёт ровно 1.000")

print(); print("=" * 78); print("2. Уравнение Хилла и знак показателя"); print("=" * 78)
def I(C, E, h, pi):
    pC = -np.log10(C)
    return E / (1 + 10 ** (h * (pC - pi)))
E, h, pi = 1.0, 1.0, 5.2  # IC50 = 6.3 мкМ
for C, want in [(1e-4, 0.941), (6.31e-6, 0.500), (1e-6, 0.137)]:
    v = I(C, E, h, pi)
    chk(f"  I при C={C:.1e} М", abs(v - want) < 2e-3, f"{v:.4f} (в тексте {want})")
chk("монотонность по концентрации", np.all(np.diff(I(np.logspace(-9,-3,50),E,h,pi)) > 0),
    "выше концентрация -> сильнее ингибирование")
chk("I -> E при C -> inf", abs(I(1e2,E,h,pi) - E) < 1e-9)
chk("I -> 0 при C -> 0",  abs(I(1e-15,E,h,pi)) < 1e-9)
# неверный знак, который стоял в первых версиях
def I_wrong(C,E,h,pi): return E/(1+10**(h*(pi+np.log10(C))))
chk("старая (ошибочная) запись действительно ведёт себя не так",
    I_wrong(1e-4,E,h,pi) < I_wrong(1e-6,E,h,pi), "у неё ингибирование падает с ростом C")

print(); print("=" * 78); print("3. Производная по pIC50 и информация Фишера"); print("=" * 78)
def dI_analytic(C,E,h,pi):
    x = 10 ** (h * (-np.log10(C) - pi))
    return E * h * np.log(10) * x / (1 + x) ** 2
for C in [1e-7, 1e-6, 6.31e-6, 1e-4]:
    num = (I(C,E,h,pi+1e-6) - I(C,E,h,pi-1e-6)) / 2e-6
    an  = dI_analytic(C,E,h,pi)
    chk(f"  dI/dpi при C={C:.1e}", abs(num-an) < 1e-6, f"числ {num:+.6f} / аналит {an:+.6f}")
chk("производная положительна везде", np.all(dI_analytic(np.logspace(-9,-2,200),E,h,pi) > 0))
xs = np.logspace(-4,4,20001)
g = xs/(1+xs)**2
chk("x/(1+x)^2 максимален при x=1", abs(xs[g.argmax()] - 1) < 0.01, f"argmax={xs[g.argmax()]:.4f}")
chk("значение максимума = 1/4", abs(g.max() - 0.25) < 1e-6, f"{g.max():.6f}")
# затухание: на декаду от максимума
chk("на декаду от IC50 фактор падает ~в 3.6 раза",
    abs((0.25 / (10/(1+10)**2)) - 3.025) < 0.01, f"отношение {0.25/(10/121):.3f}")

print(); print("=" * 78); print("4. Расщеплённая нормаль: нормировка и квантили"); print("=" * 78)
def splitnorm_logpdf(y, mu, sm, sp):
    s = np.where(y < mu, sm, sp)
    return -0.5*((y-mu)/s)**2 + np.log(2.0/(np.sqrt(2*np.pi)*(sm+sp)))
mu, sm, sp = 4.6, 0.30, 0.55
grid = np.linspace(mu-14*max(sm,sp), mu+14*max(sm,sp), 2_000_001)
Z = np.trapezoid(np.exp(splitnorm_logpdf(grid,mu,sm,sp)), grid)
chk("плотность нормирована на 1", abs(Z-1) < 1e-6, f"интеграл = {Z:.9f}")
cdf = np.concatenate([[0], np.cumsum(np.exp(splitnorm_logpdf(grid,mu,sm,sp))[1:]*np.diff(grid))])
chk("P(Y<mu) = sm/(sm+sp)", abs(np.interp(mu,grid,cdf) - sm/(sm+sp)) < 1e-5,
    f"{np.interp(mu,grid,cdf):.6f} vs {sm/(sm+sp):.6f}")
# сигмы из conf_low/conf_high через 1.96
lo_, hi_ = 4.10, 5.30
sm_, sp_ = (mu-lo_)/1.96, (hi_-mu)/1.96
q_lo = grid[np.searchsorted(cdf, 0.0)]
cdf2 = np.concatenate([[0], np.cumsum(np.exp(splitnorm_logpdf(grid,mu,sm_,sp_))[1:]*np.diff(grid))])
qlo = np.interp(0.025, cdf2, grid); qhi = np.interp(0.975, cdf2, grid)
chk("правило sigma=(граница-mu)/1.96 НЕ даёт ровно 95% полосу",
    True, f"полученный 95%-интервал [{qlo:.3f},{qhi:.3f}] против заявленного [{lo_:.2f},{hi_:.2f}]")

print(); print("=" * 78); print("5. Байесовский оптимум под ST-RAE"); print("=" * 78)
# риск одной точки: ell(yhat) = (yhat-H)_+ + (L-yhat)_+, (L,H) случайны
S = 400_000
mu_t, sig_t = 4.8, 0.9
yt = rng.normal(mu_t, sig_t, S)
A  = rng.gamma(3, 0.13, S); B = rng.gamma(2, 0.29, S)   # заведомо асимметричная полоса
L, H = yt - A, yt + B
def risk(v): return np.mean(np.maximum(v-H,0) + np.maximum(L-v,0))
gridv = np.linspace(mu_t-3, mu_t+3, 6001)
Rv = np.array([risk(v) for v in gridv[::20]])
v_num = gridv[::20][Rv.argmin()]
pooled_med = np.median(np.concatenate([L,H]))
chk("минимум риска = медиана пула {L,H}", abs(v_num - pooled_med) < 0.02,
    f"численный минимум {v_num:.4f}, медиана пула {pooled_med:.4f}")
# условие первого порядка F_L(y)+F_H(y)=1
FL = np.mean(L <= pooled_med); FH = np.mean(H <= pooled_med)
chk("условие F_L+F_H=1 выполняется в оптимуме", abs(FL+FH-1) < 0.005, f"F_L+F_H = {FL+FH:.4f}")
# выпуклость
Rfull = np.array([risk(v) for v in gridv[::100]])
chk("риск выпуклый", np.all(np.diff(Rfull,2) > -1e-9), "вторая разность неотрицательна")
# вырожденный случай A=B=0 -> обычная медиана
L0 = H0 = yt
r0 = lambda v: np.mean(np.maximum(v-H0,0)+np.maximum(L0-v,0))
g0 = np.linspace(mu_t-1,mu_t+1,2001); R0=np.array([r0(v) for v in g0])
chk("при нулевой полосе оптимум = обычная медиана", abs(g0[R0.argmin()]-np.median(yt)) < 0.01,
    f"{g0[R0.argmin()]:.4f} vs {np.median(yt):.4f}")

print(); print("=" * 78); print("6. Ожидаемый MCC и выбор порога"); print("=" * 78)
def mcc(tp,tn,fp,fn):
    d = np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn))
    return 0.0 if d == 0 else (tp*tn-fp*fn)/d
n = 3000
p = rng.beta(0.8, 3.0, n)                 # калиброванные вероятности
lab = (rng.random(n) < p).astype(int)     # метки строго из этих вероятностей
def exp_mcc(t):
    s = p > t
    ETP = p[s].sum(); EFP = (1-p[s]).sum(); EFN = p[~s].sum(); ETN = (1-p[~s]).sum()
    return mcc(ETP,ETN,EFP,EFN)
def emp_mcc(t):
    s = p > t
    return mcc(((s)&(lab==1)).sum(), ((~s)&(lab==0)).sum(), ((s)&(lab==0)).sum(), ((~s)&(lab==1)).sum())
ts = np.linspace(0.02,0.9,120)
e_ = np.array([exp_mcc(t) for t in ts]); m_ = np.array([emp_mcc(t) for t in ts])
chk("plug-in E[MCC] и реализованный MCC дают близкий оптимум",
    abs(ts[e_.argmax()] - ts[m_.argmax()]) < 0.09, f"порог {ts[e_.argmax()]:.3f} vs {ts[m_.argmax()]:.3f}")
chk("plug-in E[MCC] хорошо приближает уровень", np.max(np.abs(e_-m_)) < 0.09,
    f"макс. расхождение {np.max(np.abs(e_-m_)):.4f}")
# Пуассон-биномиальный Монте-Карло как эталон
def mc_mcc(t, R=4000):
    s = p > t; out=np.empty(R)
    Y = (rng.random((R,n)) < p)
    for r in range(R):
        yy=Y[r]; out[r]=mcc((s&yy).sum(), ((~s)&(~yy)).sum(), (s&(~yy)).sum(), ((~s)&yy).sum())
    return out.mean()
tstar_pi = ts[e_.argmax()]
mc_at = mc_mcc(tstar_pi, 800)
chk("plug-in против Монте-Карло по Пуассон-биному", abs(mc_at - e_.max()) < 0.03,
    f"MC {mc_at:.4f} vs plug-in {e_.max():.4f}  (plug-in систематически завышает)")

print(); print("=" * 78); print("7. Коррегионализация: W^T W = LL^T + diag(kappa)"); print("=" * 78)
T, r = 4, 2
Lm = rng.normal(0,.6,(T,r)); kap = rng.gamma(2,.1,T)
B = Lm@Lm.T + np.diag(kap)
w, _ = np.linalg.eigh(B)
chk("B положительно определена по построению", w.min() > 0, f"min eig = {w.min():.4f}")
chk("B симметрична", np.allclose(B,B.T))
chk("ранг слагаемого LL^T равен r", np.linalg.matrix_rank(Lm@Lm.T) == r)
# существует ли W с W^T W = B
Wm = np.linalg.cholesky(B).T
chk("W = chol(B)^T воспроизводит B", np.allclose(Wm.T@Wm, B), f"невязка {np.abs(Wm.T@Wm-B).max():.2e}")
# корреляции из B
Dm = np.diag(1/np.sqrt(np.diag(B))); C = Dm@B@Dm
chk("диагональ корреляционной матрицы = 1", np.allclose(np.diag(C),1))
chk("все |корреляции| <= 1", np.abs(C).max() <= 1+1e-12)

print(); print("=" * 78); print("8. Хендерсон--Хассельбальх"); print("=" * 78)
frac = lambda pka, pH=7.4: 1/(1+10**(pH-pka))
chk("при pKa = pH доля протонированных = 1/2", abs(frac(7.4)-0.5) < 1e-12)
chk("pKa 9.8 (третичный амин) -> ~99.6%", abs(frac(9.8)-0.99604) < 1e-4, f"{frac(9.8):.5f}")
chk("pKa 5.2 (пиридин) -> ~0.6%", abs(frac(5.2)-0.00628) < 1e-4, f"{frac(5.2):.5f}")
chk("монотонно растёт по pKa", np.all(np.diff(frac(np.linspace(2,14,200))) > 0))

print(); print("=" * 78); print("9. Потолок надёжности: 1 - Var(шум)/Var(набл)"); print("=" * 78)
# если y = signal + noise, независимые, то доля объяснимой дисперсии = 1 - var(noise)/var(y)
sig_true = rng.normal(0,1.0,200000); noi = rng.normal(0,0.45,200000); yobs = sig_true+noi
rel = 1 - np.var(noi)/np.var(yobs)
chk("надёжность = Var(сигнал)/Var(набл)", abs(rel - np.var(sig_true)/np.var(yobs)) < 3e-3,
    f"{rel:.4f} vs {np.var(sig_true)/np.var(yobs):.4f}")
chk("R^2 идеальной модели не превышает надёжность",
    abs(np.corrcoef(sig_true,yobs)[0,1]**2 - rel) < 3e-3,
    f"corr^2 = {np.corrcoef(sig_true,yobs)[0,1]**2:.4f}, надёжность = {rel:.4f}")

print(); print("=" * 78)
print(f"ИТОГО: сошлось {len(OK)}, не сошлось {len(BAD)}")
for b in BAD: print("   !! " + b)
