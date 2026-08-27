"""Насколько неверно правило sigma=(граница-mu)/1.96 для расщеплённой нормали,
и какая замена правильная. Считаем на настоящих conf_low/conf_high."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import numpy as np, pandas as pd
from scipy.stats import norm
from scipy.optimize import brentq

def quantiles(mu, sm, sp, q):
    """Точные квантили расщеплённой нормали с модой mu и полусигмами sm, sp."""
    w = sm/(sm+sp)                       # масса слева от моды
    out = np.empty_like(np.asarray(q,float))
    q = np.atleast_1d(np.asarray(q,float))
    left = q < w
    out = np.where(left, mu + sm*norm.ppf(np.clip(q/(2*w),1e-15,1-1e-15)),
                         mu + sp*norm.ppf(np.clip(0.5 + (q-w)/(2*(1-w)),1e-15,1-1e-15)))
    return out

def fit_exact(mu, lo, hi, alpha=0.025):
    """Подбираем sm, sp так, чтобы квантили alpha и 1-alpha попали ровно в lo и hi."""
    a, b = mu-lo, hi-mu
    if a <= 0 or b <= 0: return np.nan, np.nan
    def resid(logr):
        r = np.exp(logr); w = 1/(1+r)
        zl = -norm.ppf(np.clip(alpha/(2*w),1e-12,1-1e-12))     # >0
        zu =  norm.ppf(np.clip(1-alpha/(2*(1-w)),1e-12,1-1e-12))
        if not np.isfinite(zl) or not np.isfinite(zu) or zl<=0 or zu<=0: return 1e6
        return np.log((b/zu)/(a/zl)) - logr
    try: logr = brentq(resid, -6, 6, xtol=1e-10)
    except Exception: return np.nan, np.nan
    r = np.exp(logr); w = 1/(1+r)
    zl = -norm.ppf(alpha/(2*w)); zu = norm.ppf(1-alpha/(2*(1-w)))
    return a/zl, b/zu

print("Учебный пример из документа: mu=4.6, полоса [4.10, 5.30]")
mu, lo, hi = 4.6, 4.10, 5.30
sm_naive, sp_naive = (mu-lo)/1.96, (hi-mu)/1.96
q = quantiles(mu, sm_naive, sp_naive, [0.025,0.975])
print(f"  правило /1.96 :  sm={sm_naive:.4f} sp={sp_naive:.4f} -> фактическая полоса [{q[0]:.4f}, {q[1]:.4f}]")
sm_ex, sp_ex = fit_exact(mu, lo, hi)
q2 = quantiles(mu, sm_ex, sp_ex, [0.025,0.975])
print(f"  точная подгонка: sm={sm_ex:.4f} sp={sp_ex:.4f} -> фактическая полоса [{q2[0]:.4f}, {q2[1]:.4f}]")
print(f"  поправочные множители: снизу x{sm_ex/sm_naive:.4f}, сверху x{sp_ex/sp_naive:.4f}")

CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
print("\nНа настоящих метках (медианы по ферменту):")
rowsout=[]
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; d=tr[tr[col].notna()]
    y=d[col].to_numpy(); l=d[col+"_conf_low"].to_numpy(); h=d[col+"_conf_high"].to_numpy()
    m=(y-l>1e-9)&(h-y>1e-9)
    y,l,h=y[m],l[m],h[m]
    sm=(y-l)/1.96; sp=(h-y)/1.96
    qq=np.array([quantiles(y[i],sm[i],sp[i],[0.025,0.975]) for i in range(len(y))])
    cov_naive=np.mean((qq[:,0]<=l+1e-9)&(qq[:,1]>=h-1e-9))
    ex=np.array([fit_exact(y[i],l[i],h[i]) for i in range(len(y))])
    ok=np.isfinite(ex).all(1)
    rowsout.append(dict(cyp=c, n=int(m.sum()),
        err_low_med=round(float(np.median(qq[:,0]-l)),4),
        err_high_med=round(float(np.median(qq[:,1]-h)),4),
        width_ratio_med=round(float(np.median((qq[:,1]-qq[:,0])/(h-l))),4),
        corr_lo=round(float(np.median(ex[ok,0]/sm[ok])),4),
        corr_hi=round(float(np.median(ex[ok,1]/sp[ok])),4)))
print(pd.DataFrame(rowsout).to_string(index=False))
print("""
err_low_med  -- насколько нижняя граница модели смещена относительно conf_low (>0 = полоса модели уже)
width_ratio  -- ширина 95%-интервала модели / ширина заявленной полосы
corr_lo/hi   -- множители, на которые надо домножить (граница-mu)/1.96, чтобы стало точно""")
