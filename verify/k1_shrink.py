"""Усадка к среднему: настоящий эффект или геометрия метрики. Четыре блока.

1. что усадка делает с ST-RAE, MAE, R2 и Спирменом по отдельности;
2. усадка к СМЕЩЁННОМУ центру — сколько выигрыша выживает и где оптимум;
3. усадка против честной перекалибровки вне выборки (a + b*yhat) и изотоники;
4. переносится ли всё это на подвыборку, обогащённую активными.

Итог: регрессионное затухание настоящее — независимый МНК вне фолда восстанавливает
a = (1-b)*mu с точностью до третьего знака, — но выигрыш почти целиком метрический:
MAE на 3A4 ухудшается (0.5348 -> 0.5409), R2 не двигается. На верхних 25% по активности
усадка портит все четыре фермента. Отсюда в submit.py она выключена по умолчанию.

Продолжение: k3_center.py разбирает, чем на самом деле является «сдвиг центра», а
k5/k6 меряют, насколько сильно тест вообще отличается от обучения."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd, json
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from sklearn.isotonic import IsotonicRegression
from scipy.stats import spearmanr
from cypsplit import butina_folds

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
rows = pd.read_csv(D + "rows.csv")
tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
        .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
oof = json.load(open(RES + "preds/oof.json"))
fold, _ = butina_folds(list(rows.SMILES))
LAM = np.linspace(0.30, 1.30, 101)


def pack(c):
    col = f"{c}_pIC50_direct_inhibition"; m = tr[col].notna().to_numpy()
    return (tr.loc[m, col].to_numpy(), tr.loc[m, col + "_conf_low"].to_numpy(),
            tr.loc[m, col + "_conf_high"].to_numpy(),
            np.asarray(oof[f"FP+DESC+MECH|{c}"]), fold[m])


def S(y, p, lo, hi):
    return strae(y, p, y_true_upper=hi, y_true_lower=lo)


print("=== 1. Что усадка делает с разными метриками ===")
print(f"{'фермент':8s} {'l*':>5s} {'ST-RAE':>16s} {'MAE':>16s} {'R2':>16s} {'rho':>8s}")
SH = {}
for c in CYPS:
    y, lo, hi, p, f = pack(c)
    ps = np.zeros_like(p); lams = []
    for k in range(5):
        a, b = f != k, f == k
        if b.sum() == 0:
            continue
        mu = y[a].mean()
        best = LAM[np.argmin([S(y[a], mu + L * (p[a] - mu), lo[a], hi[a]) for L in LAM])]
        lams.append(best); ps[b] = mu + best * (p[b] - mu)
    s0, s1 = S(y, p, lo, hi), S(y, ps, lo, hi)
    m0, m1 = np.abs(y - p).mean(), np.abs(y - ps).mean()
    r0 = 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    r1 = 1 - ((y - ps) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    SH[c] = (y, lo, hi, p, ps, f)
    print(f"{c:8s} {np.mean(lams):5.2f} {s0:7.4f}->{s1:7.4f} {m0:7.4f}->{m1:7.4f} "
          f"{r0:7.4f}->{r1:7.4f} {spearmanr(p,y).statistic:8.3f}")

print("\n=== 2. Усадка к СМЕЩЁННОМУ центру: сколько выигрыша выживает ===")
print(f"{'сдвиг центра':14s} " + " ".join(f"{c[3:]:>8s}" for c in CYPS) + f" {'макро':>8s}")
for off in [0.0, -0.25, +0.25, -0.5, +0.5, -1.0, +1.0]:
    v = []
    for c in CYPS:
        y, lo, hi, p, f = pack(c); ps = np.zeros_like(p)
        for k in range(5):
            a, b = f != k, f == k
            if b.sum() == 0:
                continue
            mu = y[a].mean() + off
            best = LAM[np.argmin([S(y[a], mu + L * (p[a] - mu), lo[a], hi[a]) for L in LAM])]
            ps[b] = mu + best * (p[b] - mu)
        v.append(S(y, ps, lo, hi))
    print(f"{off:+14.2f} " + " ".join(f"{x:8.4f}" for x in v) + f" {np.mean(v):8.4f}")

print("\n=== 3. Усадка против честной перекалибровки (a + b*yhat) и изотоники, вне выборки ===")
print(f"{'фермент':8s} {'сырое':>8s} {'усадка':>8s} {'a+b*yh':>8s} {'изотоника':>10s} "
      f"{'b':>6s} {'a':>7s} {'(1-b)mu':>8s}")
for c in CYPS:
    y, lo, hi, p, f = pack(c)
    pl = np.zeros_like(p); pi_ = np.zeros_like(p); bs = []; as_ = []; mus = []
    for k in range(5):
        a, b = f != k, f == k
        if b.sum() == 0:
            continue
        B, A = np.polyfit(p[a], y[a], 1); bs.append(B); as_.append(A); mus.append(y[a].mean())
        pl[b] = A + B * p[b]
        pi_[b] = IsotonicRegression(out_of_bounds="clip").fit(p[a], y[a]).predict(p[b])
    print(f"{c:8s} {S(y,p,lo,hi):8.4f} {S(y,SH[c][4],lo,hi):8.4f} {S(y,pl,lo,hi):8.4f} "
          f"{S(y,pi_,lo,hi):10.4f} {np.mean(bs):6.3f} {np.mean(as_):7.3f} "
          f"{(1-np.mean(bs))*np.mean(mus):8.3f}")

print("\n=== 4. Переносится ли усадка на подвыборку, обогащённую активными ===")
print("   верхний хвост по ИСТИННОЙ активности — грубое приближение геометрии теста")
print(f"{'подвыборка':22s} " + " ".join(f"{c[3:]:>16s}" for c in CYPS))
for q, lab in [(0, "вся выборка"), (50, "верхние 50% по y"),
               (75, "верхние 25% по y"), (90, "верхние 10% по y")]:
    out = []
    for c in CYPS:
        y, lo, hi, p, ps, f = SH[c]
        m = y >= np.percentile(y, q) if q else np.ones(len(y), bool)
        out.append(f"{S(y[m],p[m],lo[m],hi[m]):7.4f}->{S(y[m],ps[m],lo[m],hi[m]):7.4f}")
    print(f"{lab:22s} " + " ".join(f"{x:>16s}" for x in out))
