"""Перевзвешивание обучения под тестовое распределение — по одной оси, и что оно даёт.

Веса в k5 вырождались, потому что оценивались в 2295 измерениях. Здесь отношение
плотностей оценивается вдоль ОДНОЙ оси — самой оси предсказаний, — и именно вдоль неё
работает аффинная перекалибровка, так что для вопроса «переносится ли (c, l)» этого
достаточно. Эффективный размер выборки выходит 800-1450 вместо 116.

Что получается: оптимальные (off, l) под тестоподобным распределением почти не двигаются,
порядок стратегий не двигается вообще. Сдвиг, измеренный в k5, мягкий, а мягкий сдвиг не
сдвигает оптимум. Это отменяет мою же более раннюю формулировку «усадка не переносится»:
стресс-тест k1 (верхние 25%) двигал ybar на +1.1...+1.3, то есть был втрое сильнее того
сдвига, который на тесте действительно есть.

Оговорка, без которой вывод неверен: взвешивание правит только маргинал предсказаний и
предполагает, что p(y | yhat) на тесте та же. Это ровно то допущение, которое сама
перекалибровка и должна проверять, так что бесплатным обедом здесь не пахнет.

Во второй половине — совпадение двух независимых порядков по ферментам."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd, json
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from scipy.stats import spearmanr

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
rows = pd.read_csv(D + "rows.csv")
tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
        .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
oof = json.load(open(RES + "preds/oof.json"))
pred_path = _pl.Path(D) / "test_pred.npz"
if not pred_path.exists():
    raise SystemExit("нет data/test_pred.npz — сначала verify/k5_shift.py")
PT = np.load(pred_path)
OFF = np.round(np.arange(-1.0, 2.01, 0.05), 2)
LAM = np.round(np.arange(0.30, 1.31, 0.02), 2)


def Sw(y, p, lo, hi, w):
    n = (w * (np.maximum(p - hi, 0) + np.maximum(lo - p, 0))).sum()
    yb = (w * y).sum() / w.sum()
    return n / (w * (np.maximum(yb - hi, 0) + np.maximum(lo - yb, 0))).sum()


print("=" * 88)
print("Веса = p_тест(yhat) / p_обуч(yhat), 20 квантильных корзин, обрезка на 10")
print("=" * 88)
print(f"{'фермент':8s} {'эфф.n':>7s} {'ybar_w':>7s} {'ybar':>6s} | {'сырое':>7s} {'к mu':>7s} "
      f"{'к mu+off':>8s} {'карта_w':>8s} | {'off_w':>6s} {'l_w':>5s}  {'off_0':>6s} {'l_0':>5s}")
T = {k: 0.0 for k in "rABW"}
for c in CYPS:
    col = f"{c}_pIC50_direct_inhibition"; m = tr[col].notna().to_numpy()
    y = tr.loc[m, col].to_numpy()
    lo = tr.loc[m, col + "_conf_low"].to_numpy(); hi = tr.loc[m, col + "_conf_high"].to_numpy()
    p = np.asarray(oof[f"FP+DESC+MECH|{c}"]); pt = PT[c]; mu = y.mean()
    e = np.quantile(p, np.linspace(0, 1, 21)); e[0] -= 1e-6; e[-1] += 1e-6
    ct, _ = np.histogram(pt, bins=e); co, _ = np.histogram(p, bins=e)
    r = (ct / max(ct.sum(), 1)) / np.maximum(co / co.sum(), 1e-9)
    w = np.clip(r[np.clip(np.digitize(p, e) - 1, 0, 19)], 0, 10); w = w / w.mean()
    lA = LAM[np.argmin([strae(y, mu + L * (p - mu), y_true_upper=hi, y_true_lower=lo)
                        for L in LAM])]
    G0 = np.array([[strae(y, mu + o + L * (p - mu - o), y_true_upper=hi, y_true_lower=lo)
                    for L in LAM] for o in OFF])
    i0, j0 = np.unravel_index(G0.argmin(), G0.shape)
    GW = np.array([[Sw(y, mu + o + L * (p - mu - o), lo, hi, w) for L in LAM] for o in OFF])
    iw, jw = np.unravel_index(GW.argmin(), GW.shape)
    v = [Sw(y, p, lo, hi, w), Sw(y, mu + lA * (p - mu), lo, hi, w),
         Sw(y, mu + OFF[i0] + LAM[j0] * (p - mu - OFF[i0]), lo, hi, w), GW[iw, jw]]
    for k, x in zip("rABW", v):
        T[k] += x / 4
    print(f"{c:8s} {w.sum()**2/(w**2).sum():7.0f} {(w*y).sum()/w.sum():7.3f} {mu:6.3f} | "
          f"{v[0]:7.4f} {v[1]:7.4f} {v[2]:8.4f} {v[3]:8.4f} | {OFF[iw]:+6.2f} {LAM[jw]:5.2f} "
          f" {OFF[i0]:+6.2f} {LAM[j0]:5.2f}")
print(f"{'МАКРО':8s} {'':7s} {'':7s} {'':6s} | {T['r']:7.4f} {T['A']:7.4f} {T['B']:8.4f} "
      f"{T['W']:8.4f}")

print()
print("=" * 88)
print("Согласие двух независимых порядков по ферментам")
print("=" * 88)
# оба числа измерены отдельно: первое в src/trunkscore.py, второе там же по прогонам
RHO_SCR = {"CYP3A4": 0.936, "CYP2C9": 0.896, "CYP1A2": 0.862, "CYP2D6": 0.828}
DSTRAE = {"CYP3A4": -0.027, "CYP2C9": -0.025, "CYP1A2": +0.047, "CYP2D6": +0.120}
shift = {c: float(PT[c].mean() - np.asarray(oof[f"FP+DESC+MECH|{c}"]).mean()) for c in CYPS}
print(f"{'фермент':8s} {'|rho| скрин~pIC50':>18s} {'сдвиг предск. тест':>19s} "
      f"{'d ST-RAE при l=3':>18s}")
for c in sorted(CYPS, key=lambda x: -RHO_SCR[x]):
    print(f"{c:8s} {RHO_SCR[c]:18.3f} {shift[c]:+19.3f} {DSTRAE[c]:+18.3f}")
a = [RHO_SCR[c] for c in CYPS]; b = [shift[c] for c in CYPS]; d = [DSTRAE[c] for c in CYPS]
print(f"\nСпирмен |rho| ~ сдвиг     = {spearmanr(a,b).statistic:+.3f}")
print(f"Спирмен сдвиг ~ d ST-RAE  = {spearmanr(b,d).statistic:+.3f}")
print(f"Спирмен |rho| ~ d ST-RAE  = {spearmanr(a,d).statistic:+.3f}")
print("\nn = 4, перестановочная p для |rho| = 1 равна 1/24 = 0.042. Это один бит, не закон:")
print("четыре точки, и |rho| скрининга коррелирует со всем, что меняется между ферментами.")
print("Проверка внутри фермента — подмешивать шум в скрининговый канал ступенями и смотреть,")
print("монотонно ли ползёт d ST-RAE, — превращает n=4 в дозовую кривую на каждом отдельно.")
