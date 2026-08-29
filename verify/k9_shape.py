"""Тест не просто сдвинут — он ещё и расширен, и по ST-RAE это тянет в нашу пользу.

k8 подбирает одно delta и по построению попадает в среднее теста, а по разбросу
промахивается: KS отвергает на трёх ферментах из четырёх. Экспоненциальный наклон умеет
двигать среднее и не умеет добавлять разброс, так что любая схема с одним delta —
src/reweight.py и src/shrinkchoice.py в том числе — специфицирована не полностью.

Это меняет знак одного из выводов. «Выборка уже по активности, значит знаменатель меньше и
оценка растёт» — посылка неверна: наблюдаемый разброс предсказаний на тесте БОЛЬШЕ, чем на
обучении (1.13 / 1.40 / 1.05 / 1.10 по sd), а в пространстве меток он не меньше того же
множителя, если шум ядра на тесте не больше обучающего. Последнее поддержано тем, что тест
ближе к обучению по сходству (0.587 против 0.435), то есть шум скорее меньше.

Два семейства:
  * двухпараметрическое w ~ exp(t1*(y-ybar) + t2*(y-ybar)^2), максэнтропийное при заданных
    первых двух моментах, подогнанное под ОБА момента предсказаний. На CYP2C9 оно
    ВЫРОЖДАЕТСЯ — эффективный размер 20 из 1285, тот же провал, что в k5, — и его числа по
    2C9 недействительны;
  * консервативное: тот же сдвиг, но расширение метки только до наблюдённого расширения
    предсказаний, то есть нижняя граница. Эффективный размер 492-1355, считается честно.

Итог считается на полном ST-RAE, а не на одном знаменателе: числитель под сдвинутым
распределением тоже растёт, и расширение возвращает лишь около четверти потери."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import numpy as np, pandas as pd, json
from sklearn.isotonic import IsotonicRegression
from scipy.optimize import fsolve, brentq
from scipy.stats import ks_2samp

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
rows = pd.read_csv(D + "rows.csv")
tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
        .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
oof = json.load(open(RES + "preds/oof.json"))
_pp = _pl.Path(D) / "test_pred.npz"
if not _pp.exists():
    raise SystemExit("нет data/test_pred.npz — сначала verify/k5_shift.py")
PT = np.load(_pp)


def W(y, t):
    z = t[0] * (y - y.mean()) + t[1] * (y - y.mean()) ** 2; z -= z.max()
    w = np.exp(z); return w / w.sum()


def S(y, p, lo, hi, w):
    yb = w @ y
    return (w @ (np.maximum(p - hi, 0) + np.maximum(lo - p, 0))) / \
           (w @ (np.maximum(yb - hi, 0) + np.maximum(lo - yb, 0)))


def pack(c):
    col = f"{c}_pIC50_direct_inhibition"; m = tr[col].notna().to_numpy()
    return (tr.loc[m, col].to_numpy(), tr.loc[m, col + "_conf_low"].to_numpy(),
            tr.loc[m, col + "_conf_high"].to_numpy(),
            np.asarray(oof[f"FP+DESC+MECH|{c}"]), PT[c])


print("=" * 104)
print("1. Подгонка двух моментов вместо одного, и где она вырождается")
print("=" * 104)
print(f"{'фермент':8s} {'d 1 мом.':>9s} {'d 2 мом.':>9s} {'sd(y) обуч':>10s} {'sd(y) 2м':>9s} "
      f"{'x':>5s} | {'эфф.n 1м':>8s} {'эфф.n 2м':>8s} {'из':>6s} {'годно?':>8s}")
FIT = {}
for c in CYPS:
    y, lo, hi, p, pt = pack(c)
    mv = IsotonicRegression(out_of_bounds="clip").fit(y, p).predict(y)
    s2 = np.var(p - mv)
    t1 = brentq(lambda a: W(y, [a, 0]) @ mv - pt.mean(), -80, 80, xtol=1e-10)
    w1 = W(y, [t1, 0]); d1 = w1 @ y - y.mean()
    t2 = fsolve(lambda t: [W(y, t) @ mv - pt.mean(),
                           W(y, t) @ (mv - W(y, t) @ mv) ** 2 + s2 - pt.var()], [0.0, 0.0])
    w2 = W(y, t2); d2 = w2 @ y - y.mean()
    sd2 = np.sqrt(w2 @ (y - d2 - y.mean()) ** 2)
    n1, n2 = 1 / np.sum(w1 ** 2), 1 / np.sum(w2 ** 2)
    ok = "да" if n2 > 0.1 * len(y) else "НЕТ"
    ratio = pt.std() / p.std()
    tc = fsolve(lambda t: [W(y, t) @ y - (y.mean() + d1),
                           np.sqrt(W(y, t) @ (y - W(y, t) @ y) ** 2) - y.std() * ratio], [0.0, 0.0])
    FIT[c] = (d1, w1, W(y, tc), 1 / np.sum(W(y, tc) ** 2), ratio)
    print(f"{c:8s} {d1:+9.3f} {d2:+9.3f} {y.std():10.3f} {sd2:9.3f} {sd2/y.std():5.2f} | "
          f"{n1:8.0f} {n2:8.0f} {len(y):6d} {ok:>8s}")
print("\nСтрока с 'НЕТ' — те же грабли, что в k5: подгонка держится на двух десятках молекул.")
print("Её числа выброшены; дальше только консервативная версия.")

print()
print("=" * 104)
print("2. Консервативное расширение: метка расширена ровно во столько, во сколько предсказания")
print("=" * 104)
print(f"{'фермент':8s} {'delta':>7s} {'sd(yhat) т/о':>12s} {'эфф.n':>7s} {'/n':>6s}")
for c in CYPS:
    y = pack(c)[0]
    d1, w1, wc, nc, ratio = FIT[c]
    print(f"{c:8s} {d1:+7.3f} {ratio:12.3f} {nc:7.0f} {nc/len(y):6.2f}")

print()
print("=" * 104)
print("3. Полный ST-RAE, а не один знаменатель: числитель тоже растёт")
print("=" * 104)
print(f"{'фермент':8s} {'delta = 0':>10s} {'только сдвиг':>13s} {'сдвиг + расширение':>19s}")
A = B = C = 0.0
for c in CYPS:
    y, lo, hi, p, pt = pack(c)
    d1, w1, wc, nc, ratio = FIT[c]
    u = np.ones(len(y)) / len(y)
    a, b, cc = S(y, p, lo, hi, u), S(y, p, lo, hi, w1), S(y, p, lo, hi, wc)
    A += a / 4; B += b / 4; C += cc / 4
    print(f"{c:8s} {a:10.4f} {b:13.4f} {cc:19.4f}")
print(f"{'МАКРО':8s} {A:10.4f} {B:13.4f} {C:19.4f}")
print(f"\nРасширение возвращает {(B-C)/(B-A)*100:.0f}% потери от сдвига, а не отменяет её.")
print(f"Ожидание на промежуточном лидерборде — порядка {C:.2f} макро, а не {A:.2f}.")
print("Проверяется одним числом 24 сентября.")
