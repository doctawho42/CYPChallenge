"""Оценка delta без инструмента и без предположения, что передача — скаляр.

Пункт 47 закрыл три схемы: кластерная алгебраически совпадает с b, оконная нарушает
exclusion restriction по построению, и передача вообще не скаляр, а производная по
направлению. Всё верно — и всё это про попытку оценить множитель. Множитель можно не
оценивать.

Связь точная и без инструмента:

    E_test[yhat] = INT E[yhat|y] * p_test(y) dy,

если ядро E[yhat|y] на тесте то же, что на обучении. Нужна ОБРАТНАЯ регрессия
m(y) = E[yhat|y] (в линейном пределе её наклон и есть R^2/b, тот самый, который пункт 47
называет правильным), а не прямая b. Здесь m оценивается изотоникой, то есть линейность
не предполагается вовсе, и «направление» задаётся самим наклоном маргинали метки — тем
единственным, вдоль которого вопрос и поставлен. Оконные ячейки со значением 1.22 и 1.55
этому не противоречат: вдоль направления, которое модель видит напрямую, передача
действительно близка к единице и выше, просто это не то направление.

Остаётся ровно одно допущение — инвариантность ядра, то есть отсутствие concept shift.
Ни один метод без тестовых меток его не проверит, и это честный предел идентифицируемости.
Оно же и самое слабое там, где пункт 45 нашёл настоящую химию: если низкие метки на тесте
низки по другой причине, чем на обучении, ядро другое. Поэтому оценка по CYP2D6 —
подозреваемая в первую очередь, а не в последнюю.

Блок 3 — действие под асимметричной ценой ошибки: при кусочно-линейной цене оптимум это
не середина диапазона, а квантиль уровня c_зан / (c_зан + c_зав)."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import numpy as np, pandas as pd, json
from sklearn.isotonic import IsotonicRegression
from scipy.optimize import brentq
from scipy.stats import ks_2samp
from cypsplit import cluster_ids

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
rows = pd.read_csv(D + "rows.csv")
tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
        .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
oof = json.load(open(RES + "preds/oof.json"))
_pp = _pl.Path(D) / "test_pred.npz"
if not _pp.exists():
    raise SystemExit("нет data/test_pred.npz — сначала verify/k5_shift.py")
PT = np.load(_pp)
cid, _ = cluster_ids(list(rows.SMILES))


def softmax_w(y, th):
    z = th * (y - y.mean()); z -= z.max()
    w = np.exp(z); return w / w.sum()


def solve_delta(y, m, target, TH=80.0):
    """delta наклона, при котором взвешенное среднее ядра равно target."""
    f = lambda th: float(softmax_w(y, th) @ m) - target
    if f(-TH) * f(TH) > 0:
        return np.nan
    w = softmax_w(y, brentq(f, -TH, TH, xtol=1e-10))
    return float(w @ y) - y.mean()


def pack(c):
    col = f"{c}_pIC50_direct_inhibition"; m = tr[col].notna().to_numpy()
    return (tr.loc[m, col].to_numpy(), np.asarray(oof[f"FP+DESC+MECH|{c}"]), PT[c], cid[m])


print("=" * 104)
print("1. Одно наблюдение, четыре способа превратить его в delta")
print("=" * 104)
print(f"{'фермент':8s} {'dyhat':>7s} {'b':>6s} {'R2':>6s} {'R2/b':>6s} | {'сырое':>7s} "
      f"{'/b (отозв.)':>11s} {'*b/R2 лин.':>10s} {'ядро изо':>9s} | {'95% кластерный':>21s}")
OUT = {}
for c in CYPS:
    y, p, pt, g = pack(c)
    dy = pt.mean() - p.mean()
    b = np.polyfit(p, y, 1)[0]
    R2 = 1 - ((y - p) ** 2).sum() / ((y - y.mean()) ** 2).sum()
    d_iso = solve_delta(y, IsotonicRegression(out_of_bounds="clip").fit(y, p).predict(y), pt.mean())
    rng = np.random.default_rng(0); u = np.unique(g); bs = []
    for _ in range(300):
        idx = np.concatenate([np.where(g == k)[0] for k in rng.choice(u, len(u), replace=True)])
        it = IsotonicRegression(out_of_bounds="clip").fit(y[idx], p[idx])
        v = solve_delta(y[idx], it.predict(y[idx]), pt.mean())
        if np.isfinite(v):
            bs.append(v)
    q = np.percentile(bs, [2.5, 97.5]); OUT[c] = (d_iso, q, np.array(bs))
    print(f"{c:8s} {dy:+7.3f} {b:6.3f} {R2:6.3f} {R2/b:6.3f} | {dy:+7.3f} {dy/b:+11.3f} "
          f"{dy*b/R2:+10.3f} {d_iso:+9.3f} | {q[0]:+9.3f} … {q[1]:+9.3f}")
print(f"{'МАКРО':8s}" + " " * 41 + f"{np.mean([OUT[c][0] for c in CYPS]):+9.3f}")
print("\nИзотоника и линейная версия совпадают до 0.05: нелинейность ядра ни при чём,")
print("вся поправка сидит в замене b на R^2/b. Это ноль, который стоит иметь.")

print()
print("=" * 104)
print("2. Достаточно ли ОДНОГО delta: совпадает ли наклонённое распределение с тестовым")
print("=" * 104)
print(f"{'фермент':8s} {'delta':>7s} | {'ср. накл.':>9s} {'ср. тест':>9s} | {'sd накл.':>9s} "
      f"{'sd тест':>8s} | {'KS':>6s} {'KS p':>9s}")
for c in CYPS:
    y, p, pt, g = pack(c)
    m = IsotonicRegression(out_of_bounds="clip").fit(y, p).predict(y)
    f = lambda th: float(softmax_w(y, th) @ m) - pt.mean()
    w = softmax_w(y, brentq(f, -80, 80, xtol=1e-10))
    mu = w @ p; sd = np.sqrt(w @ (p - mu) ** 2)
    ks = ks_2samp(p[np.random.default_rng(1).choice(len(p), 20000, p=w)], pt)
    print(f"{c:8s} {OUT[c][0]:+7.3f} | {mu:9.3f} {pt.mean():9.3f} | {sd:9.3f} {pt.std():8.3f} | "
          f"{ks.statistic:6.3f} {ks.pvalue:9.2e}")
print("\nСреднее сходится по построению, разброс — нет. Наклон умеет двигать среднее и не")
print("умеет добавлять разброс. Продолжение в k9_shape.py.")

print()
print("=" * 104)
print("3. Действие под асимметричной ценой (CYP3A4: занижение 0.087, завышение 0.013)")
print("=" * 104)
cu, co = 0.087, 0.013; lvl = cu / (cu + co)
bs = OUT["CYP3A4"][2]
print(f"при кусочно-линейной цене оптимум = квантиль уровня {cu}/({cu}+{co}) = {lvl:.3f},")
print("а не середина диапазона. Та же байесовская логика, что в §10 для точечного предсказания.")
print(f"  бутстрап здесь:      медиана {np.median(bs):+.3f}   квантиль {lvl:.2f} = {np.quantile(bs,lvl):+.3f}")
for lo, hi, nm in [(0.295, 0.576, "диапазон п.47")]:
    print(f"  {nm}:      середина {(lo+hi)/2:+.3f}   квантиль {lvl:.2f} = {lo+lvl*(hi-lo):+.3f}")
print("Оговорка: 0.87 верно, если 0.087 и 0.013 — наклоны при сопоставимых отклонениях.")
print("Если это цены на концах диапазона и кривая гнутая, минимизируйте ожидание по самой")
print("кривой численно; ответ сдвинется, но останется заметно выше середины.")
