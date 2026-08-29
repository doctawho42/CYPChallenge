"""The kernel estimate of delta, split along the one axis known to have moved.

What this fixes. verify/k8_kernel.py infers the label shift by tilting our own label marginal
until the mean of the kernel E[yhat|y] matches the mean prediction observed on the test. It
needs one assumption and only one: that the kernel is the same function on the test as on the
training set. That assumption is weakest exactly where the answer matters most.

verify/k7_2d6shift.py showed that the test's composition moved in a chemically specific
direction: it carries about a third as many compounds that are basic at pH 7.4 as the CYP2D6
label mask does, 0.104 against 0.355. CYP2D6 is the one enzyme where basic compounds are MORE
active, because it binds through a salt bridge to a protonated nitrogen. So on that enzyme the
test differs from the training set along an axis that is directly relevant to what is being
predicted, and a single kernel fitted on the pooled training data cannot know which side of
that axis a test compound sits on.

The consequence is not subtle. An unstratified kernel attributes the whole observed drop in
mean prediction to a drop in labels. But part of that drop is composition: fewer of the
compounds that score high on this enzyme. Charging composition to the label shift inflates
|delta|, and CYP2D6 is where our chosen shift is largest and least trustworthy.

The repair replaces the assumption with a measurement. Estimate the kernel separately inside
and outside the basic stratum, then mix the two by the test's OWN measured composition rather
than the training set's:

    E_test[yhat] = f_test * Int m_basic(y) p_basic(y) dy
                 + (1 - f_test) * Int m_other(y) p_other(y) dy

and solve for the tilt that makes it match. What was an assumption of invariance becomes a
measured mixing weight along the single axis known to have moved. Everything else still rests
on invariance, so this narrows the assumption rather than removing it.

Two readings come out of it separately, and the split is the point:

  the composition part, obtained at zero tilt - what the shift would be if the labels within
  each stratum were unchanged and only the mixture moved;

  the residual tilt - what is left over and must be a genuine change in the labels.

The other three enzymes serve as the control. Their masks are 0.114 to 0.122 basic against the
test's 0.104, a difference of about one point, so stratification should barely move them. If
it moves them a lot, the split is doing something other than what it claims.

Reads results/preds/oof.json, data/test_pred.npz, data/test_feats.npz, data/feats.npz.
data/test_* come from verify/k5_shift.py. Writes nothing.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()

import json

import numpy as np
import pandas as pd
from scipy.optimize import brentq
from sklearn.isotonic import IsotonicRegression

from cypsplit import cluster_ids

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
TH = 80.0


def softmax_w(y, th):
    """Exponential tilt as weights, in the same form verify/k8_kernel.py uses."""
    z = th * (y - y.mean())
    z -= z.max()
    w = np.exp(z)
    return w / w.sum()


def iso(y, p):
    return IsotonicRegression(out_of_bounds="clip").fit(y, p).predict(y)


def solve_pooled(y, m, target):
    f = lambda th: float(softmax_w(y, th) @ m) - target
    if f(-TH) * f(TH) > 0:
        return np.nan
    w = softmax_w(y, brentq(f, -TH, TH, xtol=1e-10))
    return float(w @ y) - y.mean()


def solve_strat(yb, mb, yo, mo, f_test, y_all_mean, target):
    """One tilt applied inside both strata, mixed by the test's measured composition.

    The tilt is shared: the assumption is that the labels move by the same amount within
    each stratum, and that what differs between the sets is how much of each stratum there
    is. Returns the total shift of the label mean, composition included.
    """
    def pred(th):
        return (f_test * float(softmax_w(yb, th) @ mb)
                + (1 - f_test) * float(softmax_w(yo, th) @ mo))

    g = lambda th: pred(th) - target
    if g(-TH) * g(TH) > 0:
        return np.nan, np.nan
    th = brentq(g, -TH, TH, xtol=1e-10)
    mix = (f_test * float(softmax_w(yb, th) @ yb)
           + (1 - f_test) * float(softmax_w(yo, th) @ yo))
    comp = f_test * yb.mean() + (1 - f_test) * yo.mean() - y_all_mean
    return mix - y_all_mean, comp


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    i = mn.index("is_base_74")
    base = np.load(D + "feats.npz")["MECH"][:, i].astype(float) > 0.5
    baset = np.load(D + "test_feats.npz")["MECH"][:, i].astype(float) > 0.5
    f_test = float(baset.mean())
    oof = json.load(open(RES + "preds/oof.json"))
    PT = np.load(D + "test_pred.npz")
    cid, _ = cluster_ids(list(rows.SMILES))
    cid = np.asarray(cid)
    # Тестовая сторона тоже ресэмплится, и тоже блоками. Первая версия этого скрипта
    # держала pt.mean() фиксированной внутри бутстрапа: она ловила неопределённость ядра и
    # выбрасывала неопределённость цели, а цель входит в обращение с множителем 1/наклон и
    # весит больше. Дефект унаследован от verify/k8_kernel.py и там же исправлен.
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    tcid, _ = cluster_ids(list(te.SMILES), threshold=0.50)
    tcid = np.asarray(tcid)
    tu = np.unique(tcid)

    print("=" * 100)
    print("1. Состав страт: насколько маска фермента и тест отличаются по основаниям")
    print("=" * 100)
    print(f"доля оснований на тесте: {f_test:.3f}\n")
    print(f"{'фермент':8s} {'в маске':>8s} {'n основ.':>9s} {'n прочих':>9s} "
          f"{'акт. основ.':>12s} {'акт. прочих':>12s} {'разрыв':>8s}")
    S = {}
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy()
        p = np.asarray(oof[f"FP+DESC+MECH|{c}"], float)
        b = base[m]
        S[c] = (y, p, b, cid[m], PT[c])
        print(f"{c:8s} {b.mean():8.3f} {int(b.sum()):9d} {int((~b).sum()):9d} "
              f"{y[b].mean():12.3f} {y[~b].mean():12.3f} "
              f"{y[b].mean() - y[~b].mean():+8.3f}")

    print()
    print("=" * 100)
    print("2. Различаются ли сами ядра внутри страт")
    print("=" * 100)
    print("   ядро E[yhat|y] в трёх точках шкалы; если страты дают одну функцию, стратификация")
    print("   ничего не изменит по построению\n")
    print(f"{'фермент':8s} " + " ".join(f"{('y=' + str(q)):>16s}" for q in (3.5, 4.5, 5.5)))
    for c in CYPS:
        y, p, b, g, pt = S[c]
        ib = IsotonicRegression(out_of_bounds="clip").fit(y[b], p[b])
        io = IsotonicRegression(out_of_bounds="clip").fit(y[~b], p[~b])
        cells = []
        for q in (3.5, 4.5, 5.5):
            cells.append(f"{float(ib.predict([q])[0]):6.3f}/{float(io.predict([q])[0]):<6.3f}")
        print(f"{c:8s} " + " ".join(f"{x:>16s}" for x in cells))
    print("\n   формат: основания / прочие")

    print()
    print("=" * 100)
    print("3. Delta: объединённое ядро против стратифицированного")
    print("=" * 100)
    print(f"{'фермент':8s} {'объединённое':>13s} {'стратифиц.':>11s} {'из них состав':>14s} "
          f"{'из них метки':>13s} {'разница':>9s}")
    res = {}
    for c in CYPS:
        y, p, b, g, pt = S[c]
        d_pool = solve_pooled(y, iso(y, p), pt.mean())
        d_str, d_comp = solve_strat(y[b], iso(y[b], p[b]), y[~b], iso(y[~b], p[~b]),
                                    f_test, y.mean(), pt.mean())
        res[c] = (d_pool, d_str, d_comp)
        print(f"{c:8s} {d_pool:+13.3f} {d_str:+11.3f} {d_comp:+14.3f} "
              f"{d_str - d_comp:+13.3f} {d_str - d_pool:+9.3f}")

    print()
    print("=" * 100)
    print("4. Кластерный бутстрап стратифицированной оценки")
    print("=" * 100)
    rng = np.random.default_rng(0)
    draws = {}
    print("   Ресэмплятся ОБЕ стороны: обучающие кластеры и тестовые серии.\n")
    print(f"{'фермент':8s} {'стратифиц.':>11s} {'95% интервал':>22s} {'ноль внутри':>12s}")
    for c in CYPS:
        y, p, b, g, pt = S[c]
        u = np.unique(g)
        bs = []
        for _ in range(1500):
            idx = np.concatenate([np.where(g == k)[0]
                                  for k in rng.choice(u, len(u), replace=True)])
            yb, ob = y[idx][b[idx]], p[idx][b[idx]]
            yo, oo = y[idx][~b[idx]], p[idx][~b[idx]]
            if len(yb) < 40 or len(yo) < 40:
                continue
            tid = np.concatenate([np.where(tcid == k)[0]
                                  for k in rng.choice(tu, len(tu), replace=True)])
            v, _ = solve_strat(yb, iso(yb, ob), yo, iso(yo, oo),
                               f_test, y[idx].mean(), float(pt[tid].mean()))
            if np.isfinite(v):
                bs.append(v)
        lo, hi = np.percentile(bs, [2.5, 97.5])
        draws[c] = [float(x) for x in bs]
        print(f"{c:8s} {res[c][1]:+11.3f} [{lo:+8.3f}, {hi:+8.3f}] "
              f"{('да' if lo * hi < 0 else 'нет'):>12s} "
              f"  P(delta>=0) = {float((np.asarray(bs) >= 0).mean()):.3f}")
    # Розыгрыши сохраняются: имея их, можно сравнивать правила по СРЕДНЕМУ, а не только по
    # худшему случаю внутри диапазона. Это разные цели, и выбор между ними надо делать
    # явно, а не наследовать (src/shrinkchoice.py, блок 6).
    json.dump(draws, open(RES + "preds/delta_draws.json", "w"))
    print(f"\n  розыгрыши сохранены: {RES}preds/delta_draws.json")

    print("""
Как это читать.

Столбец «из них состав» --- это сдвиг, который получился бы, если бы метки внутри каждой
страты не менялись вовсе, а поменялась только пропорция страт. Он вычисляется без всякой
подгонки: доли известны, средние по стратам известны. Столбец «из них метки» --- остаток,
который приписать составу нельзя.

На трёх ферментах, где маска и тест по основаниям почти не различаются, стратификация обязана
почти ничего не менять, и это контроль конструкции, а не результат. Смотреть надо на CYP2D6,
где доли расходятся втрое.

Чего это по-прежнему не снимает. Внутри каждой страты инвариантность ядра всё ещё
предполагается. Заменено одно допущение на более узкое: вместо «состав теста не важен» теперь
«состав важен ровно вдоль оси протонирования, и вдоль неё он измерен». Если тест смещён и по
другой оси, релевантной для 2D6, эта поправка её не поймает.""")


if __name__ == "__main__":
    main()
