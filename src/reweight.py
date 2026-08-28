"""What would we score if the test's activity mixture were ours?

The similarity geometry of the test set cannot be reproduced by any re-split of the
training data - leave-one-out tops out at median nearest-neighbour 0.450 and the test
sits at 0.587, above that ceiling. The label marginal is a different matter: it can be
reproduced by reweighting, and the label marginal is what actually breaks post-hoc
shrinkage, so it is the half worth having.

Method. Exponential tilting: w_i proportional to exp(theta * y_i), with theta chosen per
enzyme so that the weighted mean activity sits delta above the training mean. Tilting is
the minimum-relative-entropy way to move a mean, so the reweighted sample is the closest
distribution to ours that has the mean we asked for, rather than an arbitrary one.

ST-RAE is a ratio of sums, so both halves are weighted - and the denominator's constant
predictor becomes the *weighted* mean, because under the target distribution that is what
predicting one number for everything would give.

The output is a curve, not a point. We do not know the test's mean shift; we know the
anchors sit at the 92nd to 99th percentile of our own distribution per enzyme, and that
the test is those anchors plus their analogues. So the honest deliverable is "score as a
function of delta", with the plausible range marked, and the reader can see how fast the
answer moves.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()

import json

import numpy as np
import pandas as pd
from scipy.optimize import brentq

from cypsplit import butina_folds

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
GRID = np.linspace(0.2, 1.0, 41)


def soft(pred, lo, hi):
    """Per-compound ST-RAE numerator: distance outside the tolerance band."""
    return np.maximum(np.maximum(pred - hi, 0.0), np.maximum(lo - pred, 0.0))


def tilt(y, delta):
    """Weights proportional to exp(theta*y) whose weighted mean is mean(y) + delta."""
    if abs(delta) < 1e-9:
        return np.ones_like(y)
    target = y.mean() + delta

    def gap(th):
        w = np.exp(th * (y - y.mean()))
        return np.average(y, weights=w) - target

    # The reachable range is bounded by max(y); refuse rather than silently clip.
    hi_th = 50.0
    if gap(hi_th) < 0:
        raise ValueError(f"сдвиг {delta} недостижим: максимум метки {y.max():.2f}")
    th = brentq(gap, -hi_th, hi_th, xtol=1e-10)
    w = np.exp(th * (y - y.mean()))
    return w / w.mean()


def wstrae(y, pred, lo, hi, w):
    """ST-RAE under the reweighted marginal, constant predictor included."""
    mu = np.average(y, weights=w)
    return float((w * soft(pred, lo, hi)).sum() / (w * soft(np.full_like(y, mu), lo, hi)).sum())


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    oof = json.load(open(RES + "preds/oof.json"))
    fold, _ = butina_folds(list(rows.SMILES))

    def arm(c, shift):
        """Raw predictions and the shrunk ones, shrinkage fitted out of fold as usual."""
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy()
        lo = tr.loc[m, col + "_conf_low"].to_numpy()
        hi = tr.loc[m, col + "_conf_high"].to_numpy()
        p = np.array(oof[f"FP+DESC+MECH|{c}"])
        fi = fold[m]
        o = np.zeros_like(p)
        for f in range(5):
            te, trn = fi == f, fi != f
            if te.sum() == 0:
                continue
            mu = p[trn].mean() + shift
            L = min(GRID, key=lambda L: wstrae(y[trn], mu + L * (p[trn] - mu),
                                               lo[trn], hi[trn], np.ones(trn.sum())))
            o[te] = mu + L * (p[te] - mu)
        return y, lo, hi, p, o

    print("Сдвиг маргинали меток: макро ST-RAE при разных вариантах постобработки.")
    print("delta - насколько среднее активности на тесте выше обучающего.\n")
    print(f"{'delta':>7s} {'сырое':>9s} {'усадка(0)':>11s} {'усадка(+0.4)':>13s} {'лучший сдвиг':>13s}")

    cache = {sh: {c: arm(c, sh) for c in CYPS} for sh in (0.0, 0.4)}
    shifts = np.arange(0.0, 1.05, 0.1)
    for delta in np.arange(0.0, 0.85, 0.1):
        raw, s0, s4 = [], [], []
        for c in CYPS:
            y, lo, hi, p, o0 = cache[0.0][c]
            _, _, _, _, o4 = cache[0.4][c]
            w = tilt(y, delta)
            raw.append(wstrae(y, p, lo, hi, w))
            s0.append(wstrae(y, o0, lo, hi, w))
            s4.append(wstrae(y, o4, lo, hi, w))
        # which shrinkage centre would have been best at this delta
        best = None
        for sh in shifts:
            v = []
            for c in CYPS:
                y, lo, hi, p, o = arm(c, sh)
                v.append(wstrae(y, o, lo, hi, tilt(y, delta)))
            if best is None or np.mean(v) < best[1]:
                best = (sh, float(np.mean(v)))
        print(f"{delta:7.1f} {np.mean(raw):9.4f} {np.mean(s0):11.4f} {np.mean(s4):13.4f} "
              f"{best[0]:+7.1f} -> {best[1]:.4f}")

    print("""
Читать так. Правый столбец - главное: центр усадки, который был бы оптимален, если бы
тестовая маргиналь была сдвинута на delta. Он растёт вместе с delta, и это ровно то, за
что нельзя было ухватиться на обучающей маргинали: подбирая центр по своим данным, мы
подбираем его под delta = 0, а тест сидит правее.""")

    # Где на этой кривой сидит настоящий тест. Перцентили якорей измерены f7_testset и
    # устойчивы к порогу кластеризации от 0.35 до 0.50.
    anch = {"CYP1A2": 93.0, "CYP2C9": 98.0, "CYP2D6": 62.0, "CYP3A4": 90.0}
    print("\n=== где на этой кривой сидит тест ===")
    print(f"{'фермент':8s} {'среднее':>8s} {'перцентиль':>11s} {'сдвиг якоря':>12s}")
    ds = []
    for c in CYPS:
        yy = tr[f"{c}_pIC50_direct_inhibition"].dropna().to_numpy()
        d = float(np.percentile(yy, anch[c]) - yy.mean())
        ds.append(d)
        print(f"{c:8s} {yy.mean():8.3f} {anch[c]:11.0f} {d:+12.3f}")

    print(f"""
Внутренний контроль сошёлся: CYP2D6 в отборе якорей не участвовал, и сдвиг у него
{ds[2]:+.2f}, то есть практически ноль, тогда как три отборных фермента дают
{ds[0]:+.2f}, {ds[1]:+.2f} и {ds[3]:+.2f}. Метод отличает обогащённый набор от необогащённого.

Сдвиг ЯКОРЕЙ в среднем {np.mean(ds):+.2f}. Тест - это якоря плюс по десять соседей, а соседей
отбирали по сходству, не по активности, так что они регрессируют к среднему и тянут
маргиналь вниз. Значит {np.mean(ds):.2f} - верхняя граница, а не оценка. Разумный диапазон
delta примерно от 0.3 до 0.6, и на таблице выше это центр усадки около +0.8, а не +0.4,
который получается подбором на своей маргинали.

Второе, что говорит эта таблица, неприятнее и к усадке отношения не имеет. Наш заголовочный
0.767 посчитан при delta = 0. Если тест обогащён активными, то при том же качестве
предсказаний ST-RAE вырастет: у более узкой по активности выборки меньше знаменатель.
Это тот же механизм, из-за которого CYP2D6 стоит на 0.98.

Заметим, что две оси тянут в разные стороны. По СХОДСТВУ тест легче нашей кросс-валидации,
и на верхней четверти по сходству макро идёт 0.767 -> 0.739 в нашу пользу. По МАРГИНАЛИ
меток тест тяжелее. Какая ось перевесит, до промежуточного лидерборда не узнать, и до сих
пор мы отчитывались только по первой.""")


if __name__ == "__main__":
    main()
