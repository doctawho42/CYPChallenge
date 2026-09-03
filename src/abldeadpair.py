"""Are the dead zone and the pairwise objective the same idea twice, or two different ones.

Why this has to be asked before either goes into the submission. Both are attempts to stop
spending capacity on what the metric will not pay for, and they proceed from the same
observation from opposite ends:

    мёртвая зона   ошибка ВНУТРИ полосы бесплатна, значит не подгоняйся туда  (поточечно)
    попарно        пару с пересекающимися полосами метрика не упорядочит      (попарно)

Both discard information carried by the band widths. If they discard the *same* information, then
stacking them buys nothing and the ensemble should carry one of them; if different, the effects
add and both belong. Item 137 already noted that the two sets are not identical -- a pair of
wide-banded compounds can still be strictly ordered, and item 115's Jaccard overlap is lowest on
CYP2C9 at 0.738 -- but "not identical" is not the same as "additive", and nothing has measured
which.

Measured separately so far:

    мёртвая зона x1, HistGB   +0.027 ранга на сиде 0, +0.032 на сиде 1
    попарно, свой бустинг     +0.008 ранга по четырём сидам против своего же квадрата

The design is a two-by-two, and it is the only shape that answers the question. Everything runs
on the same learner -- the plain-tree booster pinned to `src/ablpairloss.py` -- because the
dead-zone number above came from HistGB and the pairwise number from this booster, so as they
stand they are not comparable to each other at all.

    квадрат                     the reference corner
    попарно                     pairwise alone
    квадрат + мёртвая зона x1   dead zone alone
    попарно + мёртвая зона x1   both

Pre-registered reading, and it is a subtraction that cannot be argued about afterwards:

    additive   ->  (попарно+МЗ - попарно)  ~=  (квадрат+МЗ - квадрат),  both belong
    overlapping->  попарно+МЗ  ~=  max(попарно, квадрат+МЗ),            carry one

The dead-zone pass is the majorise-minimise step of `src/abldead.py`, unchanged: project the
honest out-of-fold predictions onto the band with clip(p, lo, hi) and refit against that target
under L1. It is applied on top of whichever corner it is stacked on, so the pairwise corner is
projected from the pairwise predictions and not from the squared ones. One pass only -- item on
`abldead` measured passes 2 and 3 degrading (0.6931 -> 0.6951 -> 0.6999).

L1 boosting over plain trees is the sign-residual scheme: fitting each tree to sign(y - s) is
gradient descent on the absolute loss, exactly as fitting to (y - s) is on the squared one.

Same folds, masks and metric as everywhere. Writes results/preds/oof_deadpair.json. ~20 min/seed.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression
from sklearn.tree import DecisionTreeRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
NTREE, LR, DEPTH, MAXFEAT = 200, 0.06, 5, 0.3     # пины из src/ablpairloss.py


def pairwise_grad(s, y, lo, hi, chunk=512):
    """Градиент попарной логистической потери по РАЗЛИЧИМЫМ парам. Копия из ablpairloss."""
    n = len(s)
    g = np.zeros(n)
    for a in range(0, n, chunk):
        b = min(a + chunk, n)
        above = lo[a:b, None] > hi[None, :]
        if not above.any():
            continue
        d = s[a:b, None] - s[None, :]
        w = np.where(above, 1.0 / (1.0 + np.exp(np.clip(d, -30, 30))), 0.0)
        g[a:b] -= w.sum(axis=1)
        g += w.sum(axis=0)
    return g


def boost(Xtr, ytr, lotr, hitr, Xte, mode, seed):
    """mode: 'sq' квадрат, 'l1' абсолютная по метке, 'pair' попарная."""
    rng = np.random.default_rng(seed)
    s = np.full(len(ytr), float(np.mean(ytr)) if mode != "pair" else 0.0)
    pred = np.full(len(Xte), s[0])
    for _ in range(NTREE):
        if mode == "sq":
            resid = ytr - s
        elif mode == "l1":
            resid = ytr - s
        else:
            resid = -pairwise_grad(s, ytr, lotr, hitr)
            nrm = np.abs(resid).max()
            if nrm > 0:
                resid = resid / nrm
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MAXFEAT,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(Xtr, resid)
        if mode == "l1":
            # LAD-бустинг Фридмана: дерево строится по остаткам, но значение листа
            # заменяется МЕДИАНОЙ остатков в нём --- это и есть шаг под абсолютную
            # потерю. Без этой замены получается бустинг по знаку, другой оценщик,
            # и сравнивать его с HistGB(loss="absolute_error") нельзя.
            leaf = t.apply(Xtr)
            for lv in np.unique(leaf):
                t.tree_.value[lv, 0, 0] = np.median(resid[leaf == lv])
        s = s + LR * t.predict(Xtr)
        pred = pred + LR * t.predict(Xte)
    return s, pred


def cv(X, y, lo, hi, fi, mode, seed, target=None):
    """Один полный кросс-валидационный проход. target отличается от метки в шаге МЗ."""
    tgt = y if target is None else target
    p = np.zeros(len(y))
    for f in range(5):
        trn, te = fi != f, fi == f
        if te.sum() == 0:
            continue
        s_tr, s_te = boost(X[trn], tgt[trn], lo[trn], hi[trn], X[te], mode, seed * 10 + f)
        if mode == "pair":
            # Попарная потеря задаёт порядок и молчит про шкалу.
            p[te] = IsotonicRegression(out_of_bounds="clip").fit(s_tr, tgt[trn]).predict(s_te)
        else:
            p[te] = s_te
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--out", default=RES + "preds/oof_deadpair.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)

    # «L1 по метке» --- обязательный контроль: он отделяет «мёртвая зона не работает на
    # этом учителе» от «мой L1 плох». Без него разности внизу нечитаемы (пункт про
    # отозванный прогон).
    NAMES = ["квадрат", "L1 по метке", "попарно", "квадрат + МЗ x1", "попарно + МЗ x1"]
    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        acc = {nm: {"seed": seed, "рука": nm} for nm in NAMES}
        clock = {nm: 0.0 for nm in NAMES}

        for c in CYPS:
            col = f"{c}_pIC50_direct_inhibition"
            m = tr[col].notna().to_numpy()
            y = tr.loc[m, col].to_numpy()
            lo = tr.loc[m, col + "_conf_low"].to_numpy()
            hi = tr.loc[m, col + "_conf_high"].to_numpy()
            Xi, fi = X[m], fold[m]
            u = np.ones(len(y)) / len(y)

            preds = {}
            t0 = time.time(); preds["квадрат"] = cv(Xi, y, lo, hi, fi, "sq", seed)
            clock["квадрат"] += time.time() - t0
            t0 = time.time(); preds["L1 по метке"] = cv(Xi, y, lo, hi, fi, "l1", seed)
            clock["L1 по метке"] += time.time() - t0
            t0 = time.time(); preds["попарно"] = cv(Xi, y, lo, hi, fi, "pair", seed)
            clock["попарно"] += time.time() - t0
            # Шаг мажорирования поверх каждого угла: проекция ЧЕСТНЫХ OOF на полосу,
            # затем переподгонка под L1 против неё.
            for src, dst in (("квадрат", "квадрат + МЗ x1"), ("попарно", "попарно + МЗ x1")):
                t0 = time.time()
                preds[dst] = cv(Xi, y, lo, hi, fi, "l1", seed,
                                target=np.clip(preds[src], lo, hi))
                clock[dst] += time.time() - t0

            for nm, p in preds.items():
                q = fit_apply(p, lo, hi, fi, u)
                out[f"{seed}|{nm}|{c}"] = p.tolist()
                acc[nm][f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi,
                                                         y_true_lower=lo)), 4)
                acc[nm][f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)

        for nm in NAMES:
            r = acc[nm]
            for tag in ("пара", "rho"):
                r[f"MACRO {tag}"] = round(float(np.mean([r[f"{c} {tag}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {nm:18s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({clock[nm]:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    g = df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                       + [f"{c} rho" for c in CYPS]].mean()
    print()
    print(g.round(4).to_string())

    sq, pr = g.loc["квадрат", "MACRO rho"], g.loc["попарно", "MACRO rho"]
    sqd, prd = g.loc["квадрат + МЗ x1", "MACRO rho"], g.loc["попарно + МЗ x1", "MACRO rho"]
    print(f"\nэффект МЗ на квадрате : {sqd - sq:+.4f}")
    print(f"эффект МЗ на попарном : {prd - pr:+.4f}")
    print(f"эффект попарного      : {pr - sq:+.4f}")
    print(f"сложение предсказывало: {sq + (sqd - sq) + (pr - sq):.4f}, "
          f"намерено {prd:.4f}, недобор {sq + (sqd - sq) + (pr - sq) - prd:+.4f}")
    print("""
Как читать. Две разности сверху --- один и тот же шаг мёртвой зоны, приложенный к двум
разным углам. Если они близки, мёртвая зона делает своё дело независимо от того, что под
ней, и обе идеи складываются. Если на попарном углу шаг почти ничего не даёт, значит
попарная потеря уже выбросила ту же информацию, и в ансамбль идёт одна из них, а не обе.

Недобор в последней строке --- прямая величина перекрытия. Ноль означает сложение,
величина порядка самого эффекта означает, что это была одна мысль дважды.""")
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")


if __name__ == "__main__":
    main()
