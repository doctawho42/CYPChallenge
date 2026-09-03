"""The split-normal likelihood the document specifies and the code never implemented.

The divergence, verified in the files rather than argued. Section 4 of the document defines a joint
generative model in which the curve enters through a **split normal** with two scales,

    p(pi) = split-N(pi; pihat, sigma^-, sigma^+),
    sigma^- = (pihat - l)/1.96,   sigma^+ = (h - pihat)/1.96      [s04.tex:51-53]

and `src/trunk.py` implements three masked squared errors (lines 330, 333, 346). `src/abloss.py`
records the asymmetry in its own docstring -- "on CYP1A2 the median distance down is 0.997 and up
0.702" -- and optimises a square anyway. The document and the code are two different projects at
this point, and this is the loss half of the difference.

**Why that is worth acting on rather than noting.** Item 189 collected six measurements showing
physics enters through the measurement model and returns zero through the feature matrix. There is
now a seventh, and it is sharper: the **two** places where the code stopped being an MSE and became
the document are the two interventions that worked --

    мёртвая зона     квадрат -> геометрия полосы            +0.033
    calshift         свободная голова -> уравнение прибора  +0.031

-- and every other part of section 4 has never reached the code at all.

**The honest objection, which decides the design.** The dead zone already uses the band, so a
split-normal may be the same information in another costume, and rule 166 puts changes of basis in
the column with four attempts and no survivors. Measuring split-normal against squared alone would
not separate those. So the arm is a **two-by-two**:

    квадрат                       эталон
    расщ. нормаль                 асимметричный вес по сторонам полосы
    квадрат + мёртвая зона        известное: +0.033
    расщ. нормаль + мёртвая зона  добавляет ли градуированный вес к жёсткому допуску

The implementation is exact rather than approximate. The negative log-likelihood of a split normal
has gradient `(p - y)/sigma^-^2` below the point and `(p - y)/sigma^+^2` above it, so it is a
squared loss with a **per-sample, per-side weight**, and gradient boosting fits it by scaling each
residual by `1/sigma^2` on whichever side the current prediction sits. No new machinery, no custom
objective -- the same reduction that made the dead zone implementable with a stock learner.

Pre-registered, and both outcomes are useful:

    расщ. нормаль помогает, но НЕ добавляет к мёртвой зоне  ->  та же информация в другой
                                                                одежде, правило 166, закрыто
    добавляет поверх мёртвой зоны                           ->  жёсткий допуск и градуированный
                                                                вес --- разные вещи, и документ
                                                                был прав дважды

Same folds, masks and metric. Writes results/preds/oof_split.json.
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
from sklearn.tree import DecisionTreeRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
NTREE, LR, DEPTH, MAXFEAT = 200, 0.06, 5, 0.3      # пины из src/ablpairloss.py
FLOOR = {"CYP1A2": 0.0061, "CYP2C9": 0.0071, "CYP2D6": 0.0049, "CYP3A4": 0.0033}
EPS = 1e-3


def boost(Xtr, ytr, Xte, seed, sm=None, sp=None, l1=False):
    """sm, sp --- сигмы вниз и вверх; заданы -> расщеплённая нормаль, иначе квадрат или L1."""
    rng = np.random.default_rng(seed)
    s = np.full(len(ytr), float(ytr.mean()))
    pred = np.full(len(Xte), s[0])
    for _ in range(NTREE):
        if l1:
            r = ytr - s
        elif sm is None:
            r = ytr - s
        else:
            # градиент -logL расщеплённой нормали: (y-p)/sigma^2 со стороной по знаку
            side = np.where(s < ytr, sm, sp)
            r = (ytr - s) / np.maximum(side, EPS) ** 2
            r = r / max(np.abs(r).max(), 1e-9) * np.abs(ytr - s).mean()
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MAXFEAT,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(Xtr, r)
        if l1:
            # LAD-бустинг Фридмана: дерево строится по ОБЫЧНЫМ остаткам, а значение листа
            # заменяется медианой остатков в нём. Это и есть шаг под абсолютную потерю.
            # Первая версия фитила sign(y-s) со средними в листьях --- это бустинг по знаку,
            # ДРУГОЙ оценщик, и он давал -0.027 там, где мёртвая зона стоит +0.032.
            # Тот же комментарий уже стоял в src/abldeadpair.py; функцию надо было
            # переиспользовать, а не писать заново.
            leaf = t.apply(Xtr)
            for lv in np.unique(leaf):
                t.tree_.value[lv, 0, 0] = np.median(r[leaf == lv])
        s = s + LR * t.predict(Xtr)
        pred = pred + LR * t.predict(Xte)
    return s, pred


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/oof_split.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)

    print(f"{'фермент':8s} {'sigma- медиана':>15s} {'sigma+ медиана':>15s} {'отношение':>10s}")
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna()
        y = tr.loc[m, col].to_numpy()
        sm = (y - tr.loc[m, col + "_conf_low"].to_numpy()) / 1.96
        sp = (tr.loc[m, col + "_conf_high"].to_numpy() - y) / 1.96
        print(f"{c:8s} {np.median(sm):15.4f} {np.median(sp):15.4f} "
              f"{np.median(sm)/max(np.median(sp),1e-9):10.3f}")
    print()

    ARMS = ["квадрат", "расщ. нормаль", "квадрат + МЗ", "расщ. нормаль + МЗ"]
    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        acc = {n: {"seed": seed, "рука": n} for n in ARMS}
        clock = {n: 0.0 for n in ARMS}
        for c in CYPS:
            col = f"{c}_pIC50_direct_inhibition"
            m = tr[col].notna().to_numpy()
            y = tr.loc[m, col].to_numpy()
            lo = tr.loc[m, col + "_conf_low"].to_numpy()
            hi = tr.loc[m, col + "_conf_high"].to_numpy()
            sm, sp = (y - lo) / 1.96, (hi - y) / 1.96
            Xi, fi = X[m], fold[m]

            base = {}
            for nm, kw in (("квадрат", {}), ("расщ. нормаль", {"sm": sm, "sp": sp})):
                t0 = time.time()
                p = np.zeros(len(y))
                for f in range(5):
                    trn, te = fi != f, fi == f
                    kk = {k: v[trn] for k, v in kw.items()}
                    p[te] = boost(Xi[trn], y[trn], Xi[te], seed * 10 + f, **kk)[1]
                base[nm] = p
                clock[nm] += time.time() - t0
            # проход мёртвой зоны поверх каждого угла: проекция честных OOF, затем L1
            for src, dst in (("квадрат", "квадрат + МЗ"),
                             ("расщ. нормаль", "расщ. нормаль + МЗ")):
                t0 = time.time()
                tgt = np.clip(base[src], lo, hi)
                p = np.zeros(len(y))
                for f in range(5):
                    trn, te = fi != f, fi == f
                    p[te] = boost(Xi[trn], tgt[trn], Xi[te], seed * 10 + f, l1=True)[1]
                base[dst] = p
                clock[dst] += time.time() - t0

            u = np.ones(len(y)) / len(y)
            for nm, p in base.items():
                q = fit_apply(p, lo, hi, fi, u)
                out[f"{seed}|{nm}|{c}"] = p.tolist()
                acc[nm][f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi,
                                                         y_true_lower=lo)), 4)
                acc[nm][f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)

        for nm in ARMS:
            r = acc[nm]
            for t in ("пара", "rho"):
                r[f"MACRO {t}"] = round(float(np.mean([r[f"{c} {t}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {nm:20s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({clock[nm]:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    g = df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                       + [f"{c} rho" for c in CYPS]].mean()
    print("\n" + g.round(4).to_string())
    if set(ARMS) <= set(g.index):
        a1 = g.loc["расщ. нормаль", "MACRO rho"] - g.loc["квадрат", "MACRO rho"]
        a2 = g.loc["расщ. нормаль + МЗ", "MACRO rho"] - g.loc["квадрат + МЗ", "MACRO rho"]
        print(f"\nэффект расщ. нормали БЕЗ мёртвой зоны : {a1:+.4f}")
        print(f"эффект расщ. нормали ПОВЕРХ мёртвой зоны: {a2:+.4f}")
        print(f"мёртвая зона на квадрате                : "
              f"{g.loc['квадрат + МЗ','MACRO rho']-g.loc['квадрат','MACRO rho']:+.4f}")
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Решают две строки внизу. Расщеплённая нормаль помогает без мёртвой зоны и не
добавляет поверх --- значит это та же информация о полосе в другой одежде, и по правилу 166
арка закрыта. Добавляет поверх --- значит жёсткий допуск и градуированный вес несут разное,
и документ был прав в обоих местах, а код терял оба.""")


if __name__ == "__main__":
    main()
