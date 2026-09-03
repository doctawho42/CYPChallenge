"""The whole ensemble refitted against the band, not one dead-zone member added to four.

Why the earlier comparison was the wrong one. The dead-zone objective is worth -0.022 to -0.027
of pair and +0.027 to +0.037 of rank against a single per-enzyme boosting, on two seeds. Averaged
into the four-member ensemble it is worth -0.0036, and that number was read as "the ensemble
already captures it". It cannot mean that, because the comparison is malformed: one member fitted
against the band was being averaged with four fitted by squared error, so the ensemble's other
three quarters pull the average back toward the objective the dead zone exists to leave.

The honest question is what an ensemble does when **every** member is fitted against the band.

The reduction, and where it is exact. Item 122 established that the soft-threshold absolute loss
is L1 against a target reprojected onto the band, so alternating "reproject, refit" is a
majorise-minimise scheme. For the two boostings that is exact: scikit-learn's
`loss="absolute_error"` gives the L1 the reduction needs. The Gaussian process and the ridge
minimise squared error and have no L1 form, so reprojecting their target majorises a *different*
loss -- soft-threshold **squared** error rather than absolute. That is still a dead zone and still
declines to push inside the band; it is simply not the same estimator. The distinction is real
enough to be an arm rather than a footnote, so the two are reported separately:

    базовый ансамбль          all four by squared error. The control, and it must reproduce
                              0.6819 / 0.6009 or nothing below is comparable.
    мёртвая зона в бустингах  the two boostings refitted against the band, GP and ridge untouched.
                              The half where the reduction is exact.
    мёртвая зона везде        all four refitted against the band.

One pass, not three. Item 122 measured the majorise-minimise scheme at 0.6931, 0.6951, 0.6999 for
one, two and three passes: the training loss falls monotonically and the held-out rank peaks at
one. The scheme converges away from the optimum, so it is stopped where it is best rather than
where it stops moving.

Every reprojection is taken from **out-of-fold** predictions. Projecting a model onto its own
training predictions would put every training row inside its band, zero every gradient and turn
the refit into a no-op -- the trap item 122 documents.

Reads results/preds/oof_pool_all.json, oof_gp.json, oof_weak.json for the plain members and
recomputes the reprojected ones. Writes results/preds/oof_dzens.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import RidgeCV
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply
from gp import prepare as gp_prepare, gp_predict

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
DESC_MECH = 247
SRC = {"поферментно": ("oof_pool_all", "независимо"), "пул": ("oof_pool_all", "пул"),
       "GP": ("oof_gp", "GP"), "гребневая": ("oof_weak", "гребневая")}


def desc_scaled(X):
    A = np.nan_to_num(X[:, -DESC_MECH:].astype(np.float64), posinf=0.0, neginf=0.0)
    return np.clip((A - A.mean(0)) / (A.std(0) + 1e-9), -5.0, 5.0)


def refit(kind, X, target, fi, e):
    """Один проход мажорирования: та же модель, но по спроецированной мишени."""
    p = np.zeros(len(target))
    if kind == "GP":
        tf = gp_prepare(X)
        Xi = tf(X)
    elif kind == "гребневая":
        Xi = desc_scaled(X)
    elif kind == "пул":
        ind = np.zeros((len(X), len(CYPS)), np.float32)
        ind[:, e] = 1.0
        Xi = np.hstack([X, ind])
    else:
        Xi = X
    for f in range(5):
        trn, te = fi != f, fi == f
        if te.sum() == 0:
            continue
        if kind == "GP":
            p[te] = gp_predict(Xi[trn], target[trn], Xi[te])
        elif kind == "гребневая":
            p[te] = RidgeCV(alphas=np.logspace(-1, 4, 12)).fit(
                Xi[trn], target[trn]).predict(Xi[te])
        else:
            p[te] = HistGradientBoostingRegressor(
                **KW, loss="absolute_error").fit(Xi[trn], target[trn]).predict(Xi[te])
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--out", default=RES + "preds/oof_dzens.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)
    cache = {f: json.load(open(RES + f"preds/{f}.json"))["preds"]
             for f in {v[0] for v in SRC.values()}}

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        rec = {}
        for e, c in enumerate(CYPS):
            col = f"{c}_pIC50_direct_inhibition"
            m = tr[col].notna().to_numpy()
            y = tr.loc[m, col].to_numpy()
            lo = tr.loc[m, col + "_conf_low"].to_numpy()
            hi = tr.loc[m, col + "_conf_high"].to_numpy()
            Xi, fi = X[m], fold[m]

            plain, dz = {}, {}
            for nm, (fn, arm) in SRC.items():
                plain[nm] = np.asarray(cache[fn][f"{seed}|{arm}|{c}"], float)
                t0 = time.time()
                # Мишень --- проекция ЧЕСТНОГО предсказания вне фолда на полосу.
                dz[nm] = refit(nm, Xi, np.clip(plain[nm], lo, hi), fi, e)
                print(f"    {c} {nm:12s} перепод. за {time.time()-t0:.0f} с", flush=True)
            rec[c] = (y, lo, hi, fi, plain, dz)

        arms = {
            "базовый ансамбль": lambda p, d: [p[k] for k in SRC],
            "мёртвая зона в бустингах": lambda p, d: [d["поферментно"], d["пул"],
                                                      p["GP"], p["гребневая"]],
            "мёртвая зона везде": lambda p, d: [d[k] for k in SRC],
        }
        for nm, pick in arms.items():
            r = {"seed": seed, "рука": nm}
            for c in CYPS:
                y, lo, hi, fi, plain, dz = rec[c]
                q = np.mean(pick(plain, dz), axis=0)
                qq = fit_apply(q, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{nm}|{c}"] = q.tolist()
                r[f"{c} пара"] = round(float(strae(y, qq, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, q).statistic), 4)
            for tag in ("пара", "rho"):
                r[f"MACRO {tag}"] = round(float(np.mean([r[f"{c} {tag}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {nm:26s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS), flush=True)

        print("  члены поодиночке (ранг, обычный -> с мёртвой зоной):")
        for nm in SRC:
            v = [(spearmanr(rec[c][0], rec[c][4][nm]).statistic,
                  spearmanr(rec[c][0], rec[c][5][nm]).statistic) for c in CYPS]
            v = np.array(v).mean(0)
            print(f"    {nm:12s} {v[0]:.4f} -> {v[1]:.4f}  ({v[1]-v[0]:+.4f})", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                         + [f"{c} rho" for c in CYPS]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Первая рука обязана дать 0.6819 / 0.6009, иначе сравнивать не с чем.

Вторая против третьей отвечает на вопрос, ради которого это разделено: у бустингов редукция
точна (L1 против спроецированной мишени), у GP и гребневой --- нет, там получается мёртвая зона
с КВАДРАТОМ снаружи, другая оценка. Если третья лучше второй, приближение работает и на них;
если хуже --- мёртвую зону стоит ставить только туда, где редукция честная.

Строки по членам поодиночке говорят, кто именно выиграл от перепроецирования. Если выигрывают
все, а ансамбль нет, значит члены стали ПОХОЖЕ ошибаться --- и тогда мёртвая зона отнимает у
ансамбля разнообразие ровно в той мере, в какой добавляет каждому члену.""")


if __name__ == "__main__":
    main()
