"""Train on the metric the competition scores, instead of on the squared error.

The gap this closes. Section 10 derives the Bayes-optimal action under ST-RAE and gets a
condition, F_L(yhat) + F_H(yhat) = 1: the optimum is the median of the pooled sample of both
band bounds, not the conditional mean of the point estimate. That derivation has been in the
document since the first draft, and every model in this repository is nonetheless fitted by
least squares against the point label. We knew the right action and minimised something else.

The empirical analogue of the theorem is three lines: duplicate every training row, once with
`lo` as its target and once with `hi`, and fit the median. No new dependency - the pinned
scikit-learn has `loss="quantile"`.

Three arms, because a win has two possible causes and only one of them is the theorem:

  L2 on the point label     - what every other file here does. The control;
  L1 on the point label     - the median of the SAME target. Isolates the loss change,
                              which is mostly robustness to outliers;
  L1 on the pooled bounds   - the loss change AND the target change together. This is the
                              theorem's action.

Without the middle arm a win could be nothing but robustness. Item 42 found a single compound
whose prediction tripled an enzyme's ST-RAE, so this repository has direct evidence that
outliers matter here, and attributing their removal to the band geometry would be wrong.

The bands are asymmetric - on CYP1A2 the median distance down is 0.997 and up 0.702 - so the
pooled median is not the point estimate and the two arms should differ. That is an argument
rather than a measurement, which is why this is run rather than asserted.

Same folds, same masks, same metric, same learner settings as src/ablate.py. Reads
data/feats.npz and data/rows.csv, writes results/preds/oof_loss.json.
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
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
ARMS = {
    "L2 по метке":    dict(loss="squared_error"),
    "L1 по метке":    dict(loss="absolute_error"),
    "L1 по границам": dict(loss="quantile", quantile=0.5),
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--out", default=RES + "preds/oof_loss.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])

    print(f"{'фермент':8s} {'медиана вниз':>13s} {'медиана вверх':>14s} {'асимметрия':>11s}")
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna()
        d = (tr.loc[m, col] - tr.loc[m, col + "_conf_low"]).median()
        u = (tr.loc[m, col + "_conf_high"] - tr.loc[m, col]).median()
        print(f"{c:8s} {d:13.3f} {u:14.3f} {d - u:+11.3f}")
    print()

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for name, kw in ARMS.items():
            t0 = time.time()
            r = {"seed": seed, "рука": name}
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                Xi, fi = X[m], fold[m]
                p = np.zeros_like(y)
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    if name == "L1 по границам":
                        # Каждая строка входит дважды, с нижней и верхней границей.
                        # Медиана объединённой выборки --- это и есть оптимум из §10.
                        Xt = np.vstack([Xi[trn], Xi[trn]])
                        yt = np.concatenate([lo[trn], hi[trn]])
                    else:
                        Xt, yt = Xi[trn], y[trn]
                    p[te] = HistGradientBoostingRegressor(**KW, **kw).fit(Xt, yt).predict(Xi[te])
                out[f"{seed}|{name}|{c}"] = p.tolist()
                r[c] = round(float(strae(y, p, y_true_upper=hi, y_true_lower=lo)), 4)
            r["MACRO"] = round(float(np.mean([r[c] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {name:16s} макро {r['MACRO']:.4f}  ({time.time()-t0:.0f} с)",
                  flush=True)

    print()
    print(pd.DataFrame(table)[["seed", "рука", *CYPS, "MACRO"]].to_string(index=False))
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Первая рука обязана дать 0.7673, иначе сравнивать не с чем.

Если выигрывает уже ВТОРАЯ --- дело не в теореме, а в устойчивости к выбросам, и полосы
тут ни при чём. Если вторая стоит на месте, а третья двигается --- работает именно замена
цели на границы, то есть оптимальное действие из §10 стоит того, чтобы его минимизировать,
а не только выводить.

Ноль в третьей руке тоже результат: он означает, что на этих данных расхождение между
квадратичной ошибкой и интервальной L1 не стоит ничего, и раздел §10 остаётся верным
утверждением без последствий.""")


if __name__ == "__main__":
    main()
