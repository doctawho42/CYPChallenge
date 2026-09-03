"""Three weak but differently-wrong models, added for diversity rather than for quality.

The lesson being generalised. Item 92 found a Gaussian process that is the worst of three models
on its own and improves the ensemble anyway, because its errors are decorrelated from the
boosting's. That is a regime, not an accident: once the ensemble is what matters, the useful
question about a candidate stops being "is it good" and becomes "is it wrong in a new way".

Three candidates, all cheap and all wrong differently:

  random forest   bagged trees rather than boosted. Same hypothesis class, opposite bias -
                  averaging many deep trees against summing many shallow ones;
  ridge           a linear function of descriptors. With 247 columns and 1285 to 2335 rows the
                  problem is well conditioned, and a linear model errs smoothly where trees err
                  in steps;
  kNN             distance-weighted neighbours in descriptor space, which item 90 measured as
                  the space where a neighbour actually carries activity information - one
                  neighbour's label alone reaches rank 0.261 there. Maximal decorrelation from
                  a tree, since it is the method a tree cannot imitate.

Hyperparameters are chosen on the training folds, never on the held-out one: ridge by leave-one-
out over an alpha grid, kNN by an inner split over k. Fixing them by hand would add a degree of
freedom chosen on the same data, which is the objection this repository keeps raising against
itself.

The expectation, recorded before the run: each is worse than the boosting alone and the ensemble
improves with them. If a member is both worse alone and neutral in the ensemble, it is simply
worse and diversity was not the missing ingredient.

Reads data/feats.npz. Writes results/preds/oof_weak.json.
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
from sklearn.ensemble import RandomForestRegressor
from sklearn.linear_model import RidgeCV
from sklearn.neighbors import KNeighborsRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
ALPHAS = np.logspace(-1, 4, 12)
KGRID = [5, 10, 20, 40]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/oof_weak.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    XA = np.hstack([z["FP"], z["DESC"], z["MECH"]])          # для леса --- всё
    XD = np.hstack([z["DESC"], z["MECH"]]).astype(np.float64)  # для линейных --- дескрипторы
    XD = np.nan_to_num(XD, posinf=0.0, neginf=0.0)
    XD = np.clip((XD - XD.mean(0)) / (XD.std(0) + 1e-9), -5.0, 5.0)
    print(f"лес на {XA.shape[1]} колонках, гребневая и kNN на {XD.shape[1]}\n")

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        r = {"seed": seed}
        for nm in ("лес", "гребневая", "kNN"):
            t0 = time.time()
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                fi = fold[m]
                Xi = XA[m] if nm == "лес" else XD[m]
                p = np.zeros_like(y)
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    if nm == "лес":
                        mdl = RandomForestRegressor(n_estimators=300, min_samples_leaf=2,
                                                    max_features=0.3, random_state=0, n_jobs=2)
                        p[te] = mdl.fit(Xi[trn], y[trn]).predict(Xi[te])
                    elif nm == "гребневая":
                        p[te] = RidgeCV(alphas=ALPHAS).fit(Xi[trn], y[trn]).predict(Xi[te])
                    else:
                        # k выбирается вложенным разбиением по тем же фолдам, не по отложенному
                        gi = fi[trn]
                        best, bk = None, KGRID[0]
                        for k in KGRID:
                            q = np.zeros(int(trn.sum()))
                            for g in np.unique(gi):
                                it, ie = gi != g, gi == g
                                if ie.sum() == 0 or it.sum() <= k:
                                    continue
                                q[ie] = KNeighborsRegressor(n_neighbors=k, weights="distance").fit(
                                    Xi[trn][it], y[trn][it]).predict(Xi[trn][ie])
                            s = float(np.mean((q - y[trn]) ** 2))
                            if best is None or s < best:
                                best, bk = s, k
                        p[te] = KNeighborsRegressor(n_neighbors=bk, weights="distance").fit(
                            Xi[trn], y[trn]).predict(Xi[te])
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{nm}|{c}"] = p.tolist()
                r[f"{nm}|{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{nm}|{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            mp = np.mean([r[f"{nm}|{c} пара"] for c in CYPS])
            mr = np.mean([r[f"{nm}|{c} rho"] for c in CYPS])
            print(f"  сид {seed} {nm:10s} пара {mp:.4f} rho {mr:.4f}  ({time.time()-t0:.0f} с)",
                  flush=True)
        table.append(r)

    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Число каждого поодиночке почти не важно --- все три заведомо слабее бустинга.
Важно, что даёт добавление в ансамбль, и считается это отдельно, из сохранённых предсказаний.

Порог провала: член, который и слабее в одиночку, и не двигает ансамбль, --- просто слабее, и
разнообразие тут не было недостающим.""")


if __name__ == "__main__":
    main()
