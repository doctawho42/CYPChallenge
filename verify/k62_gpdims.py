"""Is the GP's isotropic kernel diluted by 247 descriptors, and does CYP3A4 pay for it?

Item 218 dropped the ensemble on CYP3A4 for the Gaussian process alone. That closes one question
and opens another that did not exist before: **on CYP3A4 there is no longer an ensemble to absorb
an improvement to the GP.** Items 176 and 182 are two measurements of single-member gains failing
to reach the ensemble -- +0.031 arriving as -0.0008, +0.029 as +0.0014, at error correlations of
0.94 to 0.97. That mechanism is gone on this enzyme: the member IS the model, so a gain transfers
one for one.

So the GP is worth looking at, and the obvious place is its kernel. `gp.py` uses an isotropic RBF,

    k(a, b) = exp(-||a - b||^2 / 2 (ls0 * lm)^2)

with one lengthscale for all 247 standardised descriptors, `ls0` from the median-distance
heuristic and `lm` fitted by marginal likelihood on training folds. Every dimension enters the
distance with equal weight, including the ones carrying nothing for this enzyme. In 247 dimensions
that is a real cost in principle -- the question is whether it is one here.

**This file does not fix it. It measures whether there is anything to fix**, which is cheaper by a
day and decides whether the fix is worth designing. Automatic relevance determination would give
each dimension its own lengthscale, 247 parameters fitted jointly, and item 178 is the standing
warning about what more flexibility does here. Before paying that, the diluted-kernel hypothesis
should make a prediction that can fail.

**The prediction.** If the kernel is diluted, then dropping the least informative dimensions should
*improve* it: the same GP on the top k descriptors by absolute Spearman correlation with the label
-- ranked on the training folds only, so the selection never sees the rows it is applied to --
should beat the GP on all 247 for some k below 247. If no k beats the full set, the isotropic
kernel is not paying a dilution cost worth chasing and ARD is closed for an hour's work rather than
a week's.

Selection by univariate correlation is a deliberately weak instrument, and that is the point: it is
247 one-dimensional fits rather than one 247-dimensional one, so it cannot overfit the way a joint
optimisation can, and anything it finds is a lower bound on what a real weighting would.

All four enzymes are run, not only CYP3A4 -- the other three still average five members, so a gain
there is diluted, but a **null** there is still informative about the kernel.

Reads data/feats.npz and data/rows.csv. Four seeds. Minutes per arm.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply
from gp import prepare as gp_prepare, gp_predict

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
DESC_MECH = 247


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--ks", default="247,128,64,32,16")
    ap.add_argument("--out", default=RES + "preds/oof_gpdims.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    ks = [int(x) for x in a.ks.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])[:, -DESC_MECH:]
    print(f"блок DESC+MECH: {X.shape}\n", flush=True)

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for k in ks:
            t0 = time.time()
            r = {"seed": seed, "k": k}
            for e, c in enumerate(CYPS):
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                Xi, fi = X[m], fold[m]
                p = np.zeros(len(y))
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    if k >= X.shape[1]:
                        keep = np.arange(X.shape[1])
                    else:
                        # Отбор ТОЛЬКО по обучающим фолдам: столбцы ранжируются по модулю
                        # спирменовой связи с меткой, и отложенный фолд в ранжировании
                        # не участвует ни одной строкой.
                        s_ = np.abs([spearmanr(Xi[trn, j], y[trn]).statistic
                                     for j in range(X.shape[1])])
                        s_ = np.nan_to_num(s_)
                        keep = np.argsort(-s_)[:k]
                    tf = gp_prepare(Xi[:, keep])
                    Z = tf(Xi[:, keep])
                    p[te] = gp_predict(Z[trn], y[trn], Z[te])
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{k}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for t in ("пара", "rho"):
                r[f"MACRO {t}"] = round(float(np.mean([r[f"{c} {t}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} k={k:4d}  макро rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.4f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    g = df.groupby("k", sort=False)[["MACRO rho"] + [f"{c} rho" for c in CYPS]].mean()
    print()
    print(g.round(4).to_string())
    full = g.loc[max(ks)]
    print()
    for k in ks:
        if k == max(ks):
            continue
        d = np.array([df[(df.seed == s) & (df.k == k)]["CYP3A4 rho"].iloc[0]
                      - df[(df.seed == s) & (df.k == max(ks))]["CYP3A4 rho"].iloc[0]
                      for s in seeds])
        print(f"k={k:4d} на CYP3A4: {d.mean():+.4f} знак {int((d > 0).sum())}/4  "
              + " ".join(f"{v:+.4f}" for v in d))
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Пол CYP3A4 --- 0.0033, и на нём ансамбля больше нет, так что выигрыш проходит
целиком, а не одной пятой.

Ни одно k не обыгрывает полный набор --- изотропное ядро не платит за разбавление ничего, что
стоило бы гнать, и ARD закрывается часом вместо недели. Какое-то обыгрывает --- разбавление
реально, и тогда правильная постройка не отбор столбцов (он груб), а поразмерные длины
корреляции; но её стоит проектировать, только когда предсказание уже подтвердилось.""")


if __name__ == "__main__":
    main()
