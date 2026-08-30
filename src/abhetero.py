"""Shrinkage that depends on the compound, which is the one correction the affine pair cannot be.

The argument, and it is a proof rather than a hope. Item 77 found that metric-aligned monotone
post-processing absorbs almost every improvement, and item 80 explained why: a monotone map moves
scale and location and cannot change rank. The affine pair is fitted once per fold, so it is
GLOBAL within that fold - every compound gets the same lambda. A correction indexed by the
individual compound's own uncertainty is therefore not expressible by it, whatever the data say.
This is the only proposal on the table that escapes the ceiling by construction.

What is NOT this. Replacing the point estimate by the conditional median of the pooled band
bounds - the action section 10 derives, F_L + F_H = 1 - is already measured: it is the third arm
of src/abloss.py and it is worth 0.001. Quantile regression is heteroscedastic on its own, since
it gives every compound its own median. That is not the question. The question is whether the
amount of SHRINKAGE should differ per compound.

Two conditioners, because the honest one and the powerful one are different. The loss is defined
by the band, and band width varies by a factor of 5 to 19 between compounds, so band width is
what the correction should key on - but it is unknown at test time and has to be predicted. The
model's own predictive spread, from quantile regression at 0.1 and 0.9, is available directly.
Both are fitted out of fold and compared against the global pair on the same folds.

The relationship is not assumed monotone, and measurement says it should not be: mean absolute
error by band-width tercile runs 0.65 / 0.47 / 0.92 on CYP1A2, so the narrowest bands are harder
than the middle ones, not easier. Binning rather than a linear coefficient is deliberate for that
reason.

Falsification is the same as everywhere now: if the per-bin pair does not beat the global pair
out of fold, the correction was global after all and the ceiling holds.

Reads data/feats.npz and data/rows.csv. Writes results/preds/oof_hetero.json.
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
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import OFFGRID, LAMGRID

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)


def best_pair(p, lo, hi):
    """The (offset, lambda) minimising band loss on the rows given. Whole grid at once."""
    mu = p.mean()
    A = OFFGRID[:, None] * (1.0 - LAMGRID[None, :])
    B = np.broadcast_to(LAMGRID[None, :], A.shape)
    q = (mu * (1.0 - B) + A)[:, :, None] + B[:, :, None] * p[None, None, :]
    pen = (np.maximum(q - hi, 0.0) + np.maximum(lo - q, 0.0)).sum(axis=2)
    i, j = np.unravel_index(pen.argmin(), pen.shape)
    return mu + OFFGRID[i], LAMGRID[j]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--bins", type=int, default=3)
    ap.add_argument("--out", default=RES + "preds/oof_hetero.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])

    table = {}
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for c in CYPS:
            t0 = time.time()
            col = f"{c}_pIC50_direct_inhibition"
            m = tr[col].notna().to_numpy()
            y = tr.loc[m, col].to_numpy()
            lo = tr.loc[m, col + "_conf_low"].to_numpy()
            hi = tr.loc[m, col + "_conf_high"].to_numpy()
            Xi, fi = X[m], fold[m]
            p = np.zeros_like(y)
            wpred = np.zeros_like(y)   # предсказанная ширина полосы
            spread = np.zeros_like(y)  # собственный разброс модели, q90 - q10
            glob = np.zeros_like(y)
            het = {"ширина": np.zeros_like(y), "разброс": np.zeros_like(y)}

            for f in range(5):
                trn, te = fi != f, fi == f
                if te.sum() == 0:
                    continue
                Xt, Xe = Xi[trn], Xi[te]
                p[te] = HistGradientBoostingRegressor(**KW).fit(Xt, y[trn]).predict(Xe)
                # Ширина полосы как отдельная цель, в логарифме --- она положительна и скошена.
                lw = np.log(np.maximum(hi[trn] - lo[trn], 1e-6))
                wpred[te] = HistGradientBoostingRegressor(**KW).fit(Xt, lw).predict(Xe)
                qs = {}
                for tau in (0.1, 0.9):
                    qs[tau] = HistGradientBoostingRegressor(
                        **KW, loss="quantile", quantile=tau).fit(Xt, y[trn]).predict(Xe)
                spread[te] = qs[0.9] - qs[0.1]

                # На обучающих фолдах те же величины нужны, чтобы задать границы корзин
                # и подогнать пару внутри каждой. Считаются вложенной кросс-валидацией.
                ptr = np.zeros(int(trn.sum()))
                wtr = np.zeros(int(trn.sum()))
                str_ = np.zeros(int(trn.sum()))
                ftr = fi[trn]
                for g in np.unique(ftr):
                    it, ie = ftr != g, ftr == g
                    if ie.sum() == 0 or it.sum() < 20:
                        continue
                    Xa, Xb = Xt[it], Xt[ie]
                    ya = y[trn][it]
                    ptr[ie] = HistGradientBoostingRegressor(**KW).fit(Xa, ya).predict(Xb)
                    lwa = np.log(np.maximum(hi[trn][it] - lo[trn][it], 1e-6))
                    wtr[ie] = HistGradientBoostingRegressor(**KW).fit(Xa, lwa).predict(Xb)
                    q1 = HistGradientBoostingRegressor(**KW, loss="quantile", quantile=0.1).fit(Xa, ya).predict(Xb)
                    q9 = HistGradientBoostingRegressor(**KW, loss="quantile", quantile=0.9).fit(Xa, ya).predict(Xb)
                    str_[ie] = q9 - q1

                cg, lg = best_pair(ptr, lo[trn], hi[trn])
                glob[te] = cg + lg * (p[te] - cg)
                for nm, vtr, vte in (("ширина", wtr, wpred[te]), ("разброс", str_, spread[te])):
                    edges = np.quantile(vtr, np.linspace(0, 1, a.bins + 1)[1:-1])
                    btr, bte = np.digitize(vtr, edges), np.digitize(vte, edges)
                    for b in range(a.bins):
                        sa, sb = btr == b, bte == b
                        if sa.sum() < 30 or sb.sum() == 0:
                            het[nm][np.where(te)[0][sb]] = glob[te][sb]
                            continue
                        cb, lb = best_pair(ptr[sa], lo[trn][sa], hi[trn][sa])
                        het[nm][np.where(te)[0][sb]] = cb + lb * (p[te][sb] - cb)

            sc = lambda v: round(float(strae(y, v, y_true_upper=hi, y_true_lower=lo)), 4)
            table.setdefault(c, []).append({"seed": seed, "сырой": sc(p), "глобальная пара": sc(glob),
                                            "по ширине": sc(het["ширина"]),
                                            "по разбросу": sc(het["разброс"])})
            print(f"  сид {seed} {c}: сырой {sc(p):.4f}  глоб {sc(glob):.4f}  "
                  f"ширина {sc(het['ширина']):.4f}  разброс {sc(het['разброс']):.4f}  "
                  f"({time.time()-t0:.0f} с)", flush=True)

    print()
    rec = []
    for c in CYPS:
        d = pd.DataFrame(table[c]).mean(numeric_only=True)
        rec.append({"фермент": c, **{k: round(float(d[k]), 4) for k in
                                     ("сырой", "глобальная пара", "по ширине", "по разбросу")}})
    df = pd.DataFrame(rec)
    print(df.to_string(index=False))
    mac = df[["сырой", "глобальная пара", "по ширине", "по разбросу"]].mean()
    print(f"\nмакро: " + "  ".join(f"{k} {v:.4f}" for k, v in mac.items()))
    print(f"выигрыш по ширине {mac['по ширине']-mac['глобальная пара']:+.4f}, "
          f"по разбросу {mac['по разбросу']-mac['глобальная пара']:+.4f}")
    json.dump(table, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. «Глобальная пара» обязана лечь примерно на 0.715 макро --- это то, что делает
конвейер сейчас. Две правые колонки --- та же пара, но подобранная отдельно внутри корзин по
предсказанной ширине полосы и по собственному разбросу модели.

Если обе лягут на глобальную --- поправка была глобальной, и потолок пункта 77 держится даже
против вмешательства, которое по построению должно было его обойти. Это сильный отрицательный
результат, а не пустой.""")


if __name__ == "__main__":
    main()
