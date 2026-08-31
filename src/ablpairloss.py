"""Ranking loss counted only over pairs the metric can actually see.

The observation this rests on. After the affine pair the metric pays for order, but not for every
order: if two compounds' confidence bands overlap, no prediction can be penalised for putting them
either way round, so that pair is invisible to the score. Measured, the share of pairs with
disjoint bands is

    CYP1A2 73.6 %    CYP2C9 56.5 %    CYP2D6 73.4 %    CYP3A4 72.6 %

so on CYP2C9 nearly half the pairs do not exist as far as the metric is concerned. A plain ranking
loss spends capacity on them at the same rate as on the rest.

Why this is not the dead zone in another wrapper, which is the first objection and the right one.
The pointwise dead zone drops a compound by the absolute width of its own band. The pairwise
criterion keeps a weak compound whenever it is *distinguishable from something*, and two
wide-banded compounds can still be strictly ordered. The two sets are not the same, and item 115's
Jaccard overlap says where they differ most: 0.738 on CYP2C9, the lowest of the four. CYP2C9 is
also the enzyme with the sharpest number above, 56.5 per cent. That coincidence is not one -- both
follow from CYP2C9's bands being wide relative to its label spread -- and it makes CYP2C9 the
diagnostic enzyme here. If the two interventions were the same thing, they would agree there most,
not least.

Implementation, and why it is not HistGradientBoostingRegressor. A pairwise objective needs a
custom gradient, which the pinned scikit-learn's histogram booster does not accept. So the
gradient boosting is written out over plain regression trees -- fifty lines, fully controlled --
and that raises a confound immediately: any difference could be the loss or could be the learner.
Hence three arms, and the middle one is the whole point.

    HistGB, квадрат       the reference. Must reproduce 0.7150 / 0.5651 at seed 0.
    свой бустинг, квадрат the control. Same trees, same depth, same rate, squared loss. If this
                          does not land near the reference, the learner is the difference and the
                          third arm says nothing.
    свой бустинг, попарно the pairwise logistic loss over distinguishable pairs, uniform weights.

Uniform weights over pairs, deliberately. The objection to lambdarank was that NDCG's position
weights are a ranking convention we have no reason to import; a plain Kendall surrogate with equal
pair weights carries no such assumption, and the only weighting that survives is the one the
metric itself imposes -- a pair counts or it does not.

Isotonic on top, because the pairwise loss fixes order and says nothing about scale, and the
metric is scored in pIC50. Fitted out of fold like everything else.

Same folds, masks and metric as src/ablate.py. Writes results/preds/oof_pairloss.json.
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
from sklearn.isotonic import IsotonicRegression
from sklearn.tree import DecisionTreeRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
NTREE = 200
LR = 0.06
DEPTH = 5
MAXFEAT = 0.3


def pairwise_grad(s, y, lo, hi, chunk=512):
    """Градиент попарной логистической потери по РАЗЛИЧИМЫМ парам.

    Пара (i, j) считается только если полосы не пересекаются: иначе метрика не может её
    упорядочить ни при каком предсказании, и тратить на неё градиент значит учиться тому,
    за что не платят.
    """
    n = len(s)
    g = np.zeros(n)
    for a in range(0, n, chunk):
        b = min(a + chunk, n)
        # i --- строки куска, j --- все. Различимость и порядок берутся из полос.
        above = lo[a:b, None] > hi[None, :]        # истина i строго выше j
        if not above.any():
            continue
        d = s[a:b, None] - s[None, :]
        w = np.where(above, 1.0 / (1.0 + np.exp(np.clip(d, -30, 30))), 0.0)
        g[a:b] -= w.sum(axis=1)      # i --- верхний в паре, его счёт надо поднять
        g += w.sum(axis=0)           # j --- нижний, его опустить
    return g


def boost(Xtr, ytr, lotr, hitr, Xte, mode, seed):
    """Свой бустинг над обычными деревьями. mode: 'sq' или 'pair'."""
    rng = np.random.default_rng(seed)
    s = np.full(len(ytr), float(np.mean(ytr)) if mode == "sq" else 0.0)
    base = s[0]
    pred = np.full(len(Xte), base)
    for t in range(NTREE):
        if mode == "sq":
            resid = ytr - s
        else:
            resid = -pairwise_grad(s, ytr, lotr, hitr)
            nrm = np.abs(resid).max()
            if nrm > 0:
                resid = resid / nrm
        tree = DecisionTreeRegressor(max_depth=DEPTH, max_features=MAXFEAT,
                                     random_state=int(rng.integers(1 << 30)))
        tree.fit(Xtr, resid)
        s = s + LR * tree.predict(Xtr)
        pred = pred + LR * tree.predict(Xte)
    return s, pred


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--only", default="", help="через запятую, какие ферменты")
    ap.add_argument("--out", default=RES + "preds/oof_pairloss.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    cyps = CYPS if not a.only else [c for c in CYPS if c in a.only.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)

    print("доля пар с непересекающимися полосами:")
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna()
        lo = tr.loc[m, col + "_conf_low"].to_numpy()
        hi = tr.loc[m, col + "_conf_high"].to_numpy()
        g = np.random.default_rng(0)
        i, j = g.integers(0, len(lo), 100000), g.integers(0, len(lo), 100000)
        k = i != j
        print(f"  {c}: {100*((lo[i[k]] > hi[j[k]]) | (lo[j[k]] > hi[i[k]])).mean():.1f} %")
    print()

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for nm in ("HistGB, квадрат", "свой бустинг, квадрат", "свой бустинг, попарно"):
            t0 = time.time()
            r = {"seed": seed, "рука": nm}
            for c in cyps:
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
                    if nm.startswith("HistGB"):
                        p[te] = HistGradientBoostingRegressor(**KW).fit(
                            Xi[trn], y[trn]).predict(Xi[te])
                    else:
                        mode = "sq" if "квадрат" in nm else "pair"
                        s_tr, s_te = boost(Xi[trn], y[trn], lo[trn], hi[trn],
                                           Xi[te], mode, seed * 10 + f)
                        if mode == "pair":
                            # Попарная потеря задаёт порядок и молчит про шкалу; изотоника
                            # переводит счёт в pIC50, подогнанная на обучающих фолдах.
                            iso = IsotonicRegression(out_of_bounds="clip").fit(s_tr, y[trn])
                            p[te] = iso.predict(s_te)
                        else:
                            p[te] = s_te
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{nm}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for tag in ("пара", "rho"):
                r[f"MACRO {tag}"] = round(float(np.mean([r[f"{c} {tag}"] for c in cyps])), 4)
            table.append(r)
            print(f"  сид {seed} {nm:24s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in cyps)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                         + [f"{c} rho" for c in cyps]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Средняя рука решает, читаемы ли остальные. Свой бустинг с квадратом обязан лечь
рядом с HistGB; если он заметно хуже, разница третьей руки --- это разница обучателей, а не
потерь, и сравнивать её с первой нельзя.

CYP2C9 --- диагностический фермент. У него различимых пар меньше всех (56.5 %), и у него же
пересечение с мёртвой зоной по пункту 115 наименьшее (Жаккар 0.738). Если попарная потеря
делает то же, что мёртвая зона, они совпадут там СИЛЬНЕЕ всего; если разное --- слабее.""")


if __name__ == "__main__":
    main()
