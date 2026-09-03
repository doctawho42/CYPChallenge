"""The exact Bayes action under ST-RAE, instead of an affine pair fitted to it.

What this is for. Item 77 found that a fitted affine pair rewrites the predictions before
anything is scored, and that five measured gains did not survive it; item 80 turned that into
the criterion this repository now uses -- an intervention survives only if it adds rank. Item 93
then refuted a whole class of post-hoc corrections with one argument: a correction that is
monotone within a group can only move rank *between* groups, so the pair absorbs it.

This is the one correction that argument provably does not reach, and the reason is measured
rather than assumed.

The loss is fully known as a function of (prediction, truth). ST-RAE's numerator on a compound
is max(0, lo - p, p - hi), and item 114 found the band is a near-deterministic function of the
label: isotonic regression from y onto the width explains 0.93 to 0.97 of its variance, and the
label sits at an almost fixed relative position inside it (the fraction (y - lo)/w has a standard
deviation under 0.09). So lo and hi follow from y, and the only unknown left is y itself.

Under a posterior for y the optimal action is therefore

    p* = argmin_p  E_y [ max(0, lo(y) - p, p - hi(y)) ]

which is a one-dimensional problem per compound. Each term is convex and piecewise linear in p,
so the sum is convex and the minimum is exact rather than approximate.

Why the affine pair cannot absorb it. p* depends on the *spread* of the posterior, not only on
its location, because the loss is asymmetric -- the forgiveness threshold is about 0.13 above a
label of 4.6 and about 1.25 below 3.5 (item 115). A pair applies one lambda and one shift per
fold, which is monotone and two-parameter. That would be enough if the spread were recoverable
from the point prediction, and it is not: the best monotone function of the point explains

    CYP1A2 0.110   CYP2C9 0.114   CYP2D6 0.100   CYP3A4 0.080

of the spread's variance. About nine tenths of it is information the point does not carry, so
item 93's rule does not reach this and the pair has nothing to undo it with. That measurement is
the precondition for this file and it was run before the file was written.

How the posterior is built. Non-parametrically, and out of fold. The ensemble's member spread is
used only to *order* compounds by difficulty; the scale comes from the empirical out-of-fold
residuals of compounds in the same spread bin. This keeps the asymmetry of the residual
distribution, which is real here -- the label distribution is skewed and a Gaussian posterior
would misplace the action systematically at the weak end, which is exactly where the forgiveness
is widest. Residual distributions for a scored fold are collected from the other folds only.

Arms:

    сырое                the ensemble mean, untouched
    аффинная пара        what the submission does now
    байесово действие    this file
    байесово + пара      the diagnostic that matters most. If the pair still finds something to
                         fix on top of the Bayes action, the action is not doing its job; if it
                         finds nothing, the post-processing has been replaced rather than
                         supplemented.

Reads results/preds/oof_pool_all.json, oof_gp.json, oof_weak.json. Writes
results/preds/oof_bayes.json.
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
from sklearn.isotonic import IsotonicRegression
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
MEMBERS = [("oof_pool_all", "независимо"), ("oof_pool_all", "пул"),
           ("oof_gp", "GP"), ("oof_weak", "гребневая")]
NBIN = 5          # корзин по разбросу
NRES = 150        # остатков на корзину, подвыборка для скорости
GRID = 401        # точек сетки по p
HALF = 2.0        # полуширина сетки вокруг точечного предсказания


def band_model(y, lo, hi):
    """Ширина полосы как изотоника от метки, и доля метки внутри полосы.

    Возвращает пару (g, alpha): ширина w = g(y), а границы восстанавливаются как
    lo = y - alpha*w и hi = y + (1-alpha)*w. Обе величины --- из пункта 114, где измерено,
    что первая объясняет 0.93-0.97 дисперсии ширины, а вторая почти постоянна.
    """
    w = hi - lo
    g = IsotonicRegression(increasing=False, out_of_bounds="clip").fit(y, w)
    alpha = float(np.mean((y - lo) / np.maximum(w, 1e-9)))
    return g, alpha


def bayes_action(point, resid, g, alpha, half=HALF, grid=GRID):
    """Точный минимум ожидаемой потери по сетке, векторизованно.

    Для каждой гипотезы об истине y = point + r строятся её собственные границы, и потеря
    max(0, lo - p, p - hi) усредняется по r. Функция выпукла по p, так что сетка с шагом
    2*half/grid даёт минимум с точностью шага; при half=2 и grid=401 это 0.01 pIC50.
    """
    n = len(point)
    out = np.empty(n)
    off = np.linspace(-half, half, grid)
    step = 200
    for a in range(0, n, step):
        b = min(a + step, n)
        p0 = point[a:b]
        Y = p0[:, None] + resid[None, :]                 # (m, R) гипотезы об истине
        W = g.predict(Y.ravel()).reshape(Y.shape)
        LO, HI = Y - alpha * W, Y + (1.0 - alpha) * W
        P = p0[:, None] + off[None, :]                   # (m, G) кандидаты
        # (m, G, R): потеря каждого кандидата против каждой гипотезы
        d = np.maximum(0.0, np.maximum(LO[:, None, :] - P[:, :, None],
                                       P[:, :, None] - HI[:, None, :])).mean(axis=2)
        out[a:b] = P[np.arange(b - a), d.argmin(axis=1)]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/oof_bayes.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    cache = {f: json.load(open(RES + f"preds/{f}.json"))["preds"]
             for f in {m[0] for m in MEMBERS}}

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        t0 = time.time()
        r = {"seed": seed}
        for c in CYPS:
            col = f"{c}_pIC50_direct_inhibition"
            m = tr[col].notna().to_numpy()
            y = tr.loc[m, col].to_numpy()
            lo = tr.loc[m, col + "_conf_low"].to_numpy()
            hi = tr.loc[m, col + "_conf_high"].to_numpy()
            fi = fold[m]
            P = np.column_stack([np.asarray(cache[fn][f"{seed}|{arm}|{c}"], float)
                                 for fn, arm in MEMBERS])
            point, spread = P.mean(1), P.std(1)

            act = np.empty(len(y))
            for f in np.unique(fi):
                te, trn = fi == f, fi != f
                if te.sum() == 0 or trn.sum() < 50:
                    act[te] = point[te]
                    continue
                # Всё оценивается на обучающих фолдах: и модель полосы, и остатки.
                g, alpha = band_model(y[trn], lo[trn], hi[trn])
                res_tr = y[trn] - point[trn]
                # Разброс задаёт только ПОРЯДОК; масштаб берётся из остатков своей корзины.
                edges = np.quantile(spread[trn], np.linspace(0, 1, NBIN + 1)[1:-1])
                btr = np.digitize(spread[trn], edges)
                bte = np.digitize(spread[te], edges)
                rng = np.random.default_rng(seed * 100 + f)
                pool = {}
                for k in range(NBIN):
                    s = res_tr[btr == k]
                    if len(s) < 20:
                        s = res_tr
                    pool[k] = s if len(s) <= NRES else rng.choice(s, NRES, replace=False)
                idx = np.where(te)[0]
                for k in range(NBIN):
                    sel = idx[bte == k]
                    if len(sel) == 0:
                        continue
                    act[sel] = bayes_action(point[sel], pool[k], g, alpha)

            u = np.ones(len(y)) / len(y)
            arms = {"сырое": point,
                    "аффинная пара": fit_apply(point, lo, hi, fi, u),
                    "байесово действие": act,
                    "байесово + пара": fit_apply(act, lo, hi, fi, u)}
            for nm, p in arms.items():
                out[f"{seed}|{nm}|{c}"] = p.tolist()
                r[f"{c}|{nm}"] = float(strae(y, p, y_true_upper=hi, y_true_lower=lo))
                r[f"rho {c}|{nm}"] = float(spearmanr(y, p).statistic)
            # Монотонно ли действие относительно входа: если да, пара его съест.
            r[f"mono {c}"] = float(spearmanr(point, act).statistic)
            r[f"сдвиг {c}"] = float(np.mean(act - point))

        for nm in ("сырое", "аффинная пара", "байесово действие", "байесово + пара"):
            r[f"MACRO|{nm}"] = float(np.mean([r[f"{c}|{nm}"] for c in CYPS]))
            r[f"MACRO rho|{nm}"] = float(np.mean([r[f"rho {c}|{nm}"] for c in CYPS]))
        table.append(r)
        print(f"  сид {seed}: " + "  ".join(
            f"{nm} {r[f'MACRO|{nm}']:.4f}/{r[f'MACRO rho|{nm}']:.4f}"
            for nm in ("сырое", "аффинная пара", "байесово действие", "байесово + пара"))
            + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print("\nСреднее по сидам:")
    print(f"{'рука':22s} {'ST-RAE':>9s} {'ранг':>9s}")
    for nm in ("сырое", "аффинная пара", "байесово действие", "байесово + пара"):
        print(f"{nm:22s} {df[f'MACRO|{nm}'].mean():9.4f} {df[f'MACRO rho|{nm}'].mean():9.4f}")

    print("\nПоферментно, байесово действие против аффинной пары:")
    for c in CYPS:
        d = df[f"{c}|байесово действие"].mean() - df[f"{c}|аффинная пара"].mean()
        print(f"  {c}: {df[f'{c}|аффинная пара'].mean():.4f} -> "
              f"{df[f'{c}|байесово действие'].mean():.4f}  ({d:+.4f})   "
              f"монотонность к входу rho={df[f'mono {c}'].mean():+.4f}, "
              f"средний сдвиг {df[f'сдвиг {c}'].mean():+.3f}")

    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Решает третья строка против второй: аффинная пара --- это то, что делает сабмит
сейчас, и байесово действие обязано её обойти, иначе принципиальность конструкции ничего не
купила.

Четвёртая строка --- главная диагностика. Если пара НА ВЕРХУ байесова действия ещё что-то
находит, значит действие построено неверно: правильное действие уже минимизирует ту же
метрику, и паре нечего исправлять. Если четвёртая равна третьей --- постобработка не дополнена,
а ЗАМЕНЕНА выводом.

Монотонность к входу печатается рядом не для красоты. Если rho = 1.000, действие есть
монотонная функция точки, пункт 93 достаёт, и весь аргумент падает независимо от чисел.""")


if __name__ == "__main__":
    main()
