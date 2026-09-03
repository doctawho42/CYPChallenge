"""The ensemble is averaged, and the metric is not minimised by an average.

The observation. `src/submit.py` line 230 and `src/abldzens.py` line 147 both combine members by
`np.mean`, and in 182 items nothing else has been tried. But the loss is

    d(p) = max(0, lo - p, p - hi) = расстояние от точки до отрезка,

and the minimiser of an expected distance-to-a-set is not an arithmetic mean. Differentiating,

    d/dp E[d(p)] = -P(lo > p) + P(hi < p),

so the optimum is the point where **the probability mass of bands lying entirely above equals the
mass lying entirely below** -- a generalised median of the band structure, not a centre of mass.

Why this is not item 126 again. That measured the exact Bayes action as a **post-hoc per-compound
correction** under a posterior, competing against the affine pair, and lost. This changes the
**combination operator**, which sits upstream of the pair: the pair is a monotone map of whatever
vector it receives, so a different vector is not something it can reproduce. Item 77's ceiling is a
ceiling on monotone post-processing of a fixed prediction, and this does not produce a fixed
prediction.

**Two preconditions were checked before building, and one nearly killed it.** On a symmetric
predictive distribution with symmetric bands the balancing point coincides with the mean --
verified numerically at four decimals, gain 0.0 per cent. So the idea has content only if something
is asymmetric. Measured on the training bands:

    фермент   вверх   вниз   средняя асимметрия   доля |асимметрии| > 0.05
    CYP1A2    0.277  0.306              -0.0209                     53.8 %
    CYP2C9    0.339  0.369              +0.0009                     64.3 %
    CYP2D6    0.227  0.265              -0.0025                     64.0 %
    CYP3A4    0.404  0.474              +0.0057                     72.4 %

Asymmetric per compound in a majority of cases, symmetric on average. That is the awkward shape: it
is not a bias to correct, it is a per-compound adjustment with no systematic direction -- and the
part of it that is systematic in the label is exactly the part a symmetric predictor already gets
right.

So the question this settles is narrow and worth settling: **do four real members, with bands
predicted from their own predictions, place their balancing point anywhere other than their mean,
and does it score better.**

The band at prediction time comes from item 114's relationship -- the band is a near-deterministic
function of the label -- fitted on training folds as `h(y)` and evaluated at each member's own
prediction, so no held-out information is used.

    среднее          то, что делается сейчас
    медиана          самый дешёвый соперник, тоже не среднее
    баланс полос     точка, где число полос целиком выше равно числу целиком ниже

Pre-registered: if the balancing point is numerically indistinguishable from the mean, the operator
is right and the question closes -- which is a useful thing to know about a line of code nobody has
questioned. If it differs and scores better, the ensemble has been combined by the wrong rule since
the beginning.

Reads oof_pool_all.json, oof_gp.json, oof_weak.json. Writes nothing.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
SRC = [("oof_pool_all", "независимо"), ("oof_pool_all", "пул"),
       ("oof_gp", "GP"), ("oof_weak", "гребневая")]


def balance(P, HL, HU):
    """Точка, где число полос целиком выше равно числу целиком ниже.

    P --- (M, n) предсказания членов, HL/HU --- их полуширины вниз и вверх. Решается
    точно: кандидатами являются только концы полос, между ними субградиент постоянен.
    """
    lo, hi = P - HL, P + HU
    cand = np.concatenate([lo, hi], axis=0)              # (2M, n)
    best = np.empty(P.shape[1])
    for j in range(P.shape[1]):
        c = np.unique(cand[:, j])
        d = np.array([np.mean(np.maximum(0.0, np.maximum(lo[:, j] - p, p - hi[:, j])))
                      for p in c])
        best[j] = c[int(d.argmin())]
    return best


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    blobs = {}
    for fn, _ in SRC:
        if fn not in blobs:
            blobs[fn] = json.load(open(RES + f"preds/{fn}.json"))["preds"]

    seeds = sorted({int(k.split("|")[0]) for k in blobs["oof_pool_all"]
                    if k.split("|")[0].isdigit()})
    print(f"сиды: {seeds}, членов: {len(SRC)}\n")

    rec = []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for c in CYPS:
            col = f"{c}_pIC50_direct_inhibition"
            m = tr[col].notna().to_numpy()
            y = tr.loc[m, col].to_numpy()
            lo = tr.loc[m, col + "_conf_low"].to_numpy()
            hi = tr.loc[m, col + "_conf_high"].to_numpy()
            fi = fold[m]
            P = []
            for fn, arm in SRC:
                k = f"{seed}|{arm}|{c}"
                if k not in blobs[fn] or len(blobs[fn][k]) != int(m.sum()):
                    P = None
                    break
                P.append(np.asarray(blobs[fn][k], float))
            if P is None:
                continue
            P = np.vstack(P)

            # Полуширины из отношения "полоса как функция метки", подогнанного НА ОБУЧАЮЩИХ
            # фолдах и вычисленного в предсказании каждого члена.
            HL = np.zeros_like(P); HU = np.zeros_like(P)
            for f in range(5):
                trn, te = fi != f, fi == f
                if te.sum() == 0:
                    continue
                il = IsotonicRegression(out_of_bounds="clip").fit(y[trn], (y - lo)[trn])
                iu = IsotonicRegression(out_of_bounds="clip").fit(y[trn], (hi - y)[trn])
                for mm in range(P.shape[0]):
                    HL[mm, te] = il.predict(P[mm, te])
                    HU[mm, te] = iu.predict(P[mm, te])

            rules = {"среднее": P.mean(0), "медиана": np.median(P, 0),
                     "баланс полос": balance(P, HL, HU)}
            for nm, p in rules.items():
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                rec.append(dict(seed=seed, фермент=c, правило=nm,
                                пара=float(strae(y, q, y_true_upper=hi, y_true_lower=lo)),
                                rho=float(spearmanr(y, p).statistic),
                                отклонение=float(np.mean(np.abs(p - P.mean(0))))))
        print(f"  сид {seed} готов", flush=True)

    df = pd.DataFrame(rec)
    g = df.groupby("правило")[["пара", "rho", "отклонение"]].mean()
    print("\n" + g.round(4).to_string())
    print("\nпо ферментам, ранг:")
    print(df.pivot_table(index="фермент", columns="правило", values="rho").round(4).to_string())
    print("""
Как читать. Столбец «отклонение» --- среднее |правило минус среднее| в единицах pIC50. Если
он около нуля, балансирующая точка И ЕСТЬ среднее, вопрос закрыт, и строка np.mean в submit.py
правильна не по случайности. Если он заметен, а ранг и пара не двигаются --- операторы
расходятся, но метрике всё равно, что тоже ответ.

Предрегистрация: выигрыш возможен только за счёт ПОЭКЗЕМПЛЯРНОЙ асимметрии полос, потому что
на симметричном случае баланс тождествен среднему --- проверено численно до четвёртого знака.""")


if __name__ == "__main__":
    main()
