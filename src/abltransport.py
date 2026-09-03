"""Optimal transport instead of importance weighting, as the second answer to one impossibility.

Why this is not a loophole around item 123. That item bounded *self-normalised importance
weights*: the chi-square divergence of 2.838 caps their effective sample size at 26 per cent, and
no choice of estimator moves that bound. Optimal transport does not weight. It **moves mass**:
training points are mapped barycentrically onto the test support and the model is fitted on the
mapped points. That is a different estimator with a different failure mode, so the bound does not
apply -- not because we found a gap in it, but because it was never a statement about this class
(Courty et al., domain adaptation with OT).

Together with `src/ablrobust.py` this gives two mathematically distinct responses to the same
proven impossibility: minimise the worst case inside the measured ball, or transport onto the
target support. They can disagree, and if they do the disagreement is informative.

The honest cost, recorded before the run. Transport moves the features and assumes the labels
travel with them. That assumption is exactly what an activity cliff violates, and this project has
measured its cliffs: item 97's non-monotonicity across similarity strata and item 133's finding
that near neighbours disagree more than any smooth model predicts. **So this arm can fail for a
reason we already know**, and a failure should be read as confirming that rather than as a verdict
on transport.

What it is computed on. The transport runs in the standardised DESC+MECH block, 247 dimensions,
not on the full 2295. Morgan counts are sparse and high-dimensional; a Euclidean ground cost there
is dominated by bit-count differences and the resulting plan is close to meaningless. The
fingerprint still enters the model -- only the ground cost is restricted.

Protocol. One plan, fitted from all training rows onto the test set, and both the training rows
and the held-out rows pass through it, so the model is trained and evaluated in one space. Fitting
it per fold would be the reflex here and it is the wrong reflex twice over: the transport reads
**features only**, and the test's features are published, so no label of a held-out row enters the
plan at any step; and barycentric mapping is defined only for the points the plan couples, so
held-out rows would have to be carried in by interpolation, adding an approximation to avoid a
leak that does not exist.

An `alpha` sweep rather than a single mapped set: alpha is the fraction of the way from a
compound's own coordinates to its barycentre on the test support, and alpha = 0 leaves the
features untouched and must reproduce the baseline exactly.

Evaluation is the same pair as the robust arm, and the pre-registration is the same:

  the mapped model **must lose** on the uniform out-of-fold rank -- it pays for the move;
  it **must win** on the test-regime-weighted rank of item 113;
  losing on both means the transport degraded the features and nothing else.

Writes results/preds/oof_ot.json.
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
from scipy.stats import spearmanr, rankdata
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from cypsplit import butina_folds
from shrinkchoice import fit_apply
import feats as F

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
DESC_MECH = 247
BINS = np.array([0.0, .25, .30, .35, .40, .45, .50, .55, .60, .70, 1.01])
REG = 1e-1          # энтропийная регуляризация Синкхорна


def scaled(A, mu=None, sd=None):
    A = np.nan_to_num(np.asarray(A, np.float64), posinf=0.0, neginf=0.0)
    if mu is None:
        mu, sd = A.mean(0), A.std(0) + 1e-9
    return np.clip((A - mu) / sd, -5.0, 5.0), mu, sd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--alphas", default="0,0.25,0.5,1.0",
                    help="доля пути к тестовой опоре: 0 --- без переноса")
    ap.add_argument("--out", default=RES + "preds/oof_ot.json")
    a = ap.parse_args()
    import ot

    seeds = [int(x) for x in a.seeds.split(",")]
    alphas = [float(x) for x in a.alphas.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)

    desc_names = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
    mech_names = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    FPt, dsc, M, ok = F.build(list(te.SMILES), desc_names, mech_names)
    Xte = np.hstack([FPt, dsc.to_numpy(np.float32), M.to_numpy(np.float32)])
    print(f"тест построен: {Xte.shape}, обучение {X.shape}")

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        s_oof = np.zeros(len(rows))
        for f in range(5):
            t_, tn = np.where(fold == f)[0], np.where(fold != f)[0]
            ref = [fps[j] for j in tn]
            for i in t_:
                s_oof[i] = max(DataStructs.BulkTanimotoSimilarity(fps[i], ref))
        rng = np.random.default_rng(0)
        n_sub = int(round(0.8 * len(fps)))
        fte = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in te.SMILES
               if Chem.MolFromSmiles(s) is not None]
        s_te = np.zeros(len(fte))
        for _ in range(5):
            ref = [fps[j] for j in rng.choice(len(fps), n_sub, replace=False)]
            for t_, f_ in enumerate(fte):
                s_te[t_] += max(DataStructs.BulkTanimotoSimilarity(f_, ref)) / 5
        h_o = np.histogram(s_oof, BINS)[0] / len(s_oof)
        h_t = np.histogram(s_te, BINS)[0] / len(s_te)
        ratio = np.where(h_o > 0, h_t / np.maximum(h_o, 1e-9), 1.0).clip(0, 8.0)
        w_test = ratio[np.clip(np.digitize(s_oof, BINS) - 1, 0, len(ratio) - 1)]

        Bt, mu, sd = scaled(Xte[:, -DESC_MECH:])
        # План строится один раз на всех обучающих строках, и это не утечка: транспорт
        # использует ТОЛЬКО признаки, а признаки теста опубликованы. Метка отложенного фолда
        # в план не входит ни на каком шаге, поэтому разбивать его по фолдам незачем ---
        # и вредно, потому что барицентрическое отображение определено лишь для точек,
        # участвовавших в плане, и отложенные пришлось бы доносить интерполяцией.
        t0 = time.time()
        Anat, _, _ = scaled(X[:, -DESC_MECH:], mu, sd)
        Ck = ot.dist(Anat, Bt)
        Ck /= max(Ck.max(), 1e-12)
        G = ot.sinkhorn(np.ones(len(Anat)) / len(Anat), np.ones(len(Bt)) / len(Bt),
                        Ck, REG, numItermax=300)
        Amap = (G / np.maximum(G.sum(1, keepdims=True), 1e-12)) @ Bt
        print(f"  сид {seed}: план {G.shape} посчитан за {time.time()-t0:.0f} с; "
              f"средний сдвиг точки {np.linalg.norm(Amap - Anat, axis=1).mean():.3f}",
              flush=True)

        for al in alphas:
            t0 = time.time()
            blend = (1 - al) * Anat + al * Amap
            XX = np.hstack([X[:, :-DESC_MECH], blend.astype(np.float32)])
            r = {"seed": seed, "alpha": al}
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                Xi, fi = XX[m], fold[m]
                p = np.zeros_like(y)
                for f in range(5):
                    trn, tei = fi != f, fi == f
                    if tei.sum() == 0:
                        continue
                    p[tei] = HistGradientBoostingRegressor(**KW).fit(
                        Xi[trn], y[trn]).predict(Xi[tei])
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{al:.2f}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
                ww = w_test[m]
                ra, rb = rankdata(y).astype(float), rankdata(p).astype(float)
                ra -= (ww * ra).sum() / ww.sum()
                rb -= (ww * rb).sum() / ww.sum()
                r[f"{c} rho_test"] = round(float((ww * ra * rb).sum() / np.sqrt(
                    (ww * ra * ra).sum() * (ww * rb * rb).sum())), 4)
            for tag in ("пара", "rho", "rho_test"):
                r[f"MACRO {tag}"] = round(float(np.mean([r[f"{c} {tag}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  alpha {al:4.2f}  пара {r['MACRO пара']:.4f}  ранг {r['MACRO rho']:.4f}  "
                  f"ранг-под-тест {r['MACRO rho_test']:.4f}  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("alpha")[["MACRO пара", "MACRO rho", "MACRO rho_test"]].mean().to_string())
    print("\nтолько 2C9 и 3A4 (ESS 36.7 и 39.3 %):")
    print(df.groupby("alpha")[["CYP2C9 rho_test", "CYP3A4 rho_test"]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. alpha = 0 --- признаки не тронуты, и строка обязана воспроизвести базу; иначе
ошибка в проводке, и остальное не значит ничего.

Предрегистрация та же, что у робастной руки, и сравнивать их надо между собой, а не только с
базой: это два разных ответа на одну и ту же доказанную невозможность. Если оба выигрывают по
ранг-под-тест --- режим действительно достижим двумя способами. Если оба проигрывают ---
проигрывает не метод, а посылка, что тестовый режим вообще достижим из обучающего.

И записанная заранее причина возможного провала: перенос двигает признаки, полагая, что метки
едут с ними. Обрывы активности это ломают, и они у нас измерены (пункты 97 и 133). Провал
именно с этой стороны надо читать как подтверждение обрывов, а не как приговор транспорту.""")


if __name__ == "__main__":
    main()
