"""Distributionally robust training over similarity strata, with the radius measured rather
than chosen.

The idea, and why it belongs to this project specifically. Item 123 proved that a validation set
in the test's regime cannot be built: the chi-square divergence between the test's
nearest-neighbour similarity distribution and our out-of-fold one is 2.838, which caps any
importance-weighted estimate at 26 per cent effective sample size, and item 129's composition
result then showed it cannot be built by construction either. That was recorded as the limit of
what is measurable.

It is also a radius. The decision-theoretic answer to "the true distribution is unknown but lies
within a ball of size rho" is not to estimate it but to minimise the worst case inside the ball:

    min_f  sup_{Q : chi2(Q || P) <= rho}  E_Q[L(f)]

What is unusual here is not the method but the constant. In the literature rho is a tuning
parameter chosen by cross-validation or by taste, which is the standard criticism of the whole
approach -- a robust answer to an arbitrary question. Ours is **measured**, from the published
test structures against our own folds, so the uncertainty set is built from data and contains the
test distribution by construction.

The objection that has to be answered before the arm is worth running, and is answered by the
design rather than by argument. An unrestricted chi-square ball of radius 2.838 permits any
reweighting at all, including putting the whole mass on the hardest compounds; minimising the
worst case over that is not robustness, it is underfitting. But the 2.838 was measured over a
*partition* -- bins of similarity -- and the honest uncertainty set is the one that partition
supports: a minimax over the stratum weights only. So the ball here is over five similarity
strata and nothing finer, and the radius is recomputed for exactly the partition used, never
borrowed from a different binning.

The inner problem has a closed form, which is what makes this a reweighting rather than a new
learner. Maximising sum_k q_k L_k subject to sum q_k = 1 and sum (q_k - p_k)^2 / p_k <= rho gives

    q_k = p_k * (1 + sqrt(rho / Var_p(L)) * (L_k - mean_p(L)))

so the per-sample weight is 1 + sqrt(rho / Var_p(L)) * (L_{k(i)} - mean_p(L)), clipped at zero.
Strata that are already easy get down-weighted, hard ones up-weighted, and the total displacement
is exactly the measured radius.

The loss driving the weights is ST-RAE's own numerator per compound -- max(0, lo - p, p - hi) --
not squared error, because the strata have to be judged on the quantity the competition scores.

Pre-registered reading, and it is two-sided on purpose:

  DRO **must lose** on the uniform out-of-fold metric. It buys robustness with average-case
  performance, and an arm that wins on both is not doing what it claims.
  DRO **must win** on the test-regime-weighted rank of item 113. That is the only estimate we
  have of the regime it is built for.
  Losing on both is underfitting and closes the arch.

A sweep over rho rather than a single value, since rho = 0 is exactly ERM and reproduces the
baseline -- a control that costs nothing and catches a wiring error immediately. Judge on CYP2C9
and CYP3A4, whose effective sample sizes under the test-regime weights are 36.7 and 39.3 per cent;
on CYP1A2 and CYP2D6 item 123 measured the same criterion at 24.8 and 20.7 per cent, four times
noisier, so those two need four seeds before they mean anything.

Writes results/preds/oof_dro.json.
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

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
NSTRAT = 5
BINS = np.array([0.0, .25, .30, .35, .40, .45, .50, .55, .60, .70, 1.01])


def nn_similarity(fps, fold):
    """Сходство с ближайшим обучающим ЧЕРЕЗ фолд --- та же величина, что в пунктах 97 и 123."""
    s = np.zeros(len(fps))
    for f in np.unique(fold):
        te, trn = np.where(fold == f)[0], np.where(fold != f)[0]
        ref = [fps[j] for j in trn]
        for i in te:
            s[i] = max(DataStructs.BulkTanimotoSimilarity(fps[i], ref))
    return s


def strata_and_radius(s_oof, s_te, nstrat=NSTRAT):
    """Разбиение по квантилям OOF-сходства и chi2 ДЛЯ ЭТОГО ЖЕ разбиения.

    Радиус пересчитывается под используемую сетку, а не берётся от другого биннинга: множество
    неопределённости обязано совпадать с тем, что измерено, иначе минимакс защищает не от того.
    """
    edges = np.quantile(s_oof, np.linspace(0, 1, nstrat + 1)[1:-1])
    k_oof = np.digitize(s_oof, edges)
    k_te = np.digitize(s_te, edges)
    p = np.bincount(k_oof, minlength=nstrat) / len(k_oof)
    q = np.bincount(k_te, minlength=nstrat) / len(k_te)
    rho = float(((q - p) ** 2 / np.maximum(p, 1e-12)).sum())
    return k_oof, p, q, rho, edges


def dro_weights(loss, k, p, rho):
    """Худший случай внутри chi2-шара, замкнутая форма. Возвращает вес на строку."""
    L = np.array([loss[k == j].mean() if (k == j).any() else 0.0 for j in range(len(p))])
    mu = float((p * L).sum())
    var = float((p * (L - mu) ** 2).sum())
    if var <= 1e-12 or rho <= 0:
        w_k = np.ones(len(p))
    else:
        w_k = np.clip(1.0 + np.sqrt(rho / var) * (L - mu), 0.0, None)
    return w_k[k], w_k


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--rhos", default="0,0.25,0.5,1.0,measured")
    ap.add_argument("--out", default=RES + "preds/oof_dro.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]
    fte = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in te.SMILES
           if Chem.MolFromSmiles(s) is not None]
    rng = np.random.default_rng(0)
    n_sub = int(round(0.8 * len(fps)))
    s_te = np.zeros(len(fte))
    for _ in range(5):
        ref = [fps[j] for j in rng.choice(len(fps), n_sub, replace=False)]
        for t, f_ in enumerate(fte):
            s_te[t] += max(DataStructs.BulkTanimotoSimilarity(f_, ref)) / 5

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        s_oof = nn_similarity(fps, fold)
        k_all, p, q, rho_m, edges = strata_and_radius(s_oof, s_te)
        print(f"сид {seed}: страты по квантилям, границы "
              + " ".join(f"{e:.3f}" for e in edges))
        print(f"  доли OOF  " + " ".join(f"{v:.3f}" for v in p))
        print(f"  доли тест " + " ".join(f"{v:.3f}" for v in q))
        print(f"  измеренный радиус для ЭТОГО разбиения: rho = {rho_m:.3f}", flush=True)

        # Веса под режим теста, для второй половины предрегистрации.
        h_o = np.histogram(s_oof, BINS)[0] / len(s_oof)
        h_t = np.histogram(s_te, BINS)[0] / len(s_te)
        ratio = np.where(h_o > 0, h_t / np.maximum(h_o, 1e-9), 1.0).clip(0, 8.0)
        w_test = ratio[np.clip(np.digitize(s_oof, BINS) - 1, 0, len(ratio) - 1)]

        rhos = [rho_m if r == "measured" else float(r) for r in a.rhos.split(",")]
        for rho in rhos:
            t0 = time.time()
            r = {"seed": seed, "rho": round(float(rho), 3)}
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                Xi, fi, ki = X[m], fold[m], k_all[m]

                # Шаг 1: обычная подгонка, чтобы измерить потерю по стратам.
                p0 = np.zeros_like(y)
                for f in range(5):
                    trn, tei = fi != f, fi == f
                    if tei.sum() == 0:
                        continue
                    p0[tei] = HistGradientBoostingRegressor(**KW).fit(
                        Xi[trn], y[trn]).predict(Xi[tei])
                num = np.maximum(0.0, np.maximum(lo - p0, p0 - hi))
                w, w_k = dro_weights(num, ki, p, rho)

                # Шаг 2: переподгонка под весами худшего случая.
                p1 = np.zeros_like(y)
                for f in range(5):
                    trn, tei = fi != f, fi == f
                    if tei.sum() == 0:
                        continue
                    p1[tei] = HistGradientBoostingRegressor(**KW).fit(
                        Xi[trn], y[trn], sample_weight=w[trn]).predict(Xi[tei])

                u = np.ones(len(y)) / len(y)
                qq = fit_apply(p1, lo, hi, fi, u)
                out[f"{seed}|{rho:.3f}|{c}"] = p1.tolist()
                r[f"{c} пара"] = round(float(strae(y, qq, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p1).statistic), 4)
                # Ранг под режим теста: взвешенный Спирмен, как в пункте 113.
                ww = w_test[m]
                ra, rb = rankdata(y).astype(float), rankdata(p1).astype(float)
                sw = ww.sum()
                ra -= (ww * ra).sum() / sw
                rb -= (ww * rb).sum() / sw
                r[f"{c} rho_test"] = round(float(
                    (ww * ra * rb).sum() / np.sqrt((ww * ra * ra).sum() * (ww * rb * rb).sum())), 4)
                if c == CYPS[0]:
                    r["веса страт"] = " ".join(f"{v:.2f}" for v in w_k)
            for tag in ("пара", "rho", "rho_test"):
                r[f"MACRO {tag}"] = round(float(np.mean([r[f"{c} {tag}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  rho {rho:5.3f}  пара {r['MACRO пара']:.4f}  ранг {r['MACRO rho']:.4f}  "
                  f"ранг-под-тест {r['MACRO rho_test']:.4f}   веса {r['веса страт']}  "
                  f"({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("rho")[["MACRO пара", "MACRO rho", "MACRO rho_test"]].mean().to_string())
    print()
    print("Только 2C9 и 3A4 (ESS 36.7 и 39.3 %):")
    print(df.groupby("rho")[["CYP2C9 rho_test", "CYP3A4 rho_test"]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Строка rho = 0 обязана воспроизвести базу (0.7150 / 0.5651 на сиде 0): при нулевом
радиусе веса единичные и это в точности обычная подгонка. Расхождение там --- ошибка проводки,
и остальные строки тогда не значат ничего.

Предрегистрация двусторонняя. DRO ОБЯЗАН проиграть по равномерной метрике: он покупает
робастность ценой среднего случая, и рука, выигрывающая по обеим, делает не то, что заявляет.
И ОБЯЗАН выиграть по ранг-под-тест --- это единственная имеющаяся оценка того режима, ради
которого он строится. Проигрыш по обеим --- недообучение, и арка закрыта.

Судить по 2C9 и 3A4. На 1A2 и 2D6 сам критерий вчетверо шумнее (пункт 123 намерил ESS 24.8 и
20.7 %), и там нужны четыре сида, прежде чем числа что-то значат.""")


if __name__ == "__main__":
    main()
