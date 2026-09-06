"""What could a better pKa estimator buy? Bound it by injecting the rule's own error.

The mechanistic block is the best feature block in the project -- +0.0163 of rank, +0.0313 in the
test regime, and three times less sensitive to the split than anything else -- and it rests on a
rule that estimates pKa. An outside reading proposed replacing that rule with a trained predictor,
observing that the most basic centre is an imidazole or pyridine for 55.8 per cent of the training
set and that the rule's known worst miss, caffeine, is 4.5 pKa units.

**The exposure is real.** pKa enters three features -- the value, the protonated fraction at pH 7.4,
and a hard indicator above 7.4 -- plus the CYP2D6 pharmacophore, which asks for a cation at pH 7.4
two to five bonds from an aryl. So an error matters exactly where it moves a compound across 7.4.
Measured: **13.1 per cent** of the training set sits within one pKa unit of the threshold and 23.5
per cent within two, against a typical predictor error of 0.5 to 1.0.

**But exposure is not the question, and this file asks the question instead of building the
answer.** Whether a better estimate helps depends on whether the current one is wrong *and* whether
the block is sensitive to being wrong. Both are settled at once by injecting noise of the rule's own
error magnitude and measuring what it costs:

    если шум величиной с ошибку правила НЕ вредит  -> блок к pKa нечувствителен,
                                                      лучший предсказатель не купит ничего
    если вредит на X                              -> X есть ВЕРХНЯЯ ГРАНИЦА выигрыша,
                                                      потому что снять ошибку --- самое большее,
                                                      что предсказатель может сделать

That is item 128's procedure: bound a class of proposals with one measurement before building any
member of it. Here the bound is unusually tight, because injecting error and removing error are the
same operation with opposite sign.

**The three derived features are recomputed from the perturbed pKa**, not merely jittered -- the
protonated fraction is a sigmoid of `pKa - 7.4` and the indicator is a hard threshold on it, so a
perturbation of the value propagates non-linearly and the hard feature flips or does not. Perturbing
the outputs directly would measure something else.

Reads data/feats.npz, data/rows.csv, data/mech_names.csv. Same folds, learner and metric as
src/ablate.py.
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
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--sigmas", default="0,1.0")
    ap.add_argument("--out", default=RES + "preds/oof_pkanoise.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    sigmas = [float(x) for x in a.sigmas.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    FP, DESC, MECH = z["FP"], z["DESC"], z["MECH"].astype(np.float64)
    nm = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    i_pka, i_frac, i_base = (nm.index("pka_max_basic"), nm.index("frac_prot_74"),
                             nm.index("is_base_74"))
    p0 = MECH[:, i_pka].copy()
    near = np.abs(p0 - 7.4) < 1.0
    print(f"pKa: медиана {np.median(p0):.2f}, в пределах 1 единицы от 7.4 --- "
          f"{near.sum()} молекул ({near.mean():.1%})")
    print(f"is_base_74 включён у {int(MECH[:, i_base].sum())} ({MECH[:, i_base].mean():.1%})\n",
          flush=True)

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for sg in sigmas:
            t0 = time.time()
            M = MECH.copy()
            if sg > 0:
                rng = np.random.default_rng(9000 + seed)
                # Шум только там, где основный центр вообще есть: у молекул без него
                # правило не ошибается, ему нечего оценивать.
                has = p0 > 0.5
                pn = p0.copy()
                pn[has] = np.maximum(p0[has] + rng.normal(0, sg, has.sum()), 0.0)
                M[:, i_pka] = pn
                # Производные признаки ПЕРЕСЧИТЫВАЮТСЯ, а не дрожат: доля протонирования ---
                # сигмоида от pKa-7.4, индикатор --- жёсткий порог, и они переворачиваются.
                M[:, i_frac] = 1.0 / (1.0 + np.power(10.0, 7.4 - pn))
                M[:, i_base] = (pn > 7.4).astype(float)
                flipped = int((M[:, i_base] != MECH[:, i_base]).sum())
                print(f"  сигма {sg}: перевернулось is_base_74 у {flipped} молекул "
                      f"({flipped/len(p0):.1%})", flush=True)
            X = np.hstack([FP, DESC, M.astype(np.float32)])
            r = {"seed": seed, "сигма": sg}
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
                    p[te] = HistGradientBoostingRegressor(**KW).fit(
                        Xi[trn], y[trn]).predict(Xi[te])
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{sg}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for t in ("пара", "rho"):
                r[f"MACRO {t}"] = round(float(np.mean([r[f"{c} {t}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} сигма {sg:<4} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("сигма", sort=False)[["MACRO пара", "MACRO rho"]
                                          + [f"{c} rho" for c in CYPS]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Величина падения при сигме 1.0 --- это ВЕРХНЯЯ ГРАНИЦА того, что может купить
обученный предсказатель pKa, потому что убрать ошибку --- самое большее, что он делает, а
внести её --- то, что здесь сделано.

Смотреть надо прежде всего на CYP2D6: блок живёт там (+0.0454 из +0.0163 макро), и именно там
фармакофор спрашивает про катион при pH 7.4. Падение меньше поферментного пола 0.0049 на CYP2D6
означает, что заменять правило незачем, каким бы неточным оно ни было.""")


if __name__ == "__main__":
    main()
