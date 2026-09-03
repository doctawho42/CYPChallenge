"""logD at assay pH, and lipophilic efficiency as a target — both from columns already present.

Two cheap things the descriptor block does not carry, and one of them is embarrassing.

**logD is missing and logP is not.** The matrix has Crippen MolLogP, twelve SlogP_VSA bins and two
BCUT logP eigenvalues, all of which describe the *neutral* molecule. The assay runs at pH 7.4, and
for a basic amine the difference between logP and logD there is two to three units, because most
of the molecule is protonated and does not partition. CYP2D6 is the amine enzyme — its
pharmacophore is a protonated basic nitrogen — so this is the one enzyme where the distinction
should matter most, and it is the one we predict worst.

Nothing new is computed. Both inputs are already in the matrix:

    logD = logP - log10(1 + 10^(pKa_basic - 7.4))

`MolLogP` is DESC column 130 and `pka_max_basic` is MECH column 16. A tree can in principle
reconstruct this from the two of them, and the question is whether it does — a monotone
combination of two columns is not something a depth-limited tree finds for free.

**Lipophilic efficiency as the target.** CYP inhibition, and CYP3A4 inhibition in particular, is
substantially explained by lipophilicity: bigger, greasier molecules stick to a big greasy pocket.
Fitting pIC50 spends capacity on that easy axis. Fitting `pIC50 - logP` instead and adding the
term back afterwards leaves the model only the specific part. The prediction is identical in
expectation and the *loss geometry* is not: an error in the specific part stops being an order of
magnitude cheaper than an error in the bulk part.

This is the arm to be most suspicious of, and the suspicion is pre-registered. Adding back
`c * logP` is a monotone function of a feature the model already has, so item 128's oracle gate
applies in spirit: if the whole effect is that the target became easier, the affine pair will take
it back. What would survive is a change in *which* compounds the model gets right, so this is
judged on rank and per-enzyme, not on the pair.

Arms:

    база              FP+DESC+MECH
    + logD            one extra column
    мишень LipE       target pIC50 - logP, prediction plus logP
    + logD и LipE     both

Pre-registered reading. logD should help CYP2D6 first, since that is the amine enzyme and the one
whose labels sit at the bottom of our table. LipE should help CYP3A4 first, since that is where
lipophilicity explains the most and therefore where the most capacity is being spent on it. A gain
that is uniform across enzymes is neither, and means an extra column helped a boosting.

Same folds, masks, metric and learner settings as src/ablate.py. Writes
results/preds/oof_lipo.json.
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
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
PH = 7.4


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--out", default=RES + "preds/oof_lipo.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    FP, DESC, MECH = z["FP"], z["DESC"], z["MECH"]
    X = np.hstack([FP, DESC, MECH]).astype(np.float32)

    dn = [l.strip() for l in open(D + "desc_names.csv")]
    mn = [l.strip() for l in open(D + "mech_names.csv")]
    logp = np.nan_to_num(DESC[:, dn.index("MolLogP")].astype(np.float64))
    pka = np.nan_to_num(MECH[:, mn.index("pka_max_basic")].astype(np.float64))
    # Только там, где основание вообще есть: pKa = 0 означает «нет основного центра»,
    # и подставлять его в формулу значило бы вычитать log10(1 + 10^-7.4), то есть ноль ---
    # но лучше сказать это явно, чем полагаться на то, что оно само занулится.
    has_base = pka > 0.5
    logd = logp - np.where(has_base, np.log10(1.0 + 10.0 ** (pka - PH)), 0.0)
    print(f"основной центр есть у {int(has_base.sum())} из {len(pka)} "
          f"({100*has_base.mean():.1f} %)")
    print(f"logP медиана {np.median(logp):.2f}; logD медиана {np.median(logd):.2f}; "
          f"сдвиг у оснований медиана {np.median((logp-logd)[has_base]):.2f}\n")

    XD = np.hstack([X, logd[:, None].astype(np.float32)])
    ARMS = {"база": (X, False), "+ logD": (XD, False),
            "мишень LipE": (X, True), "+ logD и LipE": (XD, True)}

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for nm, (XX, lipe) in ARMS.items():
            t0 = time.time()
            r = {"seed": seed, "рука": nm}
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                Xi, fi, lp = XX[m], fold[m], logp[m]
                # Мишень LipE: вычесть липофильность, обучиться, вернуть член обратно.
                tgt = y - lp if lipe else y
                p = np.zeros_like(y)
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    p[te] = HistGradientBoostingRegressor(**KW).fit(
                        Xi[trn], tgt[trn]).predict(Xi[te])
                if lipe:
                    p = p + lp
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{nm}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for tag in ("пара", "rho"):
                r[f"MACRO {tag}"] = round(float(np.mean([r[f"{c} {tag}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {nm:16s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                         + [f"{c} rho" for c in CYPS]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Первая рука обязана дать 0.7150 / 0.5651 на сиде 0.

logD ждём на 2D6: это фермент аминов, и разница logP с logD у оснований как раз там наибольшая.
LipE ждём на 3A4: там липофильность объясняет больше всего, значит там на неё и тратится
больше всего ёмкости. Равномерный прирост не подтверждает ни то ни другое.

LipE судить по рангу и поферментно, а не по паре: прибавление c*logP --- монотонная функция
уже имеющегося признака, и если весь эффект в том, что мишень стала легче, аффинная пара его
заберёт. Пережить может только смена того, КАКИЕ соединения модель угадывает.""")


if __name__ == "__main__":
    main()
