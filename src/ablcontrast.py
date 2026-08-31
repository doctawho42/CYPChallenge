"""The enzyme contrast from the screen, predicted from structure and handed to the learner.

Where the four-vector is. The curves cannot supply between-enzyme comparison: 3596 of 4905
compounds carry exactly one enzyme's pIC50 and only **41** carry all four, so a model fitted on
curves almost never sees the same molecule on two enzymes. The single-concentration screen is the
one place in the dataset where it does -- 4376 molecules, all four enzymes on every one of them,
one protocol, no gaps.

Level versus contrast, and why the distinction is the whole idea. The screen's *level* (the mean
of a molecule's four readings) is a proxy for general potency, which the pIC50 head already learns
from 6525 labels; handing it over again duplicates what is there. The *contrast*
(`l2fc_e - mean_e l2fc`, four numbers summing to zero) is the part that says which enzyme this
molecule prefers, and it is exactly what the curves cannot give. It is also the larger component
here, which was not obvious: measured over the 4375 complete rows the contrast has a standard
deviation of 0.914 against the level's 0.673.

Predicted from structure, never read off the screen. The test set has zero overlap with the
primary library -- checked, zero by name and zero by SMILES -- so a screening reading will not
exist for a test compound. The contrast therefore enters as a *prediction* from structure, fitted
out of fold on the training folds' screen and applied to the held-out fold from its structure
alone. That is what deployment will look like, so it is what is measured.

Four arms, with the control that makes the claim falsifiable:

    база                 FP+DESC+MECH
    + контраст           four predicted contrasts. The proposal.
    + уровень            one predicted level. The control -- if this does just as well, the
                         distinction between level and contrast is decoration and what helps is
                         simply "some screening information".
    + оба                both, to see whether they are additive.

Pre-registered reading, and one conflict inside it that is recorded before the run finishes. The
gain should be larger where the curve labels are most single-enzyme, because that is where the
contrast is least recoverable from the labels themselves. Measured, those shares are CYP1A2 45.9,
CYP2C9 43.4, CYP2D6 54.6 and CYP3A4 67.5 per cent, which points at CYP3A4.

But item 129 found CYP3A4's single-enzyme excess *is* the 530-compound analog campaign, and those
compounds are outside the screening library entirely -- their contrast is filled with zero because
there is nothing to predict it from. So the enzyme the prediction favours is the one where 23 per
cent of the rows cannot receive the channel at all. On CYP3A4 the prediction therefore applies to
its 1805 screen compounds and not to the enzyme as a whole, and a flat result there is not evidence
against the channel.

Either way it must survive the affine pair, since this is a *feature* and not a post-hoc
correction -- item 128 has since bounded that whole class at 0.0076, which is a further reason a
feature is the right shape for this.

Same folds, masks, metric and learner settings as src/ablate.py. Writes
results/preds/oof_contrast.json.
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
KW_SCR = dict(max_iter=200, learning_rate=0.06, max_leaf_nodes=31,
              l2_regularization=1.0, random_state=0)


def screen_channels(X, S, fold):
    """Out-of-fold predictions of the screen's contrast and level, from structure alone.

    S is (n, 4) of log2fc with NaN where the molecule is not in the primary library. A fold's
    channels are fitted on the other folds' screen only; nothing about a held-out molecule's own
    reading is used, because at test time there is none.
    """
    n = len(X)
    lvl = np.full(n, np.nan)
    con = np.full((n, 4), np.nan)
    have = ~np.isnan(S).any(axis=1)
    # Считать среднее только там, где есть все четыре: np.nanmean по пустой строке
    # предупреждает и всё равно даёт NaN, так что маска дешевле и честнее.
    level = np.full(len(S), np.nan)
    level[have] = S[have].mean(axis=1)
    contrast = S - level[:, None]
    for f in np.unique(fold):
        te, trn = fold == f, (fold != f) & have
        if trn.sum() < 50 or te.sum() == 0:
            continue
        lvl[te] = HistGradientBoostingRegressor(**KW_SCR).fit(
            X[trn], level[trn]).predict(X[te])
        for e in range(4):
            con[te, e] = HistGradientBoostingRegressor(**KW_SCR).fit(
                X[trn], contrast[trn, e]).predict(X[te])
    return con, lvl, have


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--out", default=RES + "preds/oof_contrast.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    piv = sc.pivot_table(index="Molecule_Name", columns="enzyme", values="log2fc_estimate")
    S = rows.set_index("Molecule_Name").join(piv)[CYPS].to_numpy(float)

    M = {c: tr[f"{c}_pIC50_direct_inhibition"].notna().to_numpy() for c in CYPS}
    k = sum(M[c].astype(int) for c in CYPS)
    print("Насколько односоставны метки (доля соединений фермента, у которых он единственный):")
    only = {}
    for c in CYPS:
        only[c] = float((k[M[c]] == 1).mean())
        print(f"  {c}: {100*only[c]:.1f} %")
    print("\nПредсказание: прирост больше там, где эта доля выше.\n")

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        t0 = time.time()
        con, lvl, have = screen_channels(X, S, fold)
        print(f"  сид {seed}: каналы скрина посчитаны на {int(have.sum())} строках "
              f"({time.time()-t0:.0f} с)", flush=True)
        # Пропуски заполняются нулём: ноль --- это «контраста нет», нейтральное значение,
        # а не выдуманное показание.
        conf = np.nan_to_num(con, nan=0.0)
        lvlf = np.nan_to_num(lvl, nan=0.0)[:, None]
        arms = {"база": X,
                "+ контраст": np.hstack([X, conf.astype(np.float32)]),
                "+ уровень": np.hstack([X, lvlf.astype(np.float32)]),
                "+ оба": np.hstack([X, conf.astype(np.float32), lvlf.astype(np.float32)])}

        for nm, XX in arms.items():
            t1 = time.time()
            r = {"seed": seed, "рука": nm}
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = M[c]
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                Xi, fi = XX[m], fold[m]
                p = np.zeros_like(y)
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    p[te] = HistGradientBoostingRegressor(**KW).fit(
                        Xi[trn], y[trn]).predict(Xi[te])
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{nm}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for tag in ("пара", "rho"):
                r[f"MACRO {tag}"] = round(float(np.mean([r[f"{c} {tag}"] for c in CYPS])), 4)
            table.append(r)
            print(f"    {nm:12s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t1:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                         + [f"{c} rho" for c in CYPS]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Контроль «+ уровень» решает, есть ли вообще различие между уровнем и контрастом.
Если он даёт столько же, значит помогает просто «какая-нибудь информация из скрина», и вся
конструкция про внутримолекулярный контраст --- украшение.

Прирост обязан быть неравномерным и следовать за односоставностью меток, напечатанной выше.
Равномерный прирост означал бы, что работает не контраст, а четыре лишние колонки.

И он обязан пережить аффинную пару: это ПРИЗНАК, а не поправка. Класс поправок закрыт
пунктом 128 с бюджетом 0.0076, поэтому признак --- единственная форма, в которой этот канал
вообще имеет шанс.""")


if __name__ == "__main__":
    main()
