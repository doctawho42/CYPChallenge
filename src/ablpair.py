"""Which partner enzyme does the pooling gain come from, and does contrast beat similarity.

Two mechanisms for pooling have been proposed and both are dead. Item 110 refuted borrowing
neighbours: pooling wins most where the per-enzyme model already sees the whole neighbourhood.
Item 111 refuted transfer of shared function, and refuted it with a perfect inversion --- the
enzyme least correlated with the others, CYP2D6 at a mean rank correlation of -0.002, takes the
largest gain at +0.0374, and the most correlated, CYP3A4 at 0.375, is the only one pooling
actually hurts.

The inversion suggests a third mechanism: what the pooled model gets is not similarity but
**contrast**. With the enzyme indicator in the design a tree can learn "this split matters for
CYP2D6 and not for CYP3A4", which is strictly more than "this split matters for CYP2D6", and a
per-enzyme model cannot represent it at all. On that reading a partner is useful in proportion
to how much its structure-activity relationship *differs*, and a partner that behaves almost
identically --- CYP2C9 and CYP3A4 correlate at 0.688 --- adds rows without adding information.

Item 111 cannot separate this from plain sample size, because across the four enzymes the
number of labels and the correlation with the rest move together: CYP3A4 has both the most rows
and the highest correlation. Four points cannot split two explanations that are themselves
correlated.

Pooling one partner at a time does split them. For a fixed target enzyme the three partners
contribute different numbers of rows *and* different correlations, and those two orderings do
not agree, so the gain can be regressed on both.

Predictions, stated before the run. Under contrast: CYP2D6 gains most from CYP2C9, the partner
it anti-correlates with at -0.120, and CYP3A4 gains least or loses from CYP2C9. Under sample
size: every target gains most from CYP3A4, which brings the most rows, and the ordering is the
same for all four targets. The two make different predictions for eight of the twelve cells.

Same folds, masks, metric and learner settings as src/ablpool.py, and the same design matrix ---
the enzyme indicator is kept at its full width so the arms differ only in which rows are present.
Writes results/preds/oof_pair.json.
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1")
    ap.add_argument("--out", default=RES + "preds/oof_pair.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)

    M = {c: tr[f"{c}_pIC50_direct_inhibition"].notna().to_numpy() for c in CYPS}
    Y = {c: tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS}

    print("Что каждый партнёр приносит целевому ферменту:")
    print(f"{'цель':8s} {'партнёр':8s} {'строк':>7s} {'x к своим':>10s} {'корреляция':>11s}")
    corr = {}
    for c in CYPS:
        for d in CYPS:
            if c == d:
                continue
            both = M[c] & M[d]
            r = float(spearmanr(Y[c][both], Y[d][both]).statistic)
            corr[(c, d)] = r
            print(f"{c:8s} {d:8s} {int(M[d].sum()):7d} {M[d].sum()/M[c].sum():10.2f} {r:11.3f}")
    print()

    def design(sel, e):
        """Строки sel, помеченные индикатором фермента e. Ширина индикатора всегда пятёрка,
        чтобы руки отличались только составом строк, а не формой матрицы."""
        ind = np.zeros((int(sel.sum()), 5), np.float32)
        ind[:, e] = 1.0
        return np.hstack([X[sel], ind])

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for e, c in enumerate(CYPS):
            lo = tr.loc[M[c], f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy()
            hi = tr.loc[M[c], f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy()
            y, fi = Y[c][M[c]], fold[M[c]]
            arms = {"один": [c]}
            for d in CYPS:
                if d != c:
                    arms[f"+{d[3:]}"] = [c, d]
            arms["+все"] = list(CYPS)

            for nm, members in arms.items():
                t0 = time.time()
                p = np.zeros_like(y)
                for f in range(5):
                    keep = fold != f
                    Xs, ys = [], []
                    for e2, c2 in enumerate(CYPS):
                        if c2 not in members:
                            continue
                        sel = M[c2] & keep
                        if not sel.any():
                            continue
                        Xs.append(design(sel, e2)); ys.append(Y[c2][sel])
                    te = M[c] & (fold == f)
                    if te.sum() == 0:
                        continue
                    p[fi == f] = HistGradientBoostingRegressor(**KW).fit(
                        np.vstack(Xs), np.concatenate(ys)).predict(design(te, e))
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{c}|{nm}"] = p.tolist()
                table.append({"seed": seed, "цель": c, "рука": nm,
                              "пара": float(strae(y, q, y_true_upper=hi, y_true_lower=lo)),
                              "ранг": float(spearmanr(y, p).statistic)})
                print(f"  сид {seed} {c} {nm:8s} пара {table[-1]['пара']:.4f} "
                      f"ранг {table[-1]['ранг']:.4f}  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print("\nПрирост ранга против одиночного обучения, среднее по сидам:")
    piv = df.pivot_table(index="цель", columns="рука", values="ранг", aggfunc="mean")
    gain = piv.sub(piv["один"], axis=0).drop(columns=["один"])
    print(gain.to_string(float_format=lambda x: f"{x:+.4f}"))

    print("\nСвязь прироста с корреляцией и с числом добавленных строк (12 клеток):")
    g, rr, nn = [], [], []
    for c in CYPS:
        for d in CYPS:
            if c == d:
                continue
            g.append(float(gain.loc[c, f"+{d[3:]}"])); rr.append(corr[(c, d)])
            nn.append(float(M[d].sum()) / float(M[c].sum()))
    print(f"  с корреляцией     rho = {spearmanr(g, rr).statistic:+.3f}")
    print(f"  с объёмом добавки rho = {spearmanr(g, nn).statistic:+.3f}")
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Отрицательная связь с корреляцией и слабая с объёмом --- контраст; сильная связь
с объёмом и никакой с корреляцией --- обычная регуляризация числом строк, и тогда индикатор
ни при чём. Обе сильные --- разделить не удалось, и это тоже честный исход при двенадцати
клетках.

Столбец «+все» проверяет, складываются ли партнёры: если он примерно равен лучшему одиночному
партнёру, то работает один из трёх, а не все вместе.""")


if __name__ == "__main__":
    main()
