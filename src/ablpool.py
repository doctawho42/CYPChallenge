"""One model over all four enzymes at once, instead of four models that never speak.

The gap. Every model in this repository is fitted per enzyme on that enzyme's own labels, so
CYP2D6 sees 1493 rows and nothing else. But the four label sets sit on overlapping molecules and
share most of their chemistry: lipophilicity, size, reactive motifs and steric bulk mean the same
thing to all four enzymes, and only the recognition motif differs. Pooling the four into one
table with an enzyme indicator lets the shared part be learned on 6525 rows and the specific part
through interactions with the indicator. It is the same move the source-indicator column made for
the external data, applied to a structure that was here all along.

Why it is worth a run under item 80's criterion, when five other ideas were not. Pooling changes
what the model can learn rather than how its output is scaled, so it is the kind of intervention
that can move rank. Falsification is therefore sharp and pre-declared: if out-of-fold Spearman
does not rise on any enzyme, the idea failed, whatever the raw ST-RAE says.

Three arms:

  independent   - four models, one per enzyme. The control, which must reproduce 0.7673;
  pooled        - one model, 6525 rows, four one-hot enzyme columns;
  pooled + TDI  - the same, plus every compound's pre-incubation pIC50 as an extra row with a
                  condition flag, which doubles the table to about 13000. The two conditions
                  correlate 0.906 to 0.990 as values, so as a target the second adds little, but
                  as extra supervision on the shared trunk of the problem it is not the same
                  thing, and the residual between them is 17 to 45 percent of the signal.

Folds are the Butina folds of the compound, so every copy of a molecule - four enzymes, two
conditions - lands in the same fold and nothing leaks between them. Scoring is on held-out
direct-inhibition labels only, per enzyme, exactly as everywhere else.

Reads data/feats.npz, data/rows.csv and the TDI file. Writes results/preds/oof_pool.json.
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
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--arms", default="независимо,пул,пул+TDI")
    ap.add_argument("--blocks", default="FP+DESC+MECH",
                    help="какие блоки признаков; FP+DESC отключает механистический")
    ap.add_argument("--out", default=RES + "preds/oof_pool.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    arms = a.arms.split(",")

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    td = (pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv")
            .set_index("Molecule_Name").reindex(rows.Molecule_Name).reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z[b] for b in a.blocks.split("+")]).astype(np.float32)
    n = len(X)
    print(f"признаки: {a.blocks}, ширина {X.shape[1]}")

    # Маски и метки: прямое ингибирование --- то, на чём считаем; TDI --- только надзор.
    M = {c: tr[f"{c}_pIC50_direct_inhibition"].notna().to_numpy() for c in CYPS}
    Y = {c: tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS}
    MT = {c: td[f"{c}_pIC50_TDI_condition"].notna().to_numpy() for c in CYPS}
    YT = {c: td[f"{c}_pIC50_TDI_condition"].to_numpy(float) for c in CYPS}
    print("строк: " + ", ".join(f"{c} {int(M[c].sum())}" for c in CYPS)
          + f"; всего прямых {sum(int(M[c].sum()) for c in CYPS)}"
          + f", с TDI {sum(int(M[c].sum()) + int(MT[c].sum()) for c in CYPS)}\n")

    def stack(keep, with_tdi):
        """Собрать пулированную таблицу из строк, помеченных keep (по молекулам)."""
        Xs, ys = [], []
        for e, c in enumerate(CYPS):
            for cond, (mask, lab) in enumerate(((M[c], Y[c]), (MT[c], YT[c]))):
                if cond and not with_tdi:
                    continue
                sel = mask & keep
                if not sel.any():
                    continue
                ind = np.zeros((int(sel.sum()), 5), np.float32)
                ind[:, e] = 1.0
                ind[:, 4] = cond
                Xs.append(np.hstack([X[sel], ind]))
                ys.append(lab[sel])
        return np.vstack(Xs), np.concatenate(ys)

    def block(e, sel):
        """Признаки для предсказания: тот же индикатор, условие --- прямое ингибирование."""
        ind = np.zeros((int(sel.sum()), 5), np.float32)
        ind[:, e] = 1.0
        return np.hstack([X[sel], ind])

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for arm in arms:
            t0 = time.time()
            r = {"seed": seed, "рука": arm}
            for e, c in enumerate(CYPS):
                m = M[c]
                y, fi = Y[c][m], fold[m]
                lo = tr.loc[m, f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy()
                hi = tr.loc[m, f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy()
                p = np.zeros_like(y)
                for f in range(5):
                    te_mol = fold == f
                    te = fi == f
                    if te.sum() == 0:
                        continue
                    if arm == "независимо":
                        Xt, yt = X[m][~te], y[~te]
                        Xe = X[m][te]
                    else:
                        Xt, yt = stack(~te_mol, arm.endswith("TDI"))
                        Xe = block(e, m & te_mol)
                    p[te] = HistGradientBoostingRegressor(**KW).fit(Xt, yt).predict(Xe)
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{arm}|{c}"] = p.tolist()
                r[f"{c} сырой"] = round(float(strae(y, p, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            r["MACRO сырой"] = round(float(np.mean([r[f"{c} сырой"] for c in CYPS])), 4)
            r["MACRO пара"] = round(float(np.mean([r[f"{c} пара"] for c in CYPS])), 4)
            r["MACRO rho"] = round(float(np.mean([r[f"{c} rho"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {arm:10s} сырой {r['MACRO сырой']:.4f}  пара {r['MACRO пара']:.4f}"
                  f"  rho {r['MACRO rho']:.4f}  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("рука")[["MACRO сырой", "MACRO пара", "MACRO rho"]].mean().to_string())
    print()
    print(df.groupby("рука")[[f"{c} rho" for c in CYPS]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Первая рука обязана дать 0.7673 сырым на сиде 0, иначе сравнивать не с чем.

Дальше смотреть на rho, а не на ST-RAE. Пул меняет то, что модель СПОСОБНА выучить, а не то,
как отмасштабирован её выход, поэтому он относится к тому классу вмешательств, которые
постобработка не съедает. Если ранг не вырос ни на одном ферменте --- идея не сработала, и
сырое ST-RAE тут ничего не спасает.

Ждём прироста прежде всего на CYP2D6: у него меньше всех своих меток относительно того, что
даёт пул, и худший ранг из четырёх.""")


if __name__ == "__main__":
    main()
