"""Multi-task trees: one split structure for all four enzymes, four values in every leaf.

The defect this is built against, measured rather than supposed. Item 180 found that pooling
reverses on plain depth-5 trees and recovers at depth 8, and said why: the share of root-to-leaf
paths passing through the enzyme indicator goes from 0.34 to 0.59. **A pooled model has to spend
depth isolating the enzyme before it can model anything conditional on it**, and five levels do not
leave enough. More trees do not substitute -- 400 and 800 recover a quarter of the gap -- because
the constraint is per-tree expressiveness.

A multi-task tree removes the cost rather than paying it. The enzyme is not a column and is never
split on; the tree is grown once on all 4905 molecules against a **four-column** residual, and each
leaf carries four numbers. Structure is shared across every labelled row in the table; specificity
lives in the leaf and costs no depth at all.

    пул            split(x) ... split(индикатор) ... split(x)  ->  одно число
    многозадачно   split(x) ... split(x)                       ->  вектор из четырёх

`DecisionTreeRegressor` does this natively: given a 2-D target it grows one tree minimising the
summed impurity across outputs. Nothing exotic is required and the constants stay pinned to
`src/ablpairloss.py`.

**The missing labels are the honest difficulty and are handled explicitly.** The label matrix is
sparse -- 1412, 1285, 1493 and 2335 of 4905 -- and a multi-output tree cannot take NaN. Unobserved
cells get a residual of exactly zero, so they contribute no gradient, but they do enter the split
criterion as zeros and therefore pull it toward not splitting where data is absent. That is a real
bias and it is why the second arm exists: scaling each column by the inverse spread of its observed
residuals stops the enzyme with the widest residuals from dominating the shared structure, which on
this table is CYP2D6.

Arms:

    независимо              четыре модели, эталон
    пул                     индикатор фермента колонкой, глубина 5 --- рука, которая проигрывает
    многозадачно            одна структура, четыре выхода
    многозадачно, масштаб   то же, колонки residual-ов нормированы на свой разброс

**Pre-registered.** If item 180's diagnosis is right, the multi-task arm must beat `независимо` at
depth 5 -- it gets the enzyme conditioning that pooling pays depth for, at no depth. If it merely
matches pooling's loss, the diagnosis is wrong and the problem is not depth. Read per enzyme against
item 165's floors, and expect the gain to be largest where labels are fewest (CYP2C9, 1285) since
shared structure is worth most where a per-enzyme model has least to fit on.

Writes results/preds/oof_multi.json.
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
from sklearn.tree import DecisionTreeRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
NTREE, LR, DEPTH, MAXFEAT = 200, 0.06, 5, 0.3      # пины из src/ablpairloss.py
FLOOR = {"CYP1A2": 0.0061, "CYP2C9": 0.0071, "CYP2D6": 0.0049, "CYP3A4": 0.0033}


def boost1(Xtr, ytr, Xte, seed):
    """Обычный одновыходный бустинг --- для независимой и пулированной рук."""
    rng = np.random.default_rng(seed)
    s = np.full(len(ytr), float(ytr.mean()))
    pred = np.full(len(Xte), s[0])
    for _ in range(NTREE):
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MAXFEAT,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(Xtr, ytr - s)
        s = s + LR * t.predict(Xtr)
        pred = pred + LR * t.predict(Xte)
    return pred


def boostM(Xtr, Ytr, Mtr, Xte, seed, scale):
    """Многозадачный: одна структура, четыре выхода.

    Ytr --- (n, 4) с нулями там, где метки нет; Mtr --- маска наблюдённого. Невидимые
    ячейки держатся на текущем предсказании, поэтому их residual тождественно ноль и
    градиента они не дают; в критерий сплита они всё же входят нулями, и это записано
    в докстринге как смещение, а не спрятано.
    """
    rng = np.random.default_rng(seed)
    base = np.array([Ytr[Mtr[:, e], e].mean() if Mtr[:, e].any() else 0.0
                     for e in range(Ytr.shape[1])])
    s = np.tile(base, (len(Ytr), 1))
    pred = np.tile(base, (len(Xte), 1))
    w = np.ones(Ytr.shape[1])
    if scale:
        for e in range(Ytr.shape[1]):
            v = Ytr[Mtr[:, e], e]
            w[e] = 1.0 / max(v.std(), 1e-6)
    for _ in range(NTREE):
        R = np.where(Mtr, (Ytr - s) * w, 0.0)
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MAXFEAT,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(Xtr, R)
        s = s + LR * t.predict(Xtr) / w
        pred = pred + LR * t.predict(Xte) / w
    return pred


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--arms", default="независимо|пул|многозадачно|многозадачно, масштаб")
    ap.add_argument("--out", default=RES + "preds/oof_multi.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    arms = a.arms.split("|")

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)
    n = len(rows)

    Y = np.zeros((n, 4)); M = np.zeros((n, 4), bool)
    LO, HI = {}, {}
    for e, c in enumerate(CYPS):
        col = f"{c}_pIC50_direct_inhibition"
        M[:, e] = tr[col].notna().to_numpy()
        Y[M[:, e], e] = tr.loc[M[:, e], col].to_numpy()
        LO[c] = tr[col + "_conf_low"].to_numpy(float)
        HI[c] = tr[col + "_conf_high"].to_numpy(float)
    print(f"молекул {n}, помечено по ферментам {M.sum(0).tolist()}, "
          f"хотя бы одна метка у {int(M.any(1).sum())}, все четыре у {int(M.all(1).sum())}\n")

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for arm in arms:
            t0 = time.time()
            r = {"seed": seed, "рука": arm}
            P = {c: np.zeros(n) for c in CYPS}
            for f in range(5):
                trn_mol, te_mol = fold != f, fold == f
                if arm.startswith("многозадачно"):
                    keep = trn_mol & M.any(1)
                    pm = boostM(X[keep], Y[keep], M[keep], X[te_mol],
                                seed * 10 + f, "масштаб" in arm)
                    for e, c in enumerate(CYPS):
                        P[c][te_mol] = pm[:, e]
                elif arm == "пул":
                    Xs, ys = [], []
                    for e, c in enumerate(CYPS):
                        sel = M[:, e] & trn_mol
                        ind = np.zeros((int(sel.sum()), 4), np.float32); ind[:, e] = 1.0
                        Xs.append(np.hstack([X[sel], ind])); ys.append(Y[sel, e])
                    Xtr, ytr = np.vstack(Xs), np.concatenate(ys)
                    for e, c in enumerate(CYPS):
                        te = M[:, e] & te_mol
                        if not te.any():
                            continue
                        ind = np.zeros((int(te.sum()), 4), np.float32); ind[:, e] = 1.0
                        P[c][te] = boost1(Xtr, ytr, np.hstack([X[te], ind]), seed * 10 + f)
                else:
                    for e, c in enumerate(CYPS):
                        trn, te = M[:, e] & trn_mol, M[:, e] & te_mol
                        if not te.any():
                            continue
                        P[c][te] = boost1(X[trn], Y[trn, e], X[te], seed * 10 + f)

            for e, c in enumerate(CYPS):
                m = M[:, e]
                y, lo, hi, fi = Y[m, e], LO[c][m], HI[c][m], fold[m]
                p = P[c][m]
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{arm}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for t in ("пара", "rho"):
                r[f"MACRO {t}"] = round(float(np.mean([r[f"{c} {t}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {arm:24s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    g = df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                       + [f"{c} rho" for c in CYPS]].mean()
    print("\n" + g.round(4).to_string())
    if "независимо" in g.index:
        print(f"\n{'рука':24s} " + "".join(f"{c[3:]+' Δ':>12s}" for c in CYPS))
        for arm in g.index:
            if arm == "независимо":
                continue
            line = f"{arm:24s} "
            for c in CYPS:
                d = g.loc[arm, f"{c} rho"] - g.loc["независимо", f"{c} rho"]
                line += f"{d:+11.4f}{'*' if abs(d) > FLOOR[c] else ' '}"
            print(line)
        print("  * --- больше СОБСТВЕННОГО пола фермента (пункт 165)")
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Решает многозадачная рука против «независимо» ПРИ ГЛУБИНЕ 5. Пункт 180 намерил,
что пул на этой глубине проигрывает, потому что тратит её на изоляцию фермента; многозадачное
дерево получает то же обусловливание в листе и не тратит ничего. Выигрыш подтверждает диагноз.
Проигрыш вровень с пулом опровергает его: тогда дело не в глубине.

Рука «пул» здесь --- не эталон, а именно та, которая обязана проиграть, и её присутствие
проверяет, что стенд воспроизводит пункт 180, а не меряет что-то своё.""")


if __name__ == "__main__":
    main()
