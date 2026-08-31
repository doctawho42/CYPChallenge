"""Is the fingerprint load-bearing or corrective, and can it be replaced by a nameable few.

Two numbers that do not contradict each other and together say something neither says alone.

Item 96 measured that dropping the fingerprint costs 0.027 of pair and 0.032 of rank -- the
largest feature effect in this file. Item 138 measured where splits actually go: **DESC takes 58
to 66 per cent of them while being 9.5 per cent of the columns, and FP takes 30 to 38 per cent
while being 89 per cent.** Per column, a descriptor is chosen about seventeen times more often
than a bit.

So the fingerprint matters a great deal and is consulted rarely. The reading that fits both is
that it is not carrying the signal but **correcting** it: a small number of specific bits repair
specific chemotypes the descriptors describe wrongly, and the other two thousand do nothing.

That is testable, and the test is cheap because it needs no new model class. Rank the bits by how
often they are chosen across the boosting, keep the top k, discard the rest, and see how much of
item 96's 0.032 survives.

    k = 20, 50, 100, 200 against the full 2048

If most of it survives at k = 50, the fingerprint block reduces to fifty named substructures, and
that is a chemical statement rather than a numerical one -- a list of structural alerts with an
account of what the descriptors get wrong. Of everything discussed in this file, it is the first
thing that could produce a sentence about chemistry instead of a fourth decimal.

If nothing survives at any k, the bits work in bulk, there is no interpretable core, and the
question is closed the other way.

Bit importance is counted on the training folds only, then the reduced matrix is used for a fresh
out-of-fold run -- selecting bits on the whole set and scoring on it would be selection on the
answer.

Reads data/feats.npz and rows.csv. Writes results/preds/oof_topk.json.
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
from sklearn.tree import DecisionTreeRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
KS = (0, 20, 50, 100, 200, 2048)
NTREE, LR, DEPTH, MF = 150, 0.06, 5, 0.3


def bit_counts(X, y, nfp, seed):
    """Сколько раз каждый бит выбран сплитом. Считается своим бустингом, потому что у
    HistGradientBoostingRegressor нет доступа к признакам сплитов, а важности он не даёт."""
    rng = np.random.default_rng(seed)
    s = np.full(len(y), float(y.mean()))
    cnt = np.zeros(nfp, int)
    for _ in range(NTREE):
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MF,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(X, y - s)
        s = s + LR * t.predict(X)
        f = t.tree_.feature
        f = f[(f >= 0) & (f < nfp)]
        np.add.at(cnt, f, 1)
    return cnt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--out", default=RES + "preds/oof_topk.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    FP, DESC, MECH = z["FP"], z["DESC"], z["MECH"]
    nfp = FP.shape[1]
    DM = np.hstack([DESC, MECH]).astype(np.float32)

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for k in KS:
            t0 = time.time()
            r = {"seed": seed, "k": k}
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                fi = fold[m]
                Fi, Di = FP[m].astype(np.float32), DM[m]
                p = np.zeros_like(y)
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    if k == 0:
                        A, B = Di[trn], Di[te]
                    elif k >= nfp:
                        A = np.hstack([Fi[trn], Di[trn]])
                        B = np.hstack([Fi[te], Di[te]])
                    else:
                        # Биты ранжируются ТОЛЬКО на обучающих фолдах.
                        cnt = bit_counts(np.hstack([Fi[trn], Di[trn]]), y[trn], nfp,
                                         seed * 10 + f)
                        keep = np.argsort(-cnt)[:k]
                        A = np.hstack([Fi[trn][:, keep], Di[trn]])
                        B = np.hstack([Fi[te][:, keep], Di[te]])
                    p[te] = HistGradientBoostingRegressor(**KW).fit(A, y[trn]).predict(B)
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{k}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for tag in ("пара", "rho"):
                r[f"MACRO {tag}"] = round(float(np.mean([r[f"{c} {tag}"] for c in CYPS])), 4)
            table.append(r)
            lab = "без FP" if k == 0 else ("все 2048" if k >= nfp else f"top-{k}")
            print(f"  сид {seed} {lab:10s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("k")[["MACRO пара", "MACRO rho"]].mean().to_string())
    full = df[df.k >= nfp]["MACRO rho"].mean()
    none = df[df.k == 0]["MACRO rho"].mean()
    print(f"\nвесь эффект фингерпринта по рангу: {full - none:+.4f}")
    for k in KS:
        if 0 < k < nfp:
            v = df[df.k == k]["MACRO rho"].mean()
            print(f"  при top-{k:<4d} выживает {100*(v-none)/max(full-none,1e-9):5.1f} %")
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Строки «без FP» и «все 2048» задают шкалу: их разность --- это весь эффект
фингерпринта, и пункт 96 намерил его в 0.032 ранга.

Если при top-50 выживает большая часть, фингерпринтовый блок сводится к полусотне
подструктур, и следующий шаг --- расшифровать их и посмотреть, чего именно не хватает
дескрипторам. Это уже химическое высказывание.

Если не выживает ни при каком k, биты работают массой, интерпретируемого ядра нет, и вопрос
закрыт с другой стороны --- тоже результат.""")


if __name__ == "__main__":
    main()
