"""The enzyme as a coordinate instead of a one-hot, with the coordinate measured on our own data.

The argument. Item 132 established that pooling works by contrast and that removing the enzyme
indicator drops the model below training each enzyme separately -- the indicator is load-bearing.
It is currently four arbitrary orthogonal units, all equidistant. The enzymes are not equidistant.
Measured on the screening four-vector, which is complete on 4376 molecules and carries no
selection:

              1A2    2C9    2D6    3A4
    1A2     1.000  0.478  0.161  0.398
    2C9     0.478  1.000  0.243  0.719
    2D6     0.161  0.243  1.000  0.331
    3A4     0.398  0.719  0.331  1.000

CYP2C9 and CYP3A4 are nearly the same enzyme at 0.719; CYP2D6 stands apart at 0.161 from CYP1A2.
A one-hot tells the model none of that, so every contrast it learns has to be learned separately
for each pair from scratch.

Where the coordinate comes from, and why not from the protein. The textbook move is
proteochemometrics: encode each isoform by physicochemical z-scales of its active-site residues.
That imports a geometry from a sequence model and needs the residue selection to be right, which
is a second thing to be wrong about. The geometry here is instead **taken from the data**: classical
multidimensional scaling of the screening correlation matrix, two coordinates. It says only what
our own measurements say, and it needs no external resource.

That makes the falsification clean. If a measured-geometry coordinate beats a one-hot, the
equidistance of the one-hot was costing something and the proteochemometric version is worth
building on top. If it does not, the geometry is not what the indicator carries, and the
sequence-based version -- three to four days, mostly assembling external isoforms -- should not be
started.

Arms:

    one-hot            what src/ablpool.py does now. Must reproduce 0.7062 / 0.5792 at seed 0.
    координата         two MDS coordinates in place of the four indicator columns.
    one-hot + координата
                       both, in case the coordinate adds without the indicator being redundant.

The width is kept honest: the one-hot arm carries four columns and the coordinate arm two, so the
coordinate arm is *narrower*. A win there cannot be dimensionality.

Same folds, masks, metric and learner settings as src/ablpool.py. Writes
results/preds/oof_coord.json.
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
NDIM = 2


def enzyme_coords(D_):
    """Классическое многомерное шкалирование корреляционной матрицы скрининга.

    Расстояние берётся как sqrt(2 (1 - r)) --- обычная метрика для корреляций, --- и
    двойное центрирование даёт координаты, в которых близкие ферменты стоят рядом.
    """
    n = len(D_)
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D_ ** 2) @ J
    w, V = np.linalg.eigh(B)
    idx = np.argsort(-w)[:NDIM]
    return V[:, idx] * np.sqrt(np.maximum(w[idx], 0.0))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--out", default=RES + "preds/oof_coord.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    piv = sc.pivot_table(index="Molecule_Name", columns="enzyme", values="log2fc_estimate")
    A = piv[CYPS].dropna().to_numpy(float)
    C = np.corrcoef(A.T)
    coord = enzyme_coords(np.sqrt(np.maximum(2.0 * (1.0 - C), 0.0)))
    print("корреляции по скринингу (%d молекул):" % len(A))
    print("        " + " ".join(f"{c[3:]:>6s}" for c in CYPS))
    for i, c in enumerate(CYPS):
        print(f"  {c[3:]:5s} " + " ".join(f"{C[i, j]:6.3f}" for j in range(4)))
    print("\nкоординаты ферментов:")
    for i, c in enumerate(CYPS):
        print(f"  {c}: " + " ".join(f"{v:+.3f}" for v in coord[i]))
    print()

    M = {c: tr[f"{c}_pIC50_direct_inhibition"].notna().to_numpy() for c in CYPS}
    Y = {c: tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS}

    def tag(e, n, kind):
        if kind == "one-hot":
            v = np.zeros((n, 4), np.float32)
            v[:, e] = 1.0
        elif kind == "координата":
            v = np.tile(coord[e].astype(np.float32), (n, 1))
        else:
            oh = np.zeros((n, 4), np.float32)
            oh[:, e] = 1.0
            v = np.hstack([oh, np.tile(coord[e].astype(np.float32), (n, 1))])
        return v

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for kind in ("one-hot", "координата", "one-hot + координата"):
            t0 = time.time()
            r = {"seed": seed, "рука": kind}
            for e, c in enumerate(CYPS):
                m = M[c]
                y, fi = Y[c][m], fold[m]
                lo = tr.loc[m, f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy()
                hi = tr.loc[m, f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy()
                p = np.zeros_like(y)
                for f in range(5):
                    keep = fold != f
                    te = fi == f
                    if te.sum() == 0:
                        continue
                    Xs, ys = [], []
                    for e2, c2 in enumerate(CYPS):
                        sel = M[c2] & keep
                        if not sel.any():
                            continue
                        Xs.append(np.hstack([X[sel], tag(e2, int(sel.sum()), kind)]))
                        ys.append(Y[c2][sel])
                    sel_te = M[c] & (fold == f)
                    Xe = np.hstack([X[sel_te], tag(e, int(sel_te.sum()), kind)])
                    p[te] = HistGradientBoostingRegressor(**KW).fit(
                        np.vstack(Xs), np.concatenate(ys)).predict(Xe)
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{kind}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for t in ("пара", "rho"):
                r[f"MACRO {t}"] = round(float(np.mean([r[f"{c} {t}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {kind:22s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                         + [f"{c} rho" for c in CYPS]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Первая рука обязана дать 0.7062 / 0.5792 на сиде 0 --- это пулированная строка из
ablpool. Иначе сравнивать не с чем.

Координата УЖЕ one-hot: две колонки против четырёх. Поэтому выигрыш второй руки не может быть
размерностью, и это единственный способ отличить геометрию от ширины без отдельного контроля.

Решает вторая против первой. Выигрыш --- равноудалённость one-hot чего-то стоила, и
последовательностный вариант (три-четыре дня, в основном сборка внешних изоформ) осмыслен.
Ноль --- геометрия не то, что несёт индикатор, и тот вариант начинать не надо.""")


if __name__ == "__main__":
    main()
