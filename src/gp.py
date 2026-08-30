"""Exact Gaussian process on descriptors: a different family, not another variant of boosting.

Why this and not the usual Tanimoto kernel. The proposal was a GP with a Tanimoto kernel over
Morgan fingerprints, which is the standard in cheminformatics. Item 90 measured which space
actually carries activity information between neighbours and Morgan came last of four, carrying
nothing at all on CYP2C9, while standardised RDKit descriptors came first by a factor of 2.5.
A similarity method is only as good as its similarity, so the kernel goes on descriptors.

Why a GP is worth a run when five other ideas were not. Its posterior mean is a
similarity-weighted average of neighbouring labels, which is exactly the local structure a tree
cannot represent: trees cut on axes and chemical similarity is not axis-aligned. And the setting
favours it - item 22 measured the test set sitting closer to the training set (0.587) than the
training set sits to itself (0.435), so almost every test molecule has a near relative to borrow
from. By item 80's criterion this is the kind of intervention that can move rank rather than
scale, which is the only kind that survives the affine pair.

Exact, not approximate. With 1285 to 2335 rows per enzyme the Cholesky is seconds, so there is
no reason for inducing points or their approximation error.

Hyperparameters by marginal likelihood on the training folds only - a lengthscale multiplier
around the median-distance heuristic and a noise level - never by held-out score. The
descriptors are standardised and clipped at five deviations, because Ipc alone spans fourteen
orders of magnitude (item 64) and would otherwise define the metric by itself.

The predictive variance is saved as well. The heteroscedastic decision layer needs a per-compound
spread and has so far been fed one derived from quantile regression; a GP produces it natively.

Writes results/preds/oof_gp.json with predictions and variances.
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
from scipy.linalg import cho_factor, cho_solve
from scipy.stats import spearmanr
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
# Сетка намеренно шире, чем нужно любому правдоподобному значению: при первом прогоне
# оптимум сел на нижний край длины и верхний край шума на трёх ферментах из четырёх, то есть
# выбирала сетка, а не маргинальное правдоподобие. Та же ошибка, что в пункте 59.
LS_MULT = [0.0625, 0.125, 0.25, 0.5, 1.0, 2.0, 4.0]
NOISE = [0.01, 0.03, 0.1, 0.3, 1.0, 3.0]


def sqdist(A, B):
    return np.maximum(0.0, (A * A).sum(1)[:, None] + (B * B).sum(1)[None, :] - 2.0 * A @ B.T)


def fit_predict(Xt, yt, Xe, ls0, var_y):
    """Подбор (длина, шум) по маргинальному правдоподобию на обучающих строках, затем прогноз."""
    mu = yt.mean()
    yc = yt - mu
    D2 = sqdist(Xt, Xt)
    best = None
    for lm in LS_MULT:
        K0 = np.exp(-D2 / (2.0 * (ls0 * lm) ** 2))
        for nz in NOISE:
            K = var_y * K0 + np.eye(len(yt)) * var_y * nz
            try:
                c = cho_factor(K, lower=True)
            except np.linalg.LinAlgError:
                continue
            al = cho_solve(c, yc)
            lml = -0.5 * yc @ al - np.log(np.diag(c[0])).sum()
            if best is None or lml > best[0]:
                best = (lml, lm, nz, c, al)
    _, lm, nz, c, al = best
    Ks = var_y * np.exp(-sqdist(Xe, Xt) / (2.0 * (ls0 * lm) ** 2))
    m = mu + Ks @ al
    v = var_y * (1.0 + nz) - np.einsum("ij,ij->i", Ks, cho_solve(c, Ks.T).T)
    return m, np.maximum(v, 1e-9), lm, nz


def prepare(X):
    """Стандартизация с обрезкой: Ipc иначе задаёт метрику в одиночку (пункт 64)."""
    X = np.nan_to_num(np.asarray(X, np.float64), posinf=0.0, neginf=0.0)
    mu, sd = X.mean(0), X.std(0) + 1e-9
    return lambda A: np.clip((np.nan_to_num(np.asarray(A, np.float64), posinf=0.0, neginf=0.0)
                              - mu) / sd, -5.0, 5.0)


def gp_predict(Xtr, ytr, Xte, seed=0):
    """Апостериорное среднее GP, гиперпараметры по маргинальному правдоподобию на Xtr."""
    sub = Xtr[np.random.default_rng(seed).choice(len(Xtr), min(600, len(Xtr)), replace=False)]
    ls0 = float(np.sqrt(np.median(sqdist(sub, sub)[np.triu_indices(len(sub), 1)])))
    m, _, _, _ = fit_predict(Xtr, ytr, Xte, ls0, ytr.var())
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--blocks", default="DESC+MECH")
    ap.add_argument("--out", default=RES + "preds/oof_gp.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z[b] for b in a.blocks.split("+")]).astype(np.float64)
    X = np.nan_to_num(X, posinf=0.0, neginf=0.0)
    X = (X - X.mean(0)) / (X.std(0) + 1e-9)
    X = np.clip(X, -5.0, 5.0)          # Ipc иначе задаёт метрику в одиночку
    print(f"ядро на {a.blocks}, {X.shape[1]} колонок после стандартизации\n")

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        r = {"seed": seed}
        for c in CYPS:
            t0 = time.time()
            col = f"{c}_pIC50_direct_inhibition"
            m = tr[col].notna().to_numpy()
            y = tr.loc[m, col].to_numpy()
            lo = tr.loc[m, col + "_conf_low"].to_numpy()
            hi = tr.loc[m, col + "_conf_high"].to_numpy()
            Xi, fi = X[m], fold[m]
            # Медианная эвристика для длины, на подвыборке чтобы не строить всю матрицу.
            sub = Xi[np.random.default_rng(0).choice(len(Xi), min(600, len(Xi)), replace=False)]
            ls0 = float(np.sqrt(np.median(sqdist(sub, sub)[np.triu_indices(len(sub), 1)])))
            p = np.zeros_like(y)
            v = np.zeros_like(y)
            used = []
            for f in range(5):
                trn, te = fi != f, fi == f
                if te.sum() == 0:
                    continue
                p[te], v[te], lm, nz = fit_predict(Xi[trn], y[trn], Xi[te], ls0, y[trn].var())
                used.append((lm, nz))
            q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
            out[f"{seed}|GP|{c}"] = p.tolist()
            out[f"{seed}|GPvar|{c}"] = v.tolist()
            r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
            r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            edge = any(l in (LS_MULT[0], LS_MULT[-1]) or n in (NOISE[0], NOISE[-1])
                       for l, n in used)
            print(f"  сид {seed} {c}: пара {r[f'{c} пара']:.4f} rho {r[f'{c} rho']:.4f}  "
                  f"длина x{used[0][0]} шум {used[0][1]}"
                  + ("  ВНИМАНИЕ: край сетки" if edge else "")
                  + f"  ({time.time()-t0:.0f} с)", flush=True)
        r["MACRO пара"] = round(float(np.mean([r[f"{c} пара"] for c in CYPS])), 4)
        r["MACRO rho"] = round(float(np.mean([r[f"{c} rho"] for c in CYPS])), 4)
        table.append(r)
        print(f"  сид {seed}: макро пара {r['MACRO пара']:.4f}  rho {r['MACRO rho']:.4f}\n", flush=True)

    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"сохранено: {a.out}")
    print("""
Как читать. Сравнивать с бустингом надо не столько по метрике, сколько по РАНГУ, и не столько
в одиночку, сколько в ансамбле: смысл другого семейства в том, что оно ошибается иначе. У
бустинга ранг 0.5630 раздельно и 0.5767 пулом, у их ансамбля 0.5921.

Если ранг GP заметно ниже, но ансамбль с бустингом всё равно выше 0.5921 --- метод полезен
именно тем, чем задумывался. Если и ранг ниже, и ансамбль не двигается --- дерево уже забрало
всё, что несут дескрипторы, и локальной структуры сверх этого нет.""")


if __name__ == "__main__":
    main()
