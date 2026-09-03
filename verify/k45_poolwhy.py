"""Why pooling gains on HistGB and reverses on plain trees: the routing tax hypothesis.

The problem. Pooling is the effect the submission is built around. It gains +0.0141 of rank on
HistGradientBoostingRegressor (0.5792 against the per-enzyme 0.5651) and **loses 0.027 on plain
depth-5 trees** (0.5416 against 0.5684). Item 158 refuted the first explanation -- column
subsampling starving the enzyme indicator -- by running at max_features=1.0 and getting 0.5416
against 0.5420, identical. A sign reversal between two tree ensembles on the same rows, the same
folds and the same indicator has to have a mechanism, and until it has one the largest effect in
the submission rests on a learner rather than on the data.

The hypothesis, which is about how the two learners spend their budget rather than about capacity.

    HistGB      max_leaf_nodes=31, max_depth=None  -> лучший-первым, глубина не ограничена
    свой        max_depth=5, max_leaf_nodes=None   -> вглубь, КАЖДЫЙ путь не длиннее пяти

A pooled table forces every tree to route before it can model. With a four-column one-hot,
separating the four enzymes costs up to **three** splits on the worst branch: ind0 peels off
CYP1A2, ind1 peels off CYP2C9, ind2 separates CYP2D6 from CYP3A4. Under a depth cap of five that
leaves **two** levels for chemistry on that branch. Under best-first growth with a leaf budget
there is no cap at all -- the tree spends splits where the gain is, and a branch that has just paid
three splits for routing can keep growing.

So the prediction is not "pooling is good" or "pooling is bad" but that **the sign of pooling's
effect depends on the growth policy, and on nothing else here**.

The grid, six arms, one learner difference at a time:

    независимо, глубина 5     the reference, 0.5684
    пул, глубина 5            the reversal, 0.5416
    пул, глубина 8            more depth, same policy -- separates budget from policy
    пул, листья 31            best-first, no depth cap -- HistGB's policy on plain trees
    пул, глубина 5, 400 дер.  more trees, same depth -- separates capacity-by-count
    пул, глубина 5, 800 дер.  **the capacity-matched arm.** The per-enzyme setting fits four
                              models of 200 trees, so it receives 800 trees in total against the
                              pooled setting's 200, on a table five times smaller per model. That
                              asymmetry has been in every pooling comparison in this file,
                              including HistGB's (1200 against 300), and has never been controlled.
                              If pooling recovers here the reversal is capacity accounting, and
                              HistGB's +0.0141 becomes a stronger result rather than a weaker one:
                              it wins *despite* a fourfold handicap
    независимо, листья 31     the control: the per-enzyme model has no indicator to route
                              around, so the growth policy should barely matter to it

Pre-registered readings:

  **пул, листья 31 recovers the gain** -> the routing tax is the mechanism, the reversal is a
  property of depth-limited growth, and pooling's HistGB number stands as a statement about the
  data rather than about the learner.

  **пул, глубина 8 recovers it and листья 31 adds nothing more** -> it is the depth budget in
  particular, quantitatively, and the same tax applies to HistGB whenever its trees happen to be
  shallow.

  **400 деревьев recovers it** -> it was never routing, it was capacity, and the pooled table
  simply needs more of it than 200 depth-5 trees supply.

  **none recovers it** -> the growth policy is not the mechanism either, and what is left is the
  l2 regularisation on leaf values and the binning, neither of which has an obvious route to a
  sign flip. That outcome would make the pooling result learner-specific in a way that should be
  said out loud in the document.

The tax is also measured directly rather than only inferred from the outcome: for the fitted pooled
trees, at what depth the indicator columns are used, what share of all splits they take, and how
many indicator splits sit on the average root-to-leaf path. If the hypothesis is right, the
depth-5 pooled trees should spend a large share of their shallow levels on routing.

Same folds, masks and metric as everywhere. max_features is held at 0.3 throughout, since item 158
measured it irrelevant to the reversal. Writes results/preds/oof_poolwhy.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time
from collections import Counter

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.tree import DecisionTreeRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
LR, MAXFEAT = 0.06, 0.3
NIND = 5                      # четыре фермента плюс колонка источника


def depths(t):
    """Глубина каждого узла дерева."""
    d = np.zeros(t.tree_.node_count, int)
    st = [(0, 0)]
    while st:
        n, dep = st.pop()
        d[n] = dep
        l, r = t.tree_.children_left[n], t.tree_.children_right[n]
        if l != -1:
            st.append((l, dep + 1)); st.append((r, dep + 1))
    return d


def tax(t, nchem):
    """Пошлина за маршрутизацию: сплиты по индикатору --- сколько их, на какой глубине,
    и сколько их на среднем пути от корня до листа."""
    f, dd = t.tree_.feature, depths(t)
    ind = (f >= nchem) & (f >= 0)
    chem = (f >= 0) & (f < nchem)
    # число индикаторных сплитов на пути к каждому листу
    par = np.full(t.tree_.node_count, -1)
    for n in range(t.tree_.node_count):
        for c in (t.tree_.children_left[n], t.tree_.children_right[n]):
            if c != -1:
                par[c] = n
    onpath = []
    for n in range(t.tree_.node_count):
        if t.tree_.children_left[n] != -1:
            continue
        k, p = 0, par[n]
        while p != -1:
            if ind[p]:
                k += 1
            p = par[p]
        onpath.append(k)
    return (int(ind.sum()), int(chem.sum()),
            list(dd[ind]), float(np.mean(onpath)) if onpath else 0.0,
            float(np.mean(dd[[i for i in range(len(f)) if t.tree_.children_left[i] == -1]])))


def boost(Xtr, ytr, Xte, seed, ntree, kw, nchem=None, diag=None):
    rng = np.random.default_rng(seed)
    base = float(ytr.mean())
    s = np.full(len(ytr), base)
    pred = np.full(len(Xte), base)
    for _ in range(ntree):
        t = DecisionTreeRegressor(max_features=MAXFEAT,
                                  random_state=int(rng.integers(1 << 30)), **kw)
        t.fit(Xtr, ytr - s)
        if diag is not None and nchem is not None:
            diag.append(tax(t, nchem))
        s = s + LR * t.predict(Xtr)
        pred = pred + LR * t.predict(Xte)
    return pred


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--out", default=RES + "preds/oof_poolwhy.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)
    nchem = X.shape[1]

    Y, LO, HI, M = {}, {}, {}, {}
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        M[c] = tr[col].notna().to_numpy()
        Y[c] = tr[col].to_numpy(float)
        LO[c] = tr[col + "_conf_low"].to_numpy(float)
        HI[c] = tr[col + "_conf_high"].to_numpy(float)

    ARMS = [
        ("независимо, глубина 5",      False, 200, dict(max_depth=5)),
        ("пул, глубина 5",             True,  200, dict(max_depth=5)),
        ("пул, глубина 8",             True,  200, dict(max_depth=8)),
        ("пул, листья 31",             True,  200, dict(max_leaf_nodes=31)),
        ("пул, глубина 5, 400 дер.",   True,  400, dict(max_depth=5)),
        ("пул, глубина 5, 800 дер.",   True,  800, dict(max_depth=5)),
        ("независимо, листья 31",      False, 200, dict(max_leaf_nodes=31)),
    ]

    out, table, diagrec = {}, [], {}
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for arm, pooled, ntree, kw in ARMS:
            t0 = time.time()
            r = {"seed": seed, "рука": arm}
            p_all = {c: np.zeros(len(rows)) for c in CYPS}
            diag = [] if pooled else None
            for f in range(5):
                trn_mol, te_mol = fold != f, fold == f
                if pooled:
                    Xs, ys = [], []
                    for e, c in enumerate(CYPS):
                        sel = M[c] & trn_mol
                        ind = np.zeros((int(sel.sum()), NIND), np.float32)
                        ind[:, e] = 1.0
                        Xs.append(np.hstack([X[sel], ind]))
                        ys.append(Y[c][sel])
                    Xtr, ytr = np.vstack(Xs), np.concatenate(ys)
                    for e, c in enumerate(CYPS):
                        te = M[c] & te_mol
                        if not te.any():
                            continue
                        ind = np.zeros((int(te.sum()), NIND), np.float32)
                        ind[:, e] = 1.0
                        p_all[c][te] = boost(Xtr, ytr, np.hstack([X[te], ind]),
                                             seed * 10 + f, ntree, kw, nchem,
                                             diag if f == 0 else None)
                else:
                    for c in CYPS:
                        trn, te = M[c] & trn_mol, M[c] & te_mol
                        if not te.any():
                            continue
                        p_all[c][te] = boost(X[trn], Y[c][trn], X[te],
                                             seed * 10 + f, ntree, kw)
            for c in CYPS:
                m = M[c]
                y, lo, hi, fi = Y[c][m], LO[c][m], HI[c][m], fold[m]
                p = p_all[c][m]
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{arm}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for tag in ("пара", "rho"):
                r[f"MACRO {tag}"] = round(float(np.mean([r[f"{c} {tag}"] for c in CYPS])), 4)
            table.append(r)
            if diag:
                ni = sum(d[0] for d in diag); nc = sum(d[1] for d in diag)
                dep = [x for d in diag for x in d[2]]
                diagrec[arm] = dict(доля_индикаторных=ni / max(ni + nc, 1),
                                    медиана_глубины=float(np.median(dep)) if dep else np.nan,
                                    индик_на_пути=float(np.mean([d[3] for d in diag])),
                                    средняя_глубина_листа=float(np.mean([d[4] for d in diag])))
            print(f"  сид {seed} {arm:26s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    g = df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]].mean()
    print("\n" + g.round(4).to_string())
    base = g.loc["независимо, глубина 5", "MACRO rho"]
    print(f"\n{'рука':28s} {'ранг':>8s} {'против поферментной':>20s}")
    for k in g.index:
        print(f"{k:28s} {g.loc[k, 'MACRO rho']:8.4f} {g.loc[k, 'MACRO rho'] - base:+20.4f}")
    print(f"\nдля справки, HistGB на тех же фолдах: поферментно 0.5651, пул 0.5792 (+0.0141)")

    print("\nПОШЛИНА ЗА МАРШРУТИЗАЦИЮ, замерена на деревьях фолда 0")
    print(f"{'рука':26s} {'доля сплитов':>13s} {'медиана глубины':>16s} "
          f"{'индик. на пути':>15s} {'глубина листа':>14s}")
    for k, v in diagrec.items():
        print(f"{k:26s} {100*v['доля_индикаторных']:12.1f} % {v['медиана_глубины']:16.1f} "
              f"{v['индик_на_пути']:15.2f} {v['средняя_глубина_листа']:14.2f}")
    json.dump({"table": table, "diag": diagrec, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Решает столбец «против поферментной». Пулирование на HistGB даёт +0.0141; если
какая-то из рук здесь выходит в плюс, найден тот параметр учителя, от которого зависит ЗНАК.

«листья 31» против «глубина 8» разделяет политику роста и бюджет глубины: первая снимает
ограничение на глубину вовсе, вторая только раздвигает его.

«400 деревьев» проверяет, не было ли всё это просто нехваткой ёмкости на таблице вшестеро
большей.

Таблица пошлины --- прямое измерение механизма, а не вывод из исхода. «индик. на пути» ---
сколько сплитов по индикатору фермента стоит на среднем пути от корня до листа. Если при
глубине 5 это около единицы-двух из пяти, то пятая часть или больше бюджета дерева уходит
не на химию, и гипотеза подтверждена независимо от того, что показал ранг.""")


if __name__ == "__main__":
    main()
