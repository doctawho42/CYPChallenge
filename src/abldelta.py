"""The Delta >= 0 constraint the document specifies and the code never implemented.

Section 4 defines the two arms of the assay as one latent with a monotone offset,

    pi_{m,e,tdi} = pi_{m,e,dir} + Delta_{m,e},   Delta >= 0      [s04.tex:23]

because pre-incubation can only strengthen inhibition, never weaken it. The code has no such
constraint anywhere -- every softplus and ReLU in `src/trunk.py` sits on the Hill slope or on a
network layer -- and item 125 measured the unconstrained version, pooling the two arms into one
table, at nothing.

**The precondition holds and was checked before building.** The share of compounds whose measured
Delta is negative by more than twice its own propagated error:

    фермент   обе метки   медиана Delta   доля D<0   D < -2 sigma
    CYP1A2         1412           0.024      38.5 %          0.0 %
    CYP2C9         1285           0.007      46.8 %          1.9 %
    CYP2D6         1493           0.154      13.5 %          1.1 %
    CYP3A4         2334           0.307      13.3 %          0.7 %

Under two per cent everywhere, so the 13 to 47 per cent of negative values are measurement noise
around a non-negative truth and the inequality is legitimate to impose.

**And the medians pre-register the result against the naive expectation.** On CYP1A2 and CYP2C9 the
median offset is 0.024 and 0.007 -- the pre-incubation arm is very nearly a **repeat measurement of
the same quantity**. Constraining Delta to be non-negative and letting it shrink toward zero turns
1412 and 1285 TDI labels into a second reading of direct inhibition, which is variance reduction on
the target itself. On CYP3A4 the median is 0.307 and the arm measures something else, so little
should follow. **The gain should therefore be largest where Delta is smallest**, which is the
opposite of where a "more data about inactivation" reading would put it.

Arms. The free version has to share structure, or it reduces to the reference: with two independent
models the direct head never sees a TDI label and nothing can change.

    только прямое           эталон, поферментно
    два выхода свободно     общие сплиты, выходы (pi_dir, pi_tdi) без связи между ними
    Delta >= 0              общие сплиты, выходы (pi_dir, raw), pi_tdi = pi_dir + softplus(raw)
    Delta >= 0, перемешан   TDI-метки переставлены между молекулами. Если ограничение помогает
                            и с ними, работает регуляризация, а не физика

Scored on direct inhibition only, per enzyme, against item 165's floors. Writes
results/preds/oof_delta.json.
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
LAM = 1.0          # вес TDI-плеча в потере


def softplus(x):
    return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0.0)


def sigmoid(x):
    return 0.5 * (1.0 + np.tanh(0.5 * x))


def boost1(Xtr, ytr, Xte, seed):
    rng = np.random.default_rng(seed)
    s = np.full(len(ytr), float(ytr.mean())); p = np.full(len(Xte), s[0])
    for _ in range(NTREE):
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MAXFEAT,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(Xtr, ytr - s)
        s = s + LR * t.predict(Xtr); p = p + LR * t.predict(Xte)
    return p


def boost2(Xtr, yd, md, yt, mt, Xte, seed, constrained):
    """Общие сплиты, два выхода. constrained: pi_tdi = pi_dir + softplus(raw)."""
    rng = np.random.default_rng(seed)
    a0 = float(yd[md].mean()) if md.any() else 0.0
    b0 = (float(yt[mt].mean()) if mt.any() else a0) if not constrained else 0.0
    S = np.column_stack([np.full(len(yd), a0), np.full(len(yd), b0)])
    P = np.column_stack([np.full(len(Xte), a0), np.full(len(Xte), b0)])
    for _ in range(NTREE):
        if constrained:
            dlt = softplus(S[:, 1])
            r_t = np.where(mt, yt - (S[:, 0] + dlt), 0.0)
            g0 = np.where(md, yd - S[:, 0], 0.0) + LAM * r_t
            g1 = LAM * r_t * sigmoid(S[:, 1])
        else:
            g0 = np.where(md, yd - S[:, 0], 0.0)
            g1 = LAM * np.where(mt, yt - S[:, 1], 0.0)
        G = np.column_stack([g0, g1])
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MAXFEAT,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(Xtr, G)
        S = S + LR * t.predict(Xtr)
        P = P + LR * t.predict(Xte)
    return P[:, 0]              # оценивается только прямое плечо


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/oof_delta.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    td = (pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv")
            .set_index("Molecule_Name").reindex(rows.Molecule_Name).reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)

    print(f"{'фермент':8s} {'прямых':>8s} {'TDI':>6s} {'обе':>6s} {'медиана Delta':>14s}")
    for c in CYPS:
        a_, b_ = f"{c}_pIC50_direct_inhibition", f"{c}_pIC50_TDI_condition"
        md, mt = tr[a_].notna().to_numpy(), td[b_].notna().to_numpy()
        both = md & mt
        print(f"{c:8s} {md.sum():8d} {mt.sum():6d} {both.sum():6d} "
              f"{np.median(td.loc[both, b_].to_numpy()-tr.loc[both, a_].to_numpy()):14.3f}")
    print()

    ARMS = ["только прямое", "два выхода свободно", "Delta >= 0", "Delta >= 0, перемешан"]
    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for arm in ARMS:
            t0 = time.time()
            r = {"seed": seed, "рука": arm}
            for c in CYPS:
                a_, b_ = f"{c}_pIC50_direct_inhibition", f"{c}_pIC50_TDI_condition"
                md = tr[a_].notna().to_numpy()
                yd = tr[a_].to_numpy(float)
                mt = td[b_].notna().to_numpy()
                yt = td[b_].to_numpy(float)
                if "перемешан" in arm:
                    g = np.random.default_rng(seed * 41 + 7)
                    idx = np.where(mt)[0]
                    yt = yt.copy(); yt[idx] = yt[idx][g.permutation(len(idx))]
                lo = tr[a_ + "_conf_low"].to_numpy(float)
                hi = tr[a_ + "_conf_high"].to_numpy(float)
                P = np.zeros(len(rows))
                for f in range(5):
                    trn_mol, te_mol = fold != f, fold == f
                    te = md & te_mol
                    if not te.any():
                        continue
                    if arm == "только прямое":
                        A = md & trn_mol
                        P[te] = boost1(X[A], yd[A], X[te], seed * 10 + f)
                    else:
                        keep = trn_mol & (md | mt)
                        P[te] = boost2(X[keep], np.nan_to_num(yd)[keep], md[keep],
                                       np.nan_to_num(yt)[keep], mt[keep],
                                       X[te], seed * 10 + f, "Delta" in arm)
                y, fi, p = yd[md], fold[md], P[md]
                q = fit_apply(p, lo[md], hi[md], fi, np.ones(int(md.sum())) / int(md.sum()))
                out[f"{seed}|{arm}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi[md],
                                                   y_true_lower=lo[md])), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for t in ("пара", "rho"):
                r[f"MACRO {t}"] = round(float(np.mean([r[f"{c} {t}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {arm:23s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    g = df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                       + [f"{c} rho" for c in CYPS]].mean()
    print("\n" + g.round(4).to_string())
    if "только прямое" in g.index:
        print(f"\n{'рука':23s} " + "".join(f"{c[3:]+' Δ':>12s}" for c in CYPS))
        for arm in g.index:
            if arm == "только прямое":
                continue
            line = f"{arm:23s} "
            for c in CYPS:
                d = g.loc[arm, f"{c} rho"] - g.loc["только прямое", f"{c} rho"]
                line += f"{d:+11.4f}{'*' if abs(d) > FLOOR[c] else ' '}"
            print(line)
        print("  * --- больше СОБСТВЕННОГО пола фермента (пункт 165)")
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Предрегистрация: выигрыш обязан быть НАИБОЛЬШИМ там, где медиана Delta наименьшая
--- на CYP1A2 (0.024) и CYP2C9 (0.007), где преинкубационное плечо есть почти повторный замер
той же величины, и ограничение превращает его во второй отсчёт прямого ингибирования. На
CYP3A4 (0.307) плечо меряет другое, и выигрыша ждать не надо. Обратный порядок опровергает
это прочтение.

«Два выхода свободно» отделяет ограничение от простого разделения структуры: там те же
общие сплиты и те же TDI-метки, но связи между плечами нет --- это и есть версия, которую
пункт 125 намерил в ноль.

Перемешанная рука обязана НЕ помогать, иначе работает регуляризация, а не физика.""")


if __name__ == "__main__":
    main()
