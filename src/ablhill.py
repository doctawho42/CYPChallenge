"""The Hill residual as a training weight: the first per-compound signal that is not a function of the label.

Where this comes from. Item 170 established that the Hill slope is identifiable per compound --- its
spread exceeds the full measurement error (both sources) by 2.9 to 5.8 --- and that a single pair
(shift, slope) per enzyme does **not** absorb it: the residual after that fit stays 2.7 to 4.4 times
the measurement error. So there is a real per-compound quantity in the disagreement between the
curve and the single screening point.

What it is remains undetermined, and deliberately does not matter here. One equation and one point
on the curve identify exactly one parameter; which one is a choice. Solve for the slope with the
amplitude fixed and it reads as kinetics; solve for the amplitude with the slope fixed at one and
it reads as an Emax that varies per compound. Item 170 shows the data cannot separate them --- under
the repository's own fitted E the median slope lands at 1.00 on three enzymes, which is exactly the
physically expected value, so the "impossible slope" argument for the kinetic reading fails.

**This ablation needs neither reading.** It uses only the size of the disagreement, not its
attribution: a cell where the dose-response curve and the 49.5 uM screening point contradict each
other is a cell whose label is less trustworthy, whichever measurement is at fault.

Why that is worth a run under item 166's rule. The band is 3.92 times `_std`, and `_std` is a
function of the label with R-squared 0.93 to 0.98, so **every per-compound quality signal currently
in the training set is a function of potency.** The Hill residual is by construction the part of the
screening reading the label does not explain. It is new information about label reliability, it is
not derivable from the existing block, and a training weight is not something the affine pair can
manufacture --- the three conditions item 166 found in every survivor.

Arms:

    база                    равные веса
    вес по остатку          w = 1/(1 + (|r|/s)^2), s --- медиана |r| по ферменту. Ячейка, где
                            два измерения расходятся вдвое сильнее обычного, весит впятеро меньше
    вес перемешанный        те же веса, переставленные между молекулами. Распределение весов то
                            же, связь с молекулой разорвана. ЭТО главный контроль
    вес обратный            w перевёрнут: если и он помогает, работает не надёжность, а
                            неравномерность весов сама по себе

Cells with no screening reading keep weight 1 --- on CYP1A2, CYP2C9 and CYP2D6 there are none, and
on CYP3A4 they are item 129's 530-compound campaign, which must not be down-weighted for lacking a
measurement it never had.

Pre-registered per enzyme against the floors of item 165 (0.0061, 0.0071, 0.0049, 0.0033), not the
macro. The weighted arm must beat both the shuffled and the inverted one; beating the base while
losing to either of them means the weights are not carrying reliability.

Learner pinned to src/ablpairloss.py. Writes results/preds/oof_hill.json. ~30 min for four seeds.
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
PC0 = 4.305
FIT_E = {"CYP1A2": 0.728, "CYP2C9": 0.621, "CYP2D6": 0.867, "CYP3A4": 0.931}


def hill_resid(pi, l2, E):
    """Насколько чтение скрининга отходит от предсказанного меткой через уравнение прибора.

    Считается при ПОДОГНАННОМ E, а не при единице: пункт 170 намерил, что при подогнанном E
    медианный наклон равен 1.00 на трёх ферментах из четырёх, то есть это та параметризация,
    в которой прибор физически осмыслен. Остаток --- в единицах log2fc, как само чтение.
    """
    I = E / (1 + 10 ** (1.0 * (PC0 - pi)))          # наклон фиксирован единицей
    return np.log2(np.clip(1 - I, 1e-6, None)) - l2


def boost(Xtr, ytr, wtr, Xte, seed):
    rng = np.random.default_rng(seed)
    base = float(np.average(ytr, weights=wtr))
    s = np.full(len(ytr), base)
    pred = np.full(len(Xte), base)
    for _ in range(NTREE):
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MAXFEAT,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(Xtr, ytr - s, sample_weight=wtr)
        s = s + LR * t.predict(Xtr)
        pred = pred + LR * t.predict(Xte)
    return pred


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--arms", default="база|вес по остатку|вес перемешанный|вес обратный")
    ap.add_argument("--out", default=RES + "preds/oof_hill.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    arms = a.arms.split("|")

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)

    Y, LO, HI, M, W = {}, {}, {}, {}, {}
    print(f"{'фермент':8s} {'помечено':>9s} {'со скринингом':>14s} {'медиана |r|':>12s} "
          f"{'вес: мин':>9s} {'медиана':>8s}")
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        M[c] = tr[col].notna().to_numpy()
        Y[c] = tr[col].to_numpy(float)
        LO[c] = tr[col + "_conf_low"].to_numpy(float)
        HI[c] = tr[col + "_conf_high"].to_numpy(float)
        sub = sc[sc.enzyme == c].set_index("Molecule_Name")["log2fc_estimate"]
        l2 = sub.reindex(rows.Molecule_Name).to_numpy(float)
        r = hill_resid(Y[c], l2, FIT_E[c])
        got = np.isfinite(r)
        s = np.nanmedian(np.abs(r[got & M[c]])) if (got & M[c]).any() else 1.0
        w = np.ones(len(rows))
        w[got] = 1.0 / (1.0 + (np.abs(r[got]) / max(s, 1e-9)) ** 2)
        W[c] = w
        print(f"{c:8s} {int(M[c].sum()):9d} {int((got & M[c]).sum()):14d} {s:12.3f} "
              f"{w[got & M[c]].min():9.3f} {np.median(w[got & M[c]]):8.3f}")
    print()

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for arm in arms:
            t0 = time.time()
            r = {"seed": seed, "рука": arm}
            for c in CYPS:
                m = M[c]
                y, lo, hi, fi = Y[c][m], LO[c][m], HI[c][m], fold[m]
                if arm == "база":
                    w = np.ones(int(m.sum()))
                else:
                    w = W[c][m].copy()
                    if "перемешанный" in arm:
                        w = w[np.random.default_rng(seed * 17 + 3).permutation(len(w))]
                    elif "обратный" in arm:
                        w = w.max() + w.min() - w
                p = np.zeros(len(y))
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    p[te] = boost(X[m][trn], y[trn], w[trn], X[m][te], seed * 10 + f)
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{arm}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for tag in ("пара", "rho"):
                r[f"MACRO {tag}"] = round(float(np.mean([r[f"{c} {tag}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {arm:20s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    g = df.groupby("рука", sort=False)[["MACRO rho"] + [f"{c} rho" for c in CYPS]].mean()
    print("\n" + g.round(4).to_string())
    if "база" in g.index:
        print(f"\n{'рука':20s} " + "".join(f"{c[3:]+' Δ':>12s}" for c in CYPS))
        for arm in g.index:
            if arm == "база":
                continue
            line = f"{arm:20s} "
            for c in CYPS:
                d = g.loc[arm, f"{c} rho"] - g.loc["база", f"{c} rho"]
                line += f"{d:+11.4f}{'*' if abs(d) > FLOOR[c] else ' '}"
            print(line)
        print("  * --- больше СОБСТВЕННОГО пола фермента (пункт 165)")
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. «вес по остатку» обязан обыграть И перемешанный, И обратный. Обыгрывает базу, но
не перемешанный --- работает неравномерность весов, а не надёжность. Обыгрывает перемешанный,
но не обратный --- знак не тот, и расхождение двух измерений помечает не плохие метки.

Читается поферментно против полов пункта 165. CYP3A4 --- особый случай: его 530 соединений
кампании остаются с весом 1, потому что у них нет скринингового чтения, а не потому что они
надёжны.""")


if __name__ == "__main__":
    main()
