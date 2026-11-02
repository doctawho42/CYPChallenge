"""What transform of the SHIPPED CYP2D6 column lowers its ST-RAE, scored against BAND STRUCTURE at
the test's own composition --- the one cell item 367 closed with a macro-matched argument that could
not see it.

WHY THIS IS NOT ITEM 367 AGAIN. That item killed a GLOBAL spread and recorded three errors, all of
which this file is built to avoid:
  * error one, the curve's x-axis was not the flag's argument: `--spread s` reweights our labels and
    refits, so the achieved expansion (1.769 / 1.327 / 1.569 / 1.191) was not the nominal one. Here
    the grid multiplies the prediction column DIRECTLY, so the axis means exactly what it says.
  * error two, the target was the test's LABEL SPREAD, from item 354's sd(y) >= MAE/sqrt(1-R2).
    ST-RAE does not reward matching a marginal; it pays zero inside each compound's band. Here the
    objective IS the band hinge and nothing else.
  * error three, the bands are three to seven times narrower than the label spread, so a transform
    aimed at the marginal overshoots them. Here the population is reweighted to the test's own
    composition, where item 384 measured the two band regimes: on CYP2D6 a median half-width of
    0.1306 for actives against 1.2958 for low-activity compounds, a ratio of 9.9.

WHAT THE BOARD SAYS TODAY, which is why one isoform is worth re-opening. Among the 10 entrants whose
CYP2D6 Spearman is within 0.02 of ours, every one beats our CYP2D6 ST-RAE, our MAE of 1.2207 is the
worst of the eleven and our R2 of 0.1453 the lowest --- three independent columns agreeing. The
specificity control is that the same comparison on CYP1A2 puts us at the field median (MAE 0.7584
against 0.7472, R2 0.5485 against 0.5539), so the signature is CYP2D6 and not us generally. And the
mechanism is in the shipped file: CYP2D6 carries sd 0.462 against 0.76 to 0.94 on the other three,
which item 366 already called "under-corrected, worst on CYP2D6 by a factor of nearly three" before
the fix was spent on a global flag.

THE SIMULATION. OOF composition predictions are rank-mapped onto the SHIPPED column's own marginal,
so the starting point is exactly what sits on the board; then a 2-D grid of (shift, multiplier) is
scored by the organisers' ST-RAE against each compound's own band, with the population REWEIGHTED so
the active fraction matches item 381's derived 0.377 for CYP2D6 rather than our training set's.

REGISTERED FALSIFICATION, written before the run. If the optimum for CYP2D6 lands at multiplier
within 1.00 +- 0.05 and shift within +-0.05 --- that is, where we already stand --- then the board's
0.1454 is a rank PERMUTATION we do not have and not placement, and this proposal is dead.

CONTROLS:
  * a multiplier is strictly increasing, so Spearman must not move at all: checked to 1e-12, the
    same control k131 carries;
  * at (0, 1.00) the score must equal the rank-mapped column's own score, since the transform is
    the identity there;
  * the three other isoforms are computed too, and CYP1A2 must NOT show the same headroom --- the
    board says we are at field median there, so a simulation that finds a large gain on CYP1A2 is
    measuring something other than what the board sees.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from cypsplit import butina_folds

from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
MEM = ["поферментно", "пул", "GP", "гребневая", "ствол"]
SOLO = {"CYP1A2": ["поферментно", "GP", "ствол"], "CYP2C9": ["GP", "ствол"],
        "CYP2D6": MEM, "CYP3A4": ["GP", "ствол"]}
# item 381's derived active fraction of the 750, per isoform; only CYP1A2 and CYP2D6 were derived
PACT = {"CYP1A2": 0.355, "CYP2D6": 0.377}
FLOOR_ACT = 4.0


def _hinge(x, lo, hi):
    return np.maximum(x - hi, 0.0) + np.maximum(lo - x, 0.0)


def wstrae(y, p, lo, hi, w):
    """ST-RAE with a population weight.

    The denominator is NOT sum|y - mean(y)|. The organisers' docstring is explicit: "the naive
    baseline in the denominator (a constant predictor at mean(y_true)) is put through THE SAME
    soft-thresholding as the model's predictions", so both halves are band hinges and a constant
    predictor scores zero on every compound whose band already contains the mean. Getting this
    wrong is what the identity control below caught: 0.8108 against the organisers' 1.2392."""
    mu = float(np.average(y, weights=w))
    num = (w * _hinge(p, lo, hi)).sum()
    den = (w * _hinge(np.full_like(y, mu), lo, hi)).sum()
    return float(num / den)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sub", default=RES + "submission/activity_submission_rankmap.csv")
    ap.add_argument("--seeds", type=int, default=4)
    ap.add_argument("--out", default=RES + "logs/k164_d6place.json")
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr["%s_pIC50_direct_inhibition" % c].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr["%s_pIC50_direct_inhibition_conf_low" % c].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr["%s_pIC50_direct_inhibition_conf_high" % c].to_numpy(float) for c in CYPS], 1)
    M = np.isfinite(Y)
    sub = pd.read_csv(a.sub)

    SH = np.arange(-1.2, 1.2001, 0.05)
    MU = np.arange(0.6, 3.0001, 0.05)
    out = {}
    print("  КОНТРОЛЬ соответствия ST-RAE: наша взвешенная функция против функции организаторов")
    fold, _ = butina_folds(list(rows.SMILES))
    print("\n  %-8s %9s %9s %9s %8s %9s %9s %9s"
          % ("фермент", "sd подач.", "старт", "сдвиг*", "x*", "в выборке", "ВНЕ фолда", "ВЫИГРЫШ"))
    for e, c in enumerate(CYPS):
        m = M[:, e]
        y, lo, hi = Y[m, e], LO[m, e], HI[m, e]
        mem = json.load(open(RES + "preds/members_seed0.json"))
        P = np.mean([np.asarray(mem[k][e], float) for k in SOLO[c]], 0)
        # put OOF on the SHIPPED column's marginal, so the starting point is the board's own
        col = sub["%s_pIC50_direct_inhibition" % c].to_numpy(float)
        q = np.empty(len(P))
        q[np.argsort(P, kind="stable")] = np.quantile(col, np.linspace(0, 1, len(P)))
        # control: the identity leg must agree with the organisers' unweighted function
        w1 = np.ones(len(y))
        a1, a2 = wstrae(y, q, lo, hi, w1), float(strae(y, q, y_true_upper=hi, y_true_lower=lo))
        if abs(a1 - a2) > 1e-9:
            raise SystemExit("КОНТРОЛЬ: взвешенная %.6f против организаторской %.6f" % (a1, a2))
        # reweight the population to the test's derived composition where item 381 derived it
        act = y >= FLOOR_ACT
        if c in PACT:
            w = np.where(act, PACT[c] / max(act.mean(), 1e-9),
                         (1 - PACT[c]) / max(1 - act.mean(), 1e-9))
        else:
            w = w1
        base = wstrae(y, q, lo, hi, w)
        ctr = float(q.mean())

        def opt(sel):
            """Best (shift, mult) on the rows `sel`, by the same weighted band hinge."""
            bb = (np.inf, 0.0, 1.0)
            for s in SH:
                for k in MU:
                    v = wstrae(y[sel], ctr + k * (q[sel] - ctr) + s, lo[sel], hi[sel], w[sel])
                    if v < bb[0]:
                        bb = (v, float(s), float(k))
            return bb

        best = opt(np.ones(len(y), bool))
        # The in-sample optimum of a 2-D grid is selection-contaminated by construction -- item
        # 316 called its own such number "pure selection contamination" and item 357 says the
        # floor is "in sample, so optimistic by construction". Fit the pair on four folds and
        # score it on the fifth, which is the only number that may be quoted as a gain.
        fo = np.asarray(fold)[m]
        qo = np.empty(len(y))
        par = []
        for f in range(5):
            a_, b_ = fo != f, fo == f
            if b_.sum() == 0:
                continue
            _, s_, k_ = opt(a_)
            par.append((s_, k_))
            qo[b_] = ctr + k_ * (q[b_] - ctr) + s_
        oof = wstrae(y, qo, lo, hi, w)
        r0 = float(spearmanr(y, q).statistic)
        r1 = float(spearmanr(y, ctr + best[2] * (q - ctr) + best[1]).statistic)
        if abs(r0 - r1) > 1e-12:
            raise SystemExit("КОНТРОЛЬ: преобразование сдвинуло ро на %.3e" % abs(r0 - r1))
        out[c] = {"sub_sd": float(col.std(ddof=1)), "start": base, "best": best[0],
                  "shift": best[1], "mult": best[2], "gain_insample": base - best[0],
                  "oof": oof, "gain_oof": base - oof, "folds": par,
                  "pact": PACT.get(c), "rho_moved": abs(r0 - r1)}
        print("  %-8s %9.3f %9.4f %9.2f %8.2f %9.4f %9.4f %9.4f"
              % (c, col.std(ddof=1), base, best[1], best[2], base - best[0], oof, base - oof),
              flush=True)

    print("\n  КОНТРОЛЬ монотонности: ро не сдвинулось ни на одном ферменте (допуск 1e-12)")
    d6 = out["CYP2D6"]
    dead = (abs(d6["mult"] - 1.0) <= 0.05 and abs(d6["shift"]) <= 0.05) or d6["gain_oof"] <= 0
    print("\n  ЗАРЕГИСТРИРОВАННОЕ ОПРОВЕРЖЕНИЕ: оптимум CYP2D6 при x=1.00+-0.05 и сдвиге +-0.05")
    print("  получено: x=%.2f, сдвиг=%+.2f  ->  %s"
          % (d6["mult"], d6["shift"], "ПРЕДЛОЖЕНИЕ МЕРТВО" if dead else "рычаг есть"))
    print("  КОНТРОЛЬ специфичности ВНЕ ФОЛДА: CYP1A2 %+.4f против CYP2D6 %+.4f"
          % (out["CYP1A2"]["gain_oof"], d6["gain_oof"]))
    print("  подобранные по фолдам (сдвиг, x) для CYP2D6: %s"
          % "  ".join("%+.2f/%.2f" % t for t in d6["folds"]))
    json.dump(out, open(a.out, "w"), ensure_ascii=False, indent=1)
    print("  записано %s" % a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
