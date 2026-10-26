"""Does the joint fine-tune pay as a MEMBER of the shipped composition, which is the only question
that decides its fate --- a solo rank below ours does not settle it, and a solo rank above ours
would not either.

Item 386 diagnosed our ordering deficit as activity-cliff resolution: pairwise sign accuracy decays
monotonically with similarity on CYP1A2 and CYP2D6, reaching 0.533 at the similarity the test
occupies. Fingerprint members cannot fix that, so `src/jointft.py` was built on a message-passing
representation. Its solo rank is below the shipped composition's 0.6417 and that is EXPECTED --- the
composition is five members over engineered features and briford, the best readable entrant,
averages four nets. The live question is whether a representation this different carries error the
composition does not already have.

THE RULE THIS FILE OBEYS, from the memory of being burnt. Item 373's bar --- solo rank 0.6480 --- is
NECESSARY AND NOT SUFFICIENT (item 377): a candidate clearing it still failed because it correlated
0.956 and 0.969 with the members it was assembled from. The converse has never been tested here, and
it is what this file tests: a candidate BELOW the bar whose errors are genuinely decorrelated. The
honest expectation is still failure --- items 176, 182, 191, 194 and 356 are five member-level gains
that died over the ensemble --- so the floor is 0.0036 of macro rank and the sign must hold.

Reads results/preds/members_seed{s}.json for the shipped five and the fine-tune's own OOF; prints
solo rank, error correlation against the composition, and the nested delta with a sign count.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES

import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from submit import CYPS

ALL5 = ["поферментно", "пул", "GP", "гребневая", "ствол"]
# The shipped composition is NOT a flat mean of five. `results/submission/stack374.meta.json`
# records a per-enzyme SOLO structure from items 282-285, and measuring against the flat mean is
# exactly item 339's defect -- a real gain stated against an arm that does not ship, which
# overstated a TDI result by about a third.
SOLO = {"CYP1A2": ["поферментно", "GP", "ствол"],
        "CYP2C9": ["GP", "ствол"],
        "CYP2D6": ALL5,
        "CYP3A4": ["GP", "ствол"]}
FLOOR = 0.0036


def macro(P, Y, M):
    return float(np.mean([spearmanr(P[e], Y[M[:, e], e]).statistic for e in range(4)]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ft", required=True, nargs="+", help="jointft json, по одному на сид")
    ap.add_argument("--seeds", type=int, default=4)
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr["%s_pIC50_direct_inhibition" % c].to_numpy(float) for c in CYPS], 1)
    M = np.isfinite(Y)

    print("  состав по пунктам 282-285: %s"
          % ", ".join("%s %d чл." % (c, len(SOLO[c])) for c in CYPS))
    print("  %-26s %9s %9s %11s %9s" % ("арм", "одиночн.", "корр.", "+1 членом", "знак"))
    base_all, ft_all, nest_all = [], [], []
    for si, path in enumerate(a.ft[: a.seeds]):
        mem = json.load(open(RES + "preds/members_seed%d.json" % si))
        nmemb = [len(SOLO[c]) for c in CYPS]
        P5 = [np.mean([np.asarray(mem[m][e], float) for m in SOLO[CYPS[e]]], 0) for e in range(4)]
        ft = json.load(open(path))
        # the fine-tune's folds are keyed by fold index with their own row indices
        F = [np.full(int(M[:, e].sum()), np.nan) for e in range(4)]
        pos = [np.where(M[:, e])[0] for e in range(4)]
        for k, rec in ft.items():
            idx = np.asarray(rec["idx"])
            pr = np.asarray(rec["pred"], float)
            for e in range(4):
                sel = M[idx, e]
                if not sel.any():
                    continue
                where = np.searchsorted(pos[e], idx[sel])
                F[e][where] = pr[sel, e]
        miss = [int(np.isnan(F[e]).sum()) for e in range(4)]
        if any(miss):
            print("  сид %d: НЕ ПОЛНЫЙ OOF, пропусков %s — арм пропущен" % (si, miss))
            continue
        b = macro(P5, Y, M)
        f = macro(F, Y, M)
        # each enzyme keeps its OWN member count, so the new member enters with the right weight
        n6 = macro([(nmemb[e] * P5[e] + F[e]) / (nmemb[e] + 1) for e in range(4)], Y, M)
        corr = float(np.mean([
            spearmanr(P5[e] - Y[M[:, e], e], F[e] - Y[M[:, e], e]).statistic for e in range(4)]))
        base_all.append(b); ft_all.append(f); nest_all.append(n6)
        print("  сид %-22d %9.4f %9.3f %11.4f" % (si, f, corr, n6 - b))
    if not base_all:
        raise SystemExit("ни одного полного OOF — измерять нечего")
    d = np.array(nest_all) - np.array(base_all)
    print("\n  состав как есть            %.4f (sd %.4f)" % (np.mean(base_all), np.std(base_all)))
    print("  дообучка одиночным рангом  %.4f (sd %.4f), планка пункта 373 = 0.6480"
          % (np.mean(ft_all), np.std(ft_all)))
    print("  ВЛОЖЕНА ДОП. ЧЛЕНОМ        %+.4f, знак %d/%d, пол %.4f"
          % (d.mean(), int((d > 0).sum()), len(d), FLOOR))
    print("  вердикт: %s" % ("ПЛАТИТ" if d.mean() > FLOOR and (d > 0).all() else "НЕ ПЛАТИТ"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
