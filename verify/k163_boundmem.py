"""Left-censored screening bounds inside the REAL per-enzyme member, substituted into the shipped
composition. Item 162 measured the mechanism on a single-stage reimplementation and it paid
(+0.0125 macro, null at 0/6); this file asks whether it survives the member that actually ships.

WHY A SECOND FILE. `verify/k162_bounds.py` fitted one HistGB per enzyme and measured bounds against
it. That is NOT the shipped member: the shipped one is two-stage (item 164), a base OOF followed by
a dead-zone pass that clips the prediction into the credible band and refits under absolute_error.
The single-stage version sits 0.04 to 0.05 of rank below it, correlating 0.876 with it at identical
folds and pins -- so k162's gain is a real measurement of a DIFFERENT member, and the question of
whether a stronger member has the same headroom is open. Items 269, 273 and 281 are three
member-level gains that died at exactly this step, and `verify/k87_maxfeat2.py` is the file that
measured the last one; its structure is reused here deliberately.

THE ARM. Screened-but-unmeasured compounds enter the BASE fit as a point target at pIC50 = 4.305 --
pC0 of the screen's single concentration, 4.95e-5 M, recorded in item 83 -- carrying weight w. The
dead-zone pass is left exactly as it ships: its target is the clipped OOF of the LABELLED rows only,
so nothing about the shipped second stage changes and any movement is attributable to the bounds.

Armed on CYP1A2 and CYP2D6 only. Item 381 measured the screen's Youden at 0.814 and 0.892 there
against 0.139 and 0.151 on CYP2C9 and CYP3A4, and items 282-285 leave the per-enzyme member in the
composition only on those two cells anyway (CYP1A2 via пофе+GP+ствол, CYP2D6 via all five).

FLOORS, from k87 which measured the same two cells: CYP1A2 0.0061, CYP2D6 0.0049.

CONTROLS:
  * w = 0 must reproduce the CACHED member to a tight tolerance -- it is the same code path with no
    rows added, so a mismatch means the rebuild differs from the cache for some other reason and
    nothing below can be read;
  * a NULL arm gives the same rows targets drawn from the enzyme's observed labels, and must lose;
  * CYP2C9 and CYP3A4 are carried from the cache untouched.
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

from cypsplit import butina_folds
import submit as S
from submit import CYPS, _keep, DZ_KW

MEM = ["поферментно", "пул", "GP", "гребневая", "ствол"]
BOUND = 4.305
ARMED = ["CYP1A2", "CYP2D6"]
FLOOR = {"CYP1A2": 0.0061, "CYP2D6": 0.0049}
BASE_KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
               l2_regularization=1.0, random_state=0)


def perenz(X, Y, LO, HI, mask, fold, cens, target, w):
    """The shipped two-stage per-enzyme member, with censored rows added to the BASE fit only."""
    out = []
    for e in range(4):
        m = mask[:, e]
        Xi, yy, fi = X[m], Y[m, e], fold[m]
        P = np.zeros(len(yy))
        for f in range(5):
            a, b = fi != f, fi == f
            if b.sum() == 0:
                continue
            Xtr, ytr = Xi[a], yy[a]
            sw = np.ones(int(a.sum()))
            if w > 0 and cens is not None and CYPS[e] in ARMED:
                cr = np.where(cens[CYPS[e]] & (fold != f))[0]
                Xtr = np.vstack([Xtr, X[cr]])
                ytr = np.concatenate([ytr, target[CYPS[e]][cr]])
                sw = np.concatenate([sw, np.full(len(cr), w)])
            P[b] = HistGradientBoostingRegressor(**BASE_KW).fit(
                Xtr, ytr, sample_weight=sw).predict(Xi[b])
        t = np.clip(P, LO[m, e], HI[m, e])          # dead zone: shipped, labelled rows only
        Q = np.zeros(len(yy))
        for f in range(5):
            a, b = fi != f, fi == f
            if b.sum() == 0:
                continue
            Q[b] = HistGradientBoostingRegressor(**DZ_KW).fit(Xi[a], t[a]).predict(Xi[b])
        out.append(Q)
    return out


def compose(members, e):
    return np.mean([np.asarray(members[k][e], float) for k in MEM if _keep(e, k)], axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--weights", default="0,0.3,1.0")
    ap.add_argument("--seeds", default="4,5,6,7")
    ap.add_argument("--out", default=RES + "logs/k163_boundmem.json")
    a = ap.parse_args()

    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr["%s_pIC50_direct_inhibition" % c].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr["%s_pIC50_direct_inhibition_conf_low" % c].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr["%s_pIC50_direct_inhibition_conf_high" % c].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(Y)
    S.LO, S.HI = LO, HI
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    screened = rows.Molecule_Name.isin(set(sc.Molecule_Name.unique())).to_numpy()
    CENS = {c: screened & ~mask[:, e] for e, c in enumerate(CYPS)}
    TGT = {c: np.full(len(rows), BOUND) for c in CYPS}
    rng = np.random.default_rng(0)
    NULLT = {c: rng.choice(Y[mask[:, e], e], size=len(rows), replace=True)
             for e, c in enumerate(CYPS)}
    print("  границ: %s" % "  ".join("%s %d" % (c, int(CENS[c].sum())) for c in ARMED), flush=True)

    ws = [float(x) for x in a.weights.split(",")]
    seeds = [int(x) for x in a.seeds.split(",")]
    out = {}
    print("\n  %3s %-22s %9s %9s %9s" % ("сид", "арм", "CYP1A2", "CYP2D6", "время"), flush=True)
    for s in seeds:
        mem = json.load(open(RES + "preds/members_seed%d.json" % s))
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        fold = np.asarray(fold)
        ref = {c: float(spearmanr(Y[mask[:, e], e], compose(mem, e)).statistic)
               for e, c in enumerate(CYPS)}
        print("  %3d %-22s %9.4f %9.4f   (кэш)" % (s, "подаётся", ref["CYP1A2"], ref["CYP2D6"]),
              flush=True)
        for mode, tgt in (("граница", TGT), ("НУЛЬ", NULLT)):
            for w in ws:
                if w == 0 and mode == "НУЛЬ":
                    continue
                t0 = time.time()
                new = perenz(X, Y, LO, HI, mask, fold, CENS if w > 0 else None, tgt, w)
                m2 = dict(mem); m2["поферментно"] = [list(map(float, v)) for v in new]
                got = {c: float(spearmanr(Y[mask[:, e], e], compose(m2, e)).statistic)
                       for e, c in enumerate(CYPS)}
                key = "%s w=%.2f" % (mode, w)
                out.setdefault(key, {}).setdefault(str(s), {})
                for c in ARMED:
                    out[key][str(s)][c] = got[c] - ref[c]
                print("  %3d %-22s %+9.4f %+9.4f  %6.0f с"
                      % (s, key, got["CYP1A2"] - ref["CYP1A2"], got["CYP2D6"] - ref["CYP2D6"],
                         time.time() - t0), flush=True)

    print("\n  СВОДКА: изменение СОСТАВА против кэша, среднее по %d сидам" % len(seeds))
    print("  %-22s %11s %6s %11s %6s" % ("арм", "CYP1A2", "знак", "CYP2D6", "знак"))
    for key, per in out.items():
        d1 = [per[str(s)]["CYP1A2"] for s in seeds]
        d2 = [per[str(s)]["CYP2D6"] for s in seeds]
        print("  %-22s %+11.4f %3d/%d %+11.4f %3d/%d"
              % (key, np.mean(d1), sum(x > 0 for x in d1), len(d1),
                 np.mean(d2), sum(x > 0 for x in d2), len(d2)))
    print("\n  полы: CYP1A2 %.4f, CYP2D6 %.4f (пункт 286 / k87)" % (FLOOR["CYP1A2"], FLOOR["CYP2D6"]))
    json.dump({"bound": BOUND, "floors": FLOOR, "res": out}, open(a.out, "w"),
              ensure_ascii=False, indent=1)
    print("  записано %s" % a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
