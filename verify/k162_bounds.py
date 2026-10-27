"""Add our OWN screened-but-unmeasured compounds as LEFT-CENSORED BOUNDS, on the two isoforms where
the screen can tell an active from an inactive --- the device every open-code entrant above us uses
and the one thing from item 316 this project never built on its own data.

WHY HERE AND WHY NOW. Item 316 records that all three open-code entrants ranked above us enter
inactives as left-censored bounds under a one-sided loss, "where we hold 129 CYP2D6 labels below
pIC50 4.0 against their ~20 000 bounds". Item 388 applied that device to the EXTERNAL NCGC panel
inside the joint fine-tune, and item 389 measured the whole branch dead. Our own screen has none of
the transfer problem that killed it: same laboratory, same protocol, same chemistry.

THE NUMBERS THAT MAKE THIS WORTH A DAY, verified here rather than quoted:

    фермент   pIC50 есть   ниже 4.0   скрин. БЕЗ pIC50   прирост   Youden скрининга (п. 381)
    CYP1A2          1412        232               2963     12.8x               0.814
    CYP2C9          1285        259               3090     11.9x               0.139
    CYP2D6          1493        129               2882     22.3x               0.892
    CYP3A4          2335        944               2570      2.7x               0.151

Item 353 puts 0.0609 of our 0.0909 macro deficit on CYP1A2 and CYP2D6 --- exactly the two with the
fewest inactivity labels AND exactly the two where the screen discriminates. CYP2C9 and CYP3A4 are
NOT touched: item 381 refused to derive them and item 382 explained why.

WHY A BOUND AND NOT A PSEUDO-LABEL, which is the whole point. Item 83 measured the pseudo-label
route and it collapsed on CYP2D6: only 12.6 per cent of its screen-only readings are invertible,
because 87.4 per cent show NO inhibition at all --- median -0.298, noise about zero. You cannot
compute a value from "no inhibition". You CAN state an inequality. **The rows on which item 83's
inversion failed are exactly the rows a bound is made of**, and its selection table says the same
from the other side: on CYP2D6 "curves were run according to the screen", all curve compounds are
screen-active and 5.7 per cent of the rest are. For a pseudo-label that is fatal selection; for a
bound it is the content.

The bound is the protocol's, not an estimate: item 83 records the screen's single concentration as
4.95e-5 M, which is pC0 = 4.305 exactly.

HONEST LIMIT OF THIS IMPLEMENTATION. A true one-sided loss charges only ABOVE the bound.
HistGradientBoosting has no per-row asymmetric loss, so the bound enters as a point target at 4.305
with a weight. That is an APPROXIMATION in the conservative direction --- it also charges for
predicting below 4.305, which a real bound would not. If it pays, the exact version is worth
building; if it loses, the exact one still might, and this file says so rather than claiming a
closure it did not earn.

THREE CONTROLS:
  * weight 0 IS the baseline and is run through the identical code path, so a difference at w=0
    would mean the plumbing moved rather than the data;
  * a NULL arm gives the same rows targets drawn at random from that enzyme's observed labels
    instead of the bound --- if a gain survives that, it is row count and not content;
  * CYP2C9 and CYP3A4 are carried untouched, so any macro movement must come from the two
    isoforms the arm actually touches.
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
from sklearn.ensemble import HistGradientBoostingRegressor

from cypsplit import butina_folds

# Defined here rather than imported from src/submit.py: that module calls tutorial() at import
# time and so demands the organisers' repository, which this file does not otherwise need. The
# order is fixed and matches CLAUDE.md's table.
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]

BOUND = 4.305            # pC0 of the screen's single concentration, 4.95e-5 M (item 83)
ARMED = ["CYP1A2", "CYP2D6"]      # item 381's Youden gate: 0.814 and 0.892 against 0.139 and 0.151
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31, l2_regularization=1.0)


def fit_oof(X, y, m, fold, cens, target, w, seed):
    """OOF over the Butina folds. Censored rows enter TRAINING folds only and are never scored."""
    p = np.full(int(m.sum()), np.nan)
    idx = np.where(m)[0]
    for k in np.unique(fold):
        tr_lab = idx[fold[idx] != k]
        te_lab = idx[fold[idx] == k]
        Xtr, ytr = X[tr_lab], y[tr_lab]
        sw = np.ones(len(tr_lab))
        if w > 0 and cens is not None:
            cr = np.where(cens & (fold != k))[0]
            Xtr = np.vstack([Xtr, X[cr]])
            ytr = np.concatenate([ytr, target[cr]])
            sw = np.concatenate([sw, np.full(len(cr), w)])
        mdl = HistGradientBoostingRegressor(random_state=seed, **KW).fit(Xtr, ytr, sample_weight=sw)
        p[np.searchsorted(idx, te_lab)] = mdl.predict(X[te_lab])
    if not np.isfinite(p).all():
        raise SystemExit("OOF неполон: %d пропусков" % int((~np.isfinite(p)).sum()))
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--weights", default="0,0.1,0.3,1.0")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--out", default=RES + "logs/k162_bounds.json")
    a = ap.parse_args()

    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr["%s_pIC50_direct_inhibition" % c].to_numpy(float) for c in CYPS], 1)
    M = np.isfinite(Y)
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    screened = rows.Molecule_Name.isin(set(sc.Molecule_Name.unique())).to_numpy()
    fold, ncl = butina_folds(list(rows.SMILES))
    fold = np.asarray(fold)
    print("  X %s, кластеров %d, скринировано %d" % (X.shape, ncl, int(screened.sum())), flush=True)
    print("  %-8s %10s %10s %14s" % ("фермент", "размечено", "ниже 4.0", "границ доступно"))
    CENS = {}
    for e, c in enumerate(CYPS):
        CENS[c] = screened & ~M[:, e]
        print("  %-8s %10d %10d %14d%s"
              % (c, int(M[:, e].sum()), int((M[:, e] & (Y[:, e] < 4.0)).sum()),
                 int(CENS[c].sum()), "   <- арм" if c in ARMED else ""), flush=True)

    ws = [float(x) for x in a.weights.split(",")]
    rng = np.random.default_rng(0)
    res, base = {}, {}
    print("\n  прогон (каждая строка = %d сидов)" % a.seeds, flush=True)
    for mode in ("граница", "НУЛЬ"):
        for w in ws:
            if w == 0 and mode == "НУЛЬ":
                continue
            per = {}
            for e, c in enumerate(CYPS):
                if c not in ARMED:
                    continue
                tgt = np.full(len(rows), BOUND)
                if mode == "НУЛЬ":
                    tgt = rng.choice(Y[M[:, e], e], size=len(rows), replace=True)
                per[c] = [float(spearmanr(
                    fit_oof(X, Y[:, e], M[:, e], fold, CENS[c] if w > 0 else None, tgt, w, s),
                    Y[M[:, e], e]).statistic) for s in range(a.seeds)]
            key = "%s w=%.2f" % (mode, w)
            res[key] = per
            print("  %-16s %s" % (key, "  ".join(
                "%s %.4f" % (c, np.mean(per[c])) for c in ARMED)), flush=True)
            if w == 0:
                base = {c: list(per[c]) for c in ARMED}

    print("\n  ИЗМЕНЕНИЕ ПРОТИВ БАЗЫ (вес 0, тот же код)")
    print("  %-16s %12s %12s %14s %8s" % ("арм", "CYP1A2", "CYP2D6", "вклад в макро", "знак"))
    for key, per in res.items():
        if key == "граница w=0.00":
            continue
        ds = {c: np.mean(per[c]) - np.mean(base[c]) for c in ARMED}
        sg = sum(1 for c in ARMED for i in range(a.seeds) if per[c][i] > base[c][i])
        print("  %-16s %+12.4f %+12.4f %+14.4f %5d/%d"
              % (key, ds["CYP1A2"], ds["CYP2D6"], sum(ds.values()) / 4.0, sg, 2 * a.seeds))
    json.dump({"bound": BOUND, "armed": ARMED, "res": res, "base": base},
              open(a.out, "w"), ensure_ascii=False, indent=1)
    print("\n  записано %s  (пол макро 0.0036)" % a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
