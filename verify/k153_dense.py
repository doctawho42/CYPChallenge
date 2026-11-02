"""Re-read every measurement we already have, restricted to the DENSEST part of our own data.

Item 375's geometry (k152) shows the blind test sits at a median nearest-neighbour similarity of
0.5873 while leave-one-out on our own 4905 compounds reaches only 0.4500. So a test-matched split
cannot be built from our data: the ceiling is below the target. What CAN be done costs nothing --
the out-of-fold predictions for every member and arm are committed, and k152 wrote each compound's
out-of-fold nearest-neighbour similarity, so each measurement can be re-read on the subpopulation
whose geometry is closest to the test's.

This does not create the test's regime. It asks a weaker question that is still decisive: does the
ORDERING of our methods depend on local density? If the ranking in the dense quartile differs from
the ranking overall, then conclusions drawn under the Butina split do not transfer to a denser test
and the journal's negatives have to be re-read. If the ordering is stable, the regime mismatch
costs us absolute score but not decisions.

CONTROLS, both of which must succeed:
  * the FULL-population figures must reproduce the committed numbers of items 348, 372 and 374,
    or the subsetting machinery is not measuring what those items measured;
  * a RANDOM subset of the same size must show no ordering change, or any change in the dense
    quartile is just the smaller sample.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from cypsplit import butina_folds, fold_digest
from submit import CYPS
import k106_probe as K
from k106_probe import members, probe_oof, load_emb

NN = RES + "preds/k152_nn_butina.npy"


def main():
    rows = pd.read_csv(D + "rows.csv")
    z = np.load(D + "feats.npz")
    XF = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr["%s_pIC50_direct_inhibition" % c].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)
    nn = np.load(NN)
    assert len(nn) == len(rows), (len(nn), len(rows))

    fold, _ = butina_folds(list(rows.SMILES))
    dg = fold_digest(fold)
    print("  сид 0, дайджест %s; соседство из %s" % (dg, NN.split("/")[-1]))
    pdz = members(0, XF, Y, M, fold, None, None)
    probe = [probe_oof(load_emb(rows), Y, M, fold, e)[0] for e in range(len(CYPS))]
    chem = [np.asarray(v, float) for v in
            json.load(open(RES + "preds/k122_chemeleon_probe_seed0.json"))["P"]]
    gbm = [np.asarray(v, float) for v in
           json.load(open(RES + "preds/k114_heads_seed0.json"))["heads"]["gbm"]]
    TE = [np.asarray(v, float) for v in
          json.load(open(RES + "preds/k147_tabicl_E100.json"))["preds"]["E_seed0"]]
    TS = [np.asarray(v, float) for v in
          json.load(open(RES + "preds/k147_tabicl_S100.json"))["preds"]["S_seed0"]]
    cands = [(k, Q) for k, Q in pdz] + [("зонд", probe), ("chemeleon", chem),
                                        ("gbm600", gbm), ("TabICL-E", TE), ("TabICL-S", TS)]
    kept = [[Q[e] for k, Q in pdz if K._keep(e, k)] for e in range(len(CYPS))]
    arms = {"шесть (отгружено)": [np.mean(kept[e] + [probe[e]], 0) for e in range(len(CYPS))],
            "восемь (п.346)": [np.mean(kept[e] + [probe[e], chem[e], gbm[e]], 0)
                               for e in range(len(CYPS))],
            "семь (п.374)": [np.mean(kept[e] + [TE[e], TS[e]], 0) for e in range(len(CYPS))]}

    def solo(P, sel):
        """Mean per-enzyme Spearman on the selected rows. sel is a boolean over ALL 4905."""
        r = []
        for e in range(len(CYPS)):
            m = M[:, e]
            s = sel[m]
            if s.sum() < 60:
                return None
            r.append(float(spearmanr(Y[m, e][s], P[e][s]).statistic))
        return float(np.mean(r))

    rng = np.random.default_rng(0)
    q75 = float(np.quantile(nn, 0.75))
    regimes = {"всё (контроль)": np.ones(len(nn), bool),
               "плотная четверть": nn >= q75,
               "разрежённая четверть": nn <= float(np.quantile(nn, 0.25))}
    # the control that makes a change in the dense quartile readable
    for i in range(3):
        idx = rng.choice(len(nn), int((nn >= q75).sum()), replace=False)
        b = np.zeros(len(nn), bool); b[idx] = True
        regimes["случайная четверть %d" % (i + 1)] = b

    print("  порог плотной четверти: соседство >= %.4f (медиана режима %.4f);"
          " у теста медиана 0.5873" % (q75, float(np.median(nn[nn >= q75]))))
    print("\n  %-18s %s" % ("", "".join("%12s" % r[:12] for r in regimes)))
    out = {}
    for name, P in list(arms.items()) + cands:
        vals = {r: solo(P, s) for r, s in regimes.items()}
        out[name] = vals
        print("  %-18s %s" % (name, "".join("%12.4f" % v if v is not None else "%12s" % "—"
                                            for v in vals.values())))
    # ordering: do the candidates rank differently in the dense quartile?
    nm = [k for k, _ in cands]
    base = [out[k]["всё (контроль)"] for k in nm]
    print("\n  ПОРЯДОК кандидатов (Spearman порядка против порядка на всём):")
    for r in list(regimes)[1:]:
        v = [out[k][r] for k in nm]
        print("    %-22s %+.4f" % (r, spearmanr(base, v).statistic))
    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump({"q75": q75, "regimes": {k: int(v.sum()) for k, v in regimes.items()},
               "scores": out}, open(RES + "logs/k153_dense.json", "w"),
              ensure_ascii=False, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
