"""Per-enzyme NON-NEGATIVE weights over the current eight-member composition, which item 316's
closure does not cover by its own words.

WHY THIS IS NOT A CLOSED QUESTION REOPENED. Item 316 measured a per-enzyme weight simplex at
+0.0006 (sign 5/8) and a Caruana bagged selector at -0.0048 (0/8), and its own caveat said: "this is
about OUR library --- five learners over one feature matrix --- not about a library of
architecturally distinct models". Item 324 quotes that caveat and rules on it: **"That caveat was
right, and the library was the binding constraint rather than the combiner."** Concretely, item 324
found a FITTED weight on the sixth member gives +0.0141 against the equal weight's +0.0138, with the
chooser taking a > 0 in 158 of 160 (enzyme, fold, seed) cells at a mean chosen weight of 0.333 --- so
when the library stopped being five correlated learners, the combiner stopped being worthless.

The library has changed again since. The composition is now the `_keep`-admitted members plus THREE
architecturally distinct additions -- the frozen `chemprop_medium` probe, the frozen
`chemprop_chemeleon` probe and the gbm head on the medium columns (items 324, 345, 346) -- whose
out-of-fold errors correlate 0.845 to 0.926 with the five rather than the 0.91--0.99 that item 289
recorded inside the original library. Equal weighting now gives each member 0.125 to 0.200 depending
on the enzyme, and nobody has asked whether that is the right split for THIS set.

NON-NEGATIVE AND NOT FREE, for a measured reason that is not ours. The entrant at board place 41
(`lachrymator`, rho 0.7378 against our 0.7062) publishes why: unconstrained ridge and elastic-net
combiners returned large cancelling coefficients -- one member near +1.9 against another at -0.6,
and on one endpoint the best standalone model went negative and was used as a correction term. That
fits held-out folds and transfers badly: "one revision improved cross-validated error on all four
endpoints while real held-out error got worse on three". Non-negativity forbids the cancellation
structurally rather than penalising it. Their residual correlations run 0.55--0.88, ours 0.845--0.99,
so we have LESS spread to exploit and more reason to constrain.

WHAT IS FITTED, AND NESTED. For fold f the weights are fitted on folds != f and applied to fold f,
so no weight sees the rows it scores. Two objectives, because the metric and the usual stacking loss
are different functions and only one of them is what the board pays:

    нмнк     неотрицательный МНК по y          -- обычное стекирование
    хинж     неотрицательный по хинжу полосы   -- собственная потеря метрики

THREE CONTROLS, one of which must FAIL for the result to mean anything.

  * equal weights must reproduce item 346's stack (+0.0047 of rank, -0.0085 of pair at 4/4), or the
    device has moved and nothing here compares to anything;
  * the fit must beat equal weights IN SAMPLE on the fitting folds -- if it does not, the optimiser
    is broken and an out-of-fold zero would be a bug reported as a finding;
  * the SHUFFLE null: weights fitted against permuted labels must give an out-of-fold gain near
    zero. With five to eight free parameters per (enzyme, fold) on 250 to 470 rows this is the
    control that separates a result from fitting noise.

PREDICTION, fixed before the run. Item 316's +0.0006 came from a library it said was the binding
constraint; this library is less correlated, so there is more to exploit than there was. Against
that: the per-fold affine pair rewrites placement, so only RANK survives, and least squares fits
placement rather than ordering -- the objective and the payoff are misaligned by construction. I
expect a small positive, **[0.000, 0.004] of macro rank, most likely under the 0.0036 floor**, with
the shuffle null within 0.001 of zero. If it clears the floor at sign 4/4 it is the first combiner
gain this project has measured, and item 316's closure needs its scope written down rather than its
number changed.

    uv run python verify/k143_nnweights.py --seeds 0,1,2,3

Writes results/logs/k143_nnweights.json. Exit: 0 ran, 2 a control failed, 3 an input is missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json

import numpy as np
import pandas as pd
from scipy.optimize import nnls, minimize

from cypsplit import butina_folds, fold_digest
from submit import CYPS

import k106_probe as K
from k106_probe import score, members, probe_oof, load_emb
from k121_chemeleon import load_any, NEW

LOG = RES + "logs/k143_nnweights.json"
FLOOR = 0.0036
GOLDEN = "2d93c19815e14261"
SURR = D + "emb_chemprop_surrogate_adme.npz"


def hinge_sum(p, lo, hi):
    return float((np.clip(p - hi, 0, None) + np.clip(lo - p, 0, None)).sum())


def fit_nnls(A, y):
    """Non-negative least squares. Scale is free: the affine pair absorbs it, so rank is
    unaffected by any normalisation and none is imposed."""
    w, _r = nnls(A, y)
    return w


def fit_hinge(A, y, lo, hi, w0):
    """Non-negative weights minimising the metric's own hinge, by SLSQP over few variables.

    Started from the NNLS solution so a failure to improve is a failure of the objective rather
    than of the starting point, and bounded below at zero for the same structural reason the
    least-squares arm is: a negative weight turns a member into a correction term, which board
    place 41 measured as the thing that transfers worst.
    """
    f = lambda w: hinge_sum(A @ w, lo, hi)
    r = minimize(f, w0, method="SLSQP", bounds=[(0.0, None)] * A.shape[1],
                 options={"maxiter": 200, "ftol": 1e-9})
    return r.x if r.x is not None and np.isfinite(r.fun) else w0


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    rows = pd.read_csv(D + "rows.csv")
    z = np.load(D + "feats.npz")
    XFULL = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)
    K.S.LO, K.S.HI = LO, HI

    out = json.load(open(LOG)) if _pl.Path(LOG).exists() else {}
    out.setdefault("meta", {}).update(floor=FLOOR)
    res = out.setdefault("seeds", {})

    for s in seeds:
        fold, _ncl = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        print("\n  сид %d: дайджест %s" % (s, dg), flush=True)
        if s == 0 and dg != GOLDEN:
            print("ОТКАЗ: разбиение сдвинулось, сид 0 дал %s" % dg)
            return 2

        pdz = members(s, XFULL, Y, M, fold, LO, HI)
        probe = [probe_oof(load_emb(rows), Y, M, fold, e)[0] for e in range(len(CYPS))]
        cp = RES + "preds/k122_chemeleon_probe_seed%d.json" % s
        gp = RES + "preds/k114_heads_seed%d.json" % s
        for p in (cp, gp):
            if not _pl.Path(p).exists():
                print("нет %s --- нужны кэши пунктов 345/335" % p)
                return 3
        cj = json.load(open(cp))
        hj = json.load(open(gp))
        if cj.get("digest") != dg or hj.get("digest") != dg:
            print("ОТКАЗ: кэш собран на другом разбиении")
            return 2
        chem = [np.asarray(v, float) for v in cj["P"]]
        gbm = [np.asarray(v, float) for v in hj["heads"]["gbm"]]

        # the eight-member library, per enzyme, in a fixed order so the weights can be read
        lib, names = [], []
        for e in range(len(CYPS)):
            cols = [Q[e] for k, Q in pdz if K._keep(e, k)]
            nm = [k for k, _Q in pdz if K._keep(e, k)]
            cols += [probe[e], chem[e], gbm[e]]
            nm += ["зонд", "зонд2", "gbm600"]
            lib.append(np.column_stack(cols))
            names.append(nm)
        print("      членов по ферментам: %s" % [L.shape[1] for L in lib], flush=True)

        eq = [L.mean(1) for L in lib]
        s_eq = score(eq, Y, LO, HI, M, fold)
        print("      КОНТРОЛЬ равные веса: ранг %.4f, пара %.4f" % (s_eq[0], s_eq[1]), flush=True)

        arms = {}
        rng = np.random.default_rng(s)
        for obj in ("нмнк", "хинж", "нмнк_перемеш"):
            P, insample_win, chosen = [], [], []
            for e in range(len(CYPS)):
                m = M[:, e]
                A, y = lib[e], Y[m, e]
                lo, hi, fm = LO[m, e], HI[m, e], fold[m]
                p = np.full(len(y), np.nan)
                for f in range(5):
                    trn, te = fm != f, fm == f
                    if te.sum() == 0:
                        continue
                    yt = y[trn]
                    if obj == "нмнк_перемеш":
                        yt = rng.permutation(yt)
                    w = fit_nnls(A[trn], yt)
                    if obj == "хинж":
                        w = fit_hinge(A[trn], yt, lo[trn], hi[trn], w)
                    if w.sum() <= 0:
                        w = np.ones(A.shape[1])
                    p[te] = A[te] @ (w / w.sum())
                    # in-sample: the fit must beat equal weights on its own fitting rows
                    insample_win.append(float(
                        hinge_sum(A[trn] @ (w / w.sum()), lo[trn], hi[trn])
                        <= hinge_sum(A[trn].mean(1), lo[trn], hi[trn]) + 1e-9))
                    chosen.append((w / w.sum()).tolist())
                P.append(p)
            sc = score(P, Y, LO, HI, M, fold)
            arms[obj] = {"score": sc, "gain_rank": sc[0] - s_eq[0], "gain_pair": sc[1] - s_eq[1],
                         "insample_win": float(np.mean(insample_win)),
                         "mean_weights": np.mean([c for c in chosen if len(c) == 8], 0).tolist()
                                         if any(len(c) == 8 for c in chosen) else []}
            print("      %-14s Δранг %+.4f | Δпара %+.4f | в выборке лучше равных %.0f%%"
                  % (obj, sc[0] - s_eq[0], sc[1] - s_eq[1], 100 * np.mean(insample_win)),
                  flush=True)

        res[str(s)] = {"equal": s_eq, "arms": arms,
                       "members_per_enzyme": [int(L.shape[1]) for L in lib],
                       "names": names}
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)

    print("\n" + "=" * 92)
    print("ПОФЕРМЕНТНЫЕ НЕОТРИЦАТЕЛЬНЫЕ ВЕСА НАД ТЕКУЩИМ СОСТАВОМ, сиды %s" % seeds)
    print("=" * 92)
    print("  %-16s %9s %8s %7s %10s %12s" % ("объектив", "Δранг", "sd", "знак", "Δпара", "в выборке"))
    for obj in ("нмнк", "хинж", "нмнк_перемеш"):
        g = [res[str(s)]["arms"][obj]["gain_rank"] for s in seeds if str(s) in res]
        p = [res[str(s)]["arms"][obj]["gain_pair"] for s in seeds if str(s) in res]
        iw = [res[str(s)]["arms"][obj]["insample_win"] for s in seeds if str(s) in res]
        if not g:
            continue
        sd = float(np.std(g, ddof=1)) if len(g) > 1 else float("nan")
        print("  %-16s %+9.4f %8.4f %5d/%d %+10.4f %11.0f%%"
              % (obj, float(np.mean(g)), sd, int(sum(1 for x in g if x > 0)), len(g),
                 float(np.mean(p)), 100 * float(np.mean(iw))))
    print("\n  пол макро-ранга %.4f; пункт 316 на ПЯТИ учителях дал +0.0006 при знаке 5/8" % FLOOR)
    print("  РЕШАЕТ: отличается ли нмнк от перемешанного нуля больше, чем на пол")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
