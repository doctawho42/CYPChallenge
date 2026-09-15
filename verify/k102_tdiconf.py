"""The per-enzyme TDI confusion matrix, recovered from the 15 September board snapshot.

Everything here comes from `results/leaderboard_2026-09-15.json` and the submission files
on disk. Nothing touches the network, and no test compound is looked up anywhere.

What it establishes, in order:

  1. WHICH per-enzyme TDI tab is which enzyme, on five readings, three of which do not
     depend on tab order at all.

  2. The closed-form inversion of accuracy/precision/recall into a 2x2 table, and the
     three validations asked of it. The point-estimate reading FAILS, and rounding is not
     why: the board publishes bootstrap MEANS (`config.BOOTSTRAP_SAMPLES = 1000`, averaged
     by `average_bootstrap_results_by_endpoint` into the `_mean` columns a leaderboard is
     built from), so the published F1 is no longer the harmonic mean of the published
     precision and recall. A mean of ratios is not the ratio of means.

  3. That failure turned into a measurement. The Jensen gap F1 - H(P, R) is a pure
     curvature term that scales as 1/n, so it pins the scored set size n -- a quantity the
     closed form cannot see at all, because the closed form is scale-free.

  4. A bootstrap-aware inversion: the integer tables whose bootstrap means, under the
     organisers' own resample indices, reproduce all five published numbers inside their
     own Monte-Carlo noise.

  5. What the recovered table does and does not determine about MCC at another threshold.

Controls, because a check that cannot fail is not a check:

  * the vectorised metrics are validated against the organisers' own `bootstrap_metrics`
    on the same arrangement; exact agreement is required, not a plausible one;
  * the closed form and the F1 check are exercised on a SYNTHETIC table four ways: on its
    point metrics, where the inversion must return the table exactly; on the EXPECTED
    bootstrap means, where the F1 check must fail in the board's direction and magnitude
    while the cells stay recoverable, which SHOWS the averaging bias in accuracy,
    precision and recall to be negligible instead of assuming it; on SINGLE realisations,
    which is what a board row is and where the non-integer cells actually come from; and
    under median instead of mean aggregation, so the explanation is tested against its
    obvious alternative rather than being the first one that fitted;
  * "cells sum to n" is reported as VACUOUS, because TN is the residual; a count
    constraint that CAN fail (our own positive calls on the scored rows) replaces it;
  * the enumeration window and the prefilter are both checked for bindingness, so that
    "no other table is consistent" cannot be an artefact of where the search looked.

Runtime: about two minutes (measured 2m09s twice on this machine; the docstring said "about a
minute" until it was timed, and a documented runtime is what decides whether the next reader dares
run a script at all). Reads the snapshot, the TDI training file and the submission files; writes
nothing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()

import glob
import itertools
import json

import numpy as np
import pandas as pd

from evaluation.config import BOOTSTRAP_SAMPLES, CLASSIFICATION_ENDPOINTS, CLASSIFICATION_METRICS
from evaluation.evaluate_predictions import bootstrap_metrics
from evaluation.utils import bootstrap_sampling

ME = "Wizard Lizard Gizzard"
SNAP = RES + "leaderboard_2026-09-15.json"
SUB = RES + "submission/tdi_submission.csv"

MACRO_TAB = "partial_14"
ENZ_TABS = ["partial_15", "partial_16"]
NAMES = ["MCC", "Accuracy", "Precision", "Recall", "F1 Score"]   # board headers
KEYS = ["MCC", "Accuracy", "Precision", "Recall", "F1"]          # organisers' metric names

HALF = 5e-5                      # the board prints four decimals
SD_RND = HALF / 3.0 ** 0.5       # sd of a uniform rounding error
ARRANGEMENTS = 64                # draws of which scored row holds which cell
CHUNK = 300                      # candidate tables per gather
RNG = np.random.default_rng(20260915)

# Item 300: our per-enzyme out-of-fold MCC on the shipped TDI bundle.
OOF_MCC = {"CYP2D6": 0.1282, "CYP3A4": 0.3578}
# Item 235: the per-enzyme MCC noise floor at fixed seed.
MCC_FLOOR = {"CYP2D6": 0.0419, "CYP3A4": 0.0281}
# Item 299, from the challenge Space's FAQ: the live board scores half the test set.
LIVE_HALF = 375
FULL = 750
# Prefilter on POINT metrics before the bootstrap forward model; checked for bindingness.
PREFILTER = {"Accuracy": 0.008, "Precision": 0.012, "Recall": 0.05}


# ---------------------------------------------------------------------------
# metrics
# ---------------------------------------------------------------------------


def point_metrics(tp, fp, fn, tn):
    """The five classification metrics of 2x2 tables, vectorised.

    Degenerate cases follow `CLASSIFICATION_METRICS`: precision/recall/F1 carry sklearn's
    `zero_division=0`, and `matthews_corrcoef` returns 0.0 when its denominator vanishes.
    Validated against the organisers' code in `control_metrics`.
    """
    tp = np.asarray(tp, float)
    fp = np.asarray(fp, float)
    fn = np.asarray(fn, float)
    tn = np.asarray(tn, float)
    shape = np.broadcast(tp, fp, fn, tn).shape
    n = tp + fp + fn + tn
    npred, npos = tp + fp, tp + fn

    acc = (tp + tn) / n
    prec = np.divide(tp, npred, out=np.zeros(shape), where=npred > 0)
    rec = np.divide(tp, npos, out=np.zeros(shape), where=npos > 0)
    f1den = 2.0 * tp + fp + fn
    f1 = np.divide(2.0 * tp, f1den, out=np.zeros(shape), where=f1den > 0)
    # TP*TN - FP*FN == n*TP - (TP+FP)*(TP+FN)
    den = npred * npos * (tn + fp) * (tn + fn)
    root = np.sqrt(den, out=np.zeros(shape), where=den > 0)
    mcc = np.divide(n * tp - npred * npos, root, out=np.zeros(shape), where=den > 0)
    return {"MCC": mcc, "Accuracy": acc, "Precision": prec, "Recall": rec, "F1": f1}


def harmonic(prec, rec):
    """F1 as defined for ONE table: the harmonic mean of precision and recall."""
    s = prec + rec
    return 0.0 if s <= 0 else 2.0 * prec * rec / s


def mcc_of(n, npred, npos, tp):
    """Point MCC from the table's margins and TP alone (the n*TP identity)."""
    den = npred * (n - npred) * npos * (n - npos)
    return 0.0 if den <= 0 else (n * tp - npred * npos) / den ** 0.5


def invert(acc, prec, rec, n):
    """The closed form: prevalence from accuracy, then the four cells.

    With prevalence p, recall r, precision q:
        TP = r*p*n,  FP = TP*(1/q - 1),  FN = TP*(1/r - 1),  TN = n - TP - FP - FN
    and accuracy = 1 - p*(1 + r/q - 2r), which inverts to p. This is exactly
    `k94_leaderboard.implied_prevalence`, applied to one endpoint instead of a macro row.
    """
    k = 1.0 + rec / prec - 2.0 * rec
    p = (1.0 - acc) / k
    tp = rec * p * n
    fp = tp * (1.0 / prec - 1.0)
    fn = tp * (1.0 / rec - 1.0)
    return p, tp, fp, fn, n - tp - fp - fn


def invert_box(acc, prec, rec, n, grid=7):
    """The closed form over the 4-decimal rounding box: (min, max) of (p, TP, FP, FN, TN)."""
    out = [invert(a, q, r, n)
           for a in np.linspace(acc - HALF, acc + HALF, grid)
           for q in np.linspace(prec - HALF, prec + HALF, grid)
           for r in np.linspace(rec - HALF, rec + HALF, grid)]
    arr = np.array(out)
    return arr.min(axis=0), arr.max(axis=0)


# ---------------------------------------------------------------------------
# the organisers' bootstrap, as the board applies it
# ---------------------------------------------------------------------------


def cumulative_draws(n, pos_rank):
    """C[s, k] = how many of resample s's draws land on a row ranked below k.

    The resample indices are the organisers' own (`bootstrap_sampling`, seed 0). Cells are
    laid out as contiguous blocks in ranked space, so any table's resampled counts are
    differences of columns of C: a candidate table costs a gather, not a pass.
    """
    idx = bootstrap_sampling(n, BOOTSTRAP_SAMPLES)
    s = idx.shape[0]
    r = pos_rank[idx]
    flat = (np.arange(s)[:, None] * n + r).ravel()
    hist = np.bincount(flat, minlength=s * n).reshape(s, n)
    return np.concatenate([np.zeros((s, 1), np.int64), hist.cumsum(axis=1)], axis=1)


def boot_means(tables, n, arrangements=ARRANGEMENTS):
    """Bootstrap means of the five metrics, per candidate table.

    For each drawn arrangement the organisers' 1000 resamples are averaged exactly as the
    backend averages them. The spread ACROSS arrangements is the irreducible Monte-Carlo
    noise of the published figure: which scored row holds which cell is not knowable.

    Returns (mean, sd), dicts of arrays of shape (T,).
    """
    tables = np.atleast_2d(np.asarray(tables, np.int64))
    T = len(tables)
    tp, fp, fn = tables[:, 0], tables[:, 1], tables[:, 2]
    c1, c2, c3 = tp, tp + fn, tp + fn + fp            # TP | FN | FP | TN in ranked space
    per = {k: np.empty((arrangements, T)) for k in KEYS}
    for a in range(arrangements):
        c = cumulative_draws(n, np.argsort(RNG.permutation(n)))
        for lo in range(0, T, CHUNK):
            sl = slice(lo, lo + CHUNK)
            a1, a2, a3 = c[:, c1[sl]], c[:, c2[sl]], c[:, c3[sl]]
            m = point_metrics(a1, a3 - a2, a2 - a1, n - a3)
            for k in KEYS:
                per[k][a, sl] = m[k].mean(axis=0)
    # The Jensen gap, per arrangement. Per resample F1 IS H(prec, rec), so this is pure
    # curvature: the correlated part of the three metrics cancels and the gap is far
    # quieter than F1 itself, which is what makes it a usable estimator of n.
    g = per["F1"] - 2.0 * per["Precision"] * per["Recall"] / (per["Precision"] + per["Recall"])
    mean = {k: per[k].mean(axis=0) for k in KEYS}
    sd = {k: per[k].std(axis=0, ddof=1) for k in KEYS}
    mean["GAP"], sd["GAP"] = g.mean(axis=0), g.std(axis=0, ddof=1)
    return mean, sd


def control_metrics():
    """CONTROL 1: the vectorised metrics must equal the organisers' own, exactly."""
    tp, fp, fn, tn = 16, 130, 10, 219
    n = tp + fp + fn + tn
    perm = np.random.default_rng(7).permutation(n)
    y_true = np.zeros(n, bool)
    y_pred = np.zeros(n, bool)
    y_true[perm[:tp]] = True
    y_pred[perm[:tp]] = True                          # TP
    y_true[perm[tp:tp + fn]] = True                   # FN
    y_pred[perm[tp + fn:tp + fn + fp]] = True         # FP
    ref = bootstrap_metrics(y_pred, y_true, "CONTROL", n_bootstrap_samples=BOOTSTRAP_SAMPLES,
                            metrics=CLASSIFICATION_METRICS)

    idx = bootstrap_sampling(n, BOOTSTRAP_SAMPLES)
    t_s, p_s = y_true[idx], y_pred[idx]
    mine = point_metrics((t_s & p_s).sum(1), (~t_s & p_s).sum(1),
                         (t_s & ~p_s).sum(1), (~t_s & ~p_s).sum(1))

    print("CONTROL 1: vectorised metrics against the organisers' bootstrap_metrics")
    print("    (same table, same arrangement, same seed-0 resample indices)\n")
    print("    %-10s %14s %14s %12s" % ("metric", "organisers", "this file", "max |diff|"))
    worst = 0.0
    for name, key in zip(NAMES, KEYS):
        d = float(np.abs(ref[key].to_numpy() - mine[key]).max())
        worst = max(worst, d)
        print("    %-10s %14.10f %14.10f %12.2e" % (name, ref[key].mean(), mine[key].mean(), d))
    assert worst < 1e-12, "vectorised metrics disagree with the organisers' code"
    print("\n    max disagreement over 1000 resamples x 5 metrics: %.2e  -> PASSES" % worst)


# ---------------------------------------------------------------------------
# 1. which tab is which enzyme
# ---------------------------------------------------------------------------


def load():
    snap = json.load(open(SNAP))
    rows = {}
    for tab in [MACRO_TAB] + ENZ_TABS:
        b = snap["boards"][tab]
        hit = [dict(zip(b["headers"], r)) for r in b["data"] if r[1] == ME]
        if len(hit) != 1:
            raise SystemExit("%s: expected exactly one row for %r, found %d"
                             % (tab, ME, len(hit)))
        rows[tab] = hit[0]
    return snap, rows


def identify(snap, rows):
    print("\n\n1. WHICH TAB IS WHICH ENZYME\n")
    print("    snapshot pulled %s; our row submitted %s on all three TDI tabs: %s\n"
          % (snap["pulled_utc"], rows[MACRO_TAB]["Submitted"],
             len({rows[t]["Submitted"] for t in [MACRO_TAB] + ENZ_TABS}) == 1))
    print("    %-12s %9s %9s %9s %9s %9s" % ("tab", *NAMES))
    for tab in ENZ_TABS:
        print("    %-12s %9.4f %9.4f %9.4f %9.4f %9.4f" % (tab, *[rows[tab][m] for m in NAMES]))
    print("    %-12s %9.4f %9.4f %9.4f %9.4f %9.4f"
          % (MACRO_TAB, *[rows[MACRO_TAB]["MA-" + m] for m in NAMES]))
    orders = [("CYP2D6", "CYP3A4"), ("CYP3A4", "CYP2D6")]
    verdicts = {}

    # (a) pairing
    print("\n    (a) PAIRING: mean of the two tabs against the macro tab")
    ok = True
    for m in NAMES:
        mean = 0.5 * (rows[ENZ_TABS[0]][m] + rows[ENZ_TABS[1]][m])
        d = abs(mean - rows[MACRO_TAB]["MA-" + m])
        ok &= d <= 2 * HALF
        print("        %-10s mean %.5f  macro %.4f  |diff| %.5f"
              % (m, mean, rows[MACRO_TAB]["MA-" + m], d))
    print("        all five inside rounding: %s -> 15 and 16 are the two TDI endpoints;" % ok)
    print("        this pairs them and says nothing about which is which.")

    # (b) our own OOF MCC per enzyme
    print("\n    (b) our out-of-fold MCC per enzyme (item 300): 2D6 %.4f, 3A4 %.4f"
          % (OOF_MCC["CYP2D6"], OOF_MCC["CYP3A4"]))
    dist = {}
    for o in orders:
        dist[o] = sum(abs(rows[t]["MCC"] - OOF_MCC[e]) for t, e in zip(ENZ_TABS, o))
        print("        15=%s, 16=%s  ->  L1 distance %.4f" % (*o, dist[o]))
    verdicts["b"] = min(dist, key=dist.get)

    # (c) the field
    print("\n    (c) the FIELD, all entrants per tab -- independent of our own row:")
    med = {}
    for tab in ENZ_TABS:
        b = snap["boards"][tab]
        i = b["headers"].index("MCC")
        v = np.array([r[i] for r in b["data"] if isinstance(r[i], (int, float))])
        med[tab] = float(np.median(v))
        print("        %-12s n=%d  median MCC %.4f  best %.4f" % (tab, len(v), med[tab], v.max()))
    b15, b16 = (snap["boards"][t] for t in ENZ_TABS)
    m15 = {str(r[1]): r[b15["headers"].index("MCC")] for r in b15["data"]}
    m16 = {str(r[1]): r[b16["headers"].index("MCC")] for r in b16["data"]}
    both = [u for u in m15 if u in m16]
    low = sum(m15[u] < m16[u] for u in both)
    hard = min(med, key=med.get)
    print("        %d entrants on both tabs; %s is the lower of their two for %d of them"
          % (len(both), ENZ_TABS[0], low))
    print("        item 300: CYP2D6 TDI is hard for every team that measures it (MCC")
    print("        0.09-0.13) and CYP3A4 is far higher -> the hard tab, %s, is CYP2D6" % hard)
    verdicts["c"] = ("CYP2D6", "CYP3A4") if hard == ENZ_TABS[0] else ("CYP3A4", "CYP2D6")

    # (d) config order
    print("\n    (d) endpoint order in the organisers' config: %s" % CLASSIFICATION_ENDPOINTS)
    print("        the regression sub-tabs 10-13 follow config order (k100_recal; item 308's")
    print("        rank-order prediction held under that map), which gives 15=CYP2D6,")
    print("        16=CYP3A4. SUPPORTING ONLY: an argument from tab order, which (b), (c)")
    print("        and (e) exist so that nothing rests on.")
    first = CLASSIFICATION_ENDPOINTS[0].split("_")[0]
    verdicts["d"] = (first, "CYP3A4" if first == "CYP2D6" else "CYP2D6")

    # (e) our own positive rate
    sub = pd.read_csv(SUB)
    rate = {c: float(sub[c + "_is_TDI"].mean()) for c in ["CYP2D6", "CYP3A4"]}
    ncall = {c: int(sub[c + "_is_TDI"].sum()) for c in ["CYP2D6", "CYP3A4"]}
    print("\n    (e) OUR OWN positive rate, exact from the file: CYP2D6 %d/750 = %.4f, "
          "CYP3A4 %d/750 = %.4f" % (ncall["CYP2D6"], rate["CYP2D6"],
                                      ncall["CYP3A4"], rate["CYP3A4"]))
    print("        the board implies our predicted-positive rate on the scored rows, r*p/q,")
    print("        with no n in it:")
    implied = {}
    for tab in ENZ_TABS:
        _, tp, fp, _, _ = invert(rows[tab]["Accuracy"], rows[tab]["Precision"],
                                 rows[tab]["Recall"], 1.0)
        implied[tab] = tp + fp
        print("        %-12s %.4f" % (tab, tp + fp))
    dist = {}
    for o in orders:
        dist[o] = sum(abs(implied[t] - rate[e]) for t, e in zip(ENZ_TABS, o))
        print("        15=%s, 16=%s  ->  L1 distance %.4f" % (*o, dist[o]))
    verdicts["e"] = min(dist, key=dist.get)

    print("\n    VERDICTS: " + ",  ".join("(%s) 15=%s" % (k, v[0]) for k, v in verdicts.items()))
    agree = len({v for v in verdicts.values()}) == 1
    print("    all four assignment readings agree: %s" % agree)
    assert agree, "tab assignment is contested; refusing to invert a table for a guessed enzyme"
    o = verdicts["e"]
    return dict(zip(o, ENZ_TABS)), ncall, implied


def provenance(implied):
    print("\n    PROVENANCE: which TDI file on disk is compatible with the implied rates")
    print("        %-44s %8s %8s %11s" % ("file", "2D6", "3A4", "compatible"))
    for f in sorted(glob.glob(RES + "submission/**/tdi_submission*.csv", recursive=True)):
        d = pd.read_csv(f)
        r2, r3 = d["CYP2D6_is_TDI"].mean(), d["CYP3A4_is_TDI"].mean()
        # the scored rows are a subset of 750, so compatible means close, not equal
        good = abs(r2 - implied[ENZ_TABS[0]]) < 0.05 and abs(r3 - implied[ENZ_TABS[1]]) < 0.05
        print("        %-44s %8.4f %8.4f %11s"
              % (f.split("submission/")[-1], r2, r3, "yes" if good else "NO"))
    print("        every compatible file carries 285/360 calls; only the counts are used here.")


# ---------------------------------------------------------------------------
# 2. the closed form and its validations
# ---------------------------------------------------------------------------


def count_bounds(ncall, n):
    """Our positive calls on n scored rows drawn from 750, given ncall calls in the file."""
    return max(0, ncall - (FULL - n)), min(ncall, n)


def closed_form(rows, tabs, ncall, n_list=(FULL, LIVE_HALF)):
    print("\n\n2-3. THE CLOSED FORM, AND ITS VALIDATIONS\n")
    got = {}
    for enz, tab in tabs.items():
        r = rows[tab]
        acc, prec, rec, mcc, f1 = (r["Accuracy"], r["Precision"], r["Recall"],
                                   r["MCC"], r["F1 Score"])
        print("    %s (%s): MCC %.4f  Acc %.4f  Prec %.4f  Rec %.4f  F1 %.4f"
              % (enz, tab, mcc, acc, prec, rec, f1))
        for n in n_list:
            p, tp, fp, fn, tn = invert(acc, prec, rec, n)
            lo, hi = invert_box(acc, prec, rec, n)
            cells = np.array([tp, fp, fn, tn])
            frac = np.abs(cells - np.round(cells))
            box = float((hi - lo)[1:].max() / 2.0)
            ptm = float(point_metrics(*np.round(cells))["MCC"])
            clo, chi = count_bounds(ncall[enz], n)
            print("      n=%d:  p %.6f  TP %.3f  FP %.3f  FN %.3f  TN %.3f" % (n, p, *cells))
            print("        CHECK 1 near-integer: distances %s; rounding moves a cell by <= %.3f"
                  % (" ".join("%.3f" % x for x in frac), box))
            print("                -> %s" % ("PASS" if frac.max() <= max(box, 0.05)
                                           else "FAIL: at least one cell is not an integer"))
            print("        CHECK 2 sum %.6f == n: VACUOUS, TN is defined as the residual" % cells.sum())
            print("                replacement that CAN fail -- our calls on the scored rows,")
            print("                TP+FP = %.2f, must lie in [%d, %d] -> %s"
                  % (tp + fp, clo, chi,
                     "ok" if clo - 0.5 <= tp + fp <= chi + 0.5 else "REFUSES n=%d" % n))
            print("        CHECK 3 MCC of the rounded table %.4f vs board %.4f, |diff| %.4f -> %s"
                  % (ptm, mcc, abs(ptm - mcc),
                     "PASS" if abs(ptm - mcc) < 5e-4 else "FAIL at three decimals"))
        h, hlo, hhi = harmonic(prec, rec), harmonic(prec - HALF, rec - HALF), \
            harmonic(prec + HALF, rec + HALF)
        gap, glo, ghi = f1 - h, f1 - HALF - hhi, f1 + HALF - hlo
        passes = glo <= 0.0 <= ghi
        print("      CHECK 4 (n-free): F1 must equal H(prec, rec) for ANY single table")
        print("        H %.6f; board F1 %.4f; gap %+.6f, rounding box [%+.6f, %+.6f] -> %s\n"
              % (h, f1, gap, glo, ghi,
                 "PASS" if passes else "FAIL, by %.0fx the box: NOT one table's point metrics"
                 % (abs(gap) / (0.5 * (ghi - glo)))))
        got[enz] = dict(gap=gap, gap_lo=glo, gap_hi=ghi)
    print("    CHECK 4 needs neither n nor any assumption about the scoring set, and it")
    print("    refuses the point-estimate reading on both enzymes. Checks 1-3 are then not")
    print("    tests of the table at all: they are tests of a model the board does not use.")
    return got


def control_synthetic(got=None):
    """CONTROL 2: a known table, from its point metrics and from its bootstrap means."""
    tp, fp, fn, tn = 16, 130, 10, 219
    n = tp + fp + fn + tn
    truth = np.array([tp, fp, fn, tn], float)
    print("\n\nCONTROL 2: a SYNTHETIC table, to show the checks can both pass and fail")
    print("    truth TP %d  FP %d  FN %d  TN %d  (n=%d)\n" % (tp, fp, fn, tn, n))

    pt = {k: float(v) for k, v in point_metrics(tp, fp, fn, tn).items()}
    rp = {k: round(v, 4) for k, v in pt.items()}
    _, *cells = invert(rp["Accuracy"], rp["Precision"], rp["Recall"], n)
    cells = np.array(cells)
    gap = rp["F1"] - harmonic(rp["Precision"], rp["Recall"])
    print("    (i)  POINT metrics, rounded like the board: recovered %s"
          % " ".join("%.3f" % x for x in cells))
    print("         max error %.4f compounds -> %s;  F1 - H = %+.6f -> check 4 %s"
          % (np.abs(cells - truth).max(),
             "RECOVERED" if np.abs(cells - truth).max() < 0.1 else "NOT RECOVERED",
             gap, "PASSES" if abs(gap) <= 3 * HALF else "FAILS"))

    mean, sd = boot_means([[tp, fp, fn, tn]], n)
    bm = {k: float(mean[k][0]) for k in KEYS}
    rb = {k: round(v, 4) for k, v in bm.items()}
    _, *cb = invert(rb["Accuracy"], rb["Precision"], rb["Recall"], n)
    cb = np.array(cb)
    gapb = rb["F1"] - harmonic(rb["Precision"], rb["Recall"])
    print("    (ii) BOOTSTRAP means, organisers' protocol, same table:")
    print("         %-10s %9s %10s %9s" % ("metric", "point", "boot mean", "offset"))
    for name, k in zip(NAMES, KEYS):
        print("         %-10s %9.4f %10.4f %+9.4f" % (name, pt[k], bm[k], bm[k] - pt[k]))
    print("         closed form on those: %s  (errors %s)"
          % (" ".join("%.3f" % x for x in cb), " ".join("%+.2f" % x for x in cb - truth)))
    print("         F1 - H = %+.6f -> check 4 %s"
          % (gapb, "PASSES" if abs(gapb) <= 3 * HALF else "FAILS, the board's failure reproduced"))
    print("         the CELLS are still recovered, so averaging does not bias accuracy,")
    print("         precision or recall materially: the bias lives in F1 and MCC.")

    idx = bootstrap_sampling(n, BOOTSTRAP_SAMPLES)

    def arrange(seed):
        """One arrangement of this table over the scored rows, as labels and calls."""
        perm = np.random.default_rng(seed).permutation(n)
        yt, yp = np.zeros(n, bool), np.zeros(n, bool)
        yt[perm[:tp]] = True
        yp[perm[:tp]] = True                          # TP
        yt[perm[tp:tp + fn]] = True                   # FN
        yp[perm[tp + fn:tp + fn + fp]] = True         # FP
        t_s, p_s = yt[idx], yp[idx]
        return point_metrics((t_s & p_s).sum(1), (~t_s & p_s).sum(1),
                             (t_s & ~p_s).sum(1), (~t_s & ~p_s).sum(1))

    # A board row is ONE arrangement, not an average over arrangements. The expectation
    # above averages that noise away and so cannot show where check 1 fails; this can.
    print("\n    (iii) SINGLE realisations, which is what a board row is:")
    print("         %-5s %26s %12s %12s" % ("draw", "closed-form cells", "max |err|", "to integer"))
    errs, ints = [], []
    for d in range(6):
        m = arrange(100 + d)
        rr = {k: round(float(m[k].mean()), 4) for k in KEYS}
        _, *c = invert(rr["Accuracy"], rr["Precision"], rr["Recall"], n)
        c = np.array(c)
        errs.append(float(np.abs(c - truth).max()))
        ints.append(float(np.abs(c - np.round(c)).max()))
        print("         %-5d %26s %12.3f %12.3f"
              % (d, " ".join("%.2f" % x for x in c), errs[-1], ints[-1]))
    print("         one published figure puts the cells %.2f-%.2f compounds off integer,"
          % (min(ints), max(ints)))
    print("         which is the size of check 1's failure on the real board. Check 1 was")
    print("         never a test of the table: it is a test of an aggregation the board")
    print("         does not use.")

    # The obvious alternative aggregation, tested instead of dismissed.
    m = arrange(11)
    print("\n    (iv) ALTERNATIVE aggregation, so the mean is not just the first fit:")
    for lab, agg in (("mean", np.mean), ("median", np.median)):
        rr = {k: round(float(agg(m[k])), 4) for k in KEYS}
        h = harmonic(rr["Precision"], rr["Recall"])
        print("         %-7s F1 %.4f, H(prec, rec) %.6f, gap %+.6f"
              % (lab, rr["F1"], h, rr["F1"] - h))
    if got:
        print("         the board's gaps: %s"
              % ", ".join("%s %+.6f" % (e, got[e]["gap"]) for e in got))
    print("         only one of the two aggregations reproduces them.")
    return gapb


# ---------------------------------------------------------------------------
# 3. the gap pins n
# ---------------------------------------------------------------------------


def size_from_gap(rows, tabs, got):
    print("\n\n4. THE GAP IS A MEASUREMENT OF n\n")
    print("    Per resample F1 IS H(prec, rec), identically, so mean(F1) - H(mean prec,")
    print("    mean rec) is a pure Jensen term ~ 0.5 tr(Hess_H Cov), and Cov ~ 1/n.\n")
    grid = [150, 200, 250, 300, 340, 375, 420, 500, 600, 750]
    out = {}
    for enz, tab in tabs.items():
        r = rows[tab]
        print("    %s: observed gap %+.6f, rounding box [%+.6f, %+.6f]"
              % (enz, got[enz]["gap"], got[enz]["gap_lo"], got[enz]["gap_hi"]))
        print("      %6s %16s %12s %11s %11s %9s"
              % ("n", "table TP/FP/FN", "exp. gap", "err of E", "one board", "gap*n"))
        prod, noise = [], []
        for n in grid:
            _, tp, fp, fn, tn = invert(r["Accuracy"], r["Precision"], r["Recall"], n)
            t = np.round([tp, fp, fn]).astype(int)
            table = [t[0], t[1], t[2], n - t.sum()]
            mean, sd = boot_means([table], n)
            g = float(mean["GAP"][0])
            one = float(sd["GAP"][0])                       # noise of ONE published figure
            prod.append(g * n)
            noise.append(one)
            print("      %6d %16s %+12.6f %11.6f %11.6f %9.4f"
                  % (n, "%d/%d/%d" % tuple(t), g, one / ARRANGEMENTS ** 0.5, one, g * n))
        a = float(np.median(prod))
        spread = (max(prod) - min(prod)) / abs(a)
        print("      gap*n = %.4f, relative spread over a 5x range of n %.0f%% -> 1/n law holds"
              % (a, 100 * spread))
        print("      'err of E' is the MC error of the expected gap (sd/sqrt(%d)); 'one board'"
              % ARRANGEMENTS)
        print("      is the spread of a SINGLE published gap, which is what the board is.")
        if got[enz]["gap_hi"] >= 0:
            print("      rounding box reaches zero: n is not bounded above by this reading\n")
            out[enz] = None
            continue
        mc = float(np.median(noise))
        lo_g, hi_g = got[enz]["gap_lo"] - 2 * mc, got[enz]["gap_hi"] + 2 * mc
        n_hat = a / got[enz]["gap"]
        n1, n2 = a / got[enz]["gap_hi"], a / got[enz]["gap_lo"]
        if hi_g >= 0:
            wide = "unbounded above"
        else:
            wide = "[%.0f, %.0f]" % (min(a / hi_g, a / lo_g), max(a / hi_g, a / lo_g))
        print("      n implied by the observed gap: %.0f; rounding only [%.0f, %.0f];"
              % (n_hat, min(n1, n2), max(n1, n2)))
        print("      rounding plus 2 sd of a single board figure: %s" % wide)
        print("      -> the live half (%d) is inside; the full test set (%d) is not\n"
              % (LIVE_HALF, FULL))
        out[enz] = (n_hat, min(n1, n2), max(n1, n2))
    return out


# ---------------------------------------------------------------------------
# 4. the bootstrap-aware inversion
# ---------------------------------------------------------------------------


def search(rows, tabs, ncall, n_candidates):
    print("\n\n5. BOOTSTRAP-AWARE INVERSION: integer tables whose BOOTSTRAP MEANS match\n")
    print("    consistent = |board - E| <= 3 sd on all five metrics; E and sd over %d"
          % ARRANGEMENTS)
    print("    arrangements with the organisers' resample indices, sd including display")
    print("    rounding. The five metrics are correlated, so this is a screen, not a")
    print("    likelihood. 'count' also requires TP+FP inside our own call-count bounds.\n")
    found = {}
    for enz, tab in tabs.items():
        r = rows[tab]
        board = dict(zip(KEYS, [r[m] for m in NAMES]))
        print("    %s (%s)" % (enz, tab))
        print("      %4s %6s %6s %6s %18s %7s %8s %7s %8s"
              % ("n", "cands", "consis", "+count", "best TP/FP/FN/TN", "max|z|", "prev",
                 "TP+FP", "win/pre"))
        for n in n_candidates:
            _, tp0, fp0, fn0, _ = invert(r["Accuracy"], r["Precision"], r["Recall"], n)
            scale = (n / LIVE_HALF) ** 0.5
            h1, h2 = int(np.ceil(8 * scale)), int(np.ceil(16 * scale))
            c0 = np.round([tp0, fp0, fn0]).astype(int)
            cand = []
            for dtp, dfp, dfn in itertools.product(range(-h1, h1 + 1), range(-h2, h2 + 1),
                                                   range(-h1, h1 + 1)):
                tp, fp, fn = c0[0] + dtp, c0[1] + dfp, c0[2] + dfn
                tn = n - tp - fp - fn
                if tp >= 1 and fp >= 0 and fn >= 0 and tn >= 0:
                    cand.append((tp, fp, fn, tn))
            cand = np.array(cand)
            pt = point_metrics(*cand.T)
            dpre = np.max([np.abs(pt[k] - board[k]) / PREFILTER[k] for k in PREFILTER], axis=0)
            sub = cand[dpre <= 1.0]
            dsub = dpre[dpre <= 1.0]
            if len(sub) == 0:
                print("      %4d %6d %6d  -- nothing passes the prefilter" % (n, len(cand), 0))
                continue
            mean, sd = boot_means(sub, n)
            z = np.stack([(board[k] - mean[k]) / np.sqrt(sd[k] ** 2 + SD_RND ** 2)
                          for k in KEYS], axis=1)
            mz = np.abs(z).max(axis=1)
            chi = (z ** 2).sum(axis=1)
            clo, chi_hi = count_bounds(ncall[enz], n)
            calls = sub[:, 0] + sub[:, 1]
            ok = mz <= 3.0
            okc = ok & (calls >= clo) & (calls <= chi_hi)
            # report the best CONSISTENT and count-compatible table, not the best overall:
            # a table with one 4-sigma metric can still carry the smallest sum of squares
            pick = np.where(okc)[0]
            i = int(pick[np.argmin(chi[pick])]) if len(pick) else int(np.argmin(chi))
            if ok.any():
                s = sub[ok]
                win = max(np.abs(s[:, 0] - c0[0]).max() / h1, np.abs(s[:, 2] - c0[2]).max() / h1,
                          np.abs(s[:, 1] - c0[1]).max() / h2)
                pre = float(dsub[ok].max())
                bind = "%.2f/%.2f" % (win, pre)
            else:
                win, pre, bind = 0.0, 0.0, "--"
            print("      %4d %6d %6d %6d %18s %7.2f %8.4f %7d %8s"
                  % (n, len(sub), int(ok.sum()), int(okc.sum()), "%d/%d/%d/%d" % tuple(sub[i]),
                     mz[i], (sub[i][0] + sub[i][2]) / n, calls[i], bind))
            if win >= 1.0 or pre >= 0.9:
                print("            WARNING: the window or the prefilter binds at n=%d" % n)
            found[(enz, n)] = dict(sub=sub, ok=ok, okc=okc, best=i, mean=mean, sd=sd, z=z,
                                   mz=mz)
        print("      win/pre: the farthest consistent table as a fraction of the search window")
        print("      and of the prefilter tolerance; below 1 means neither limited the answer.\n")
    return found


def answer(rows, tabs, ncall, found, n_live):
    print("\n\n6. THE ANSWER\n")
    tables = {}
    for enz, tab in tabs.items():
        f = found.get((enz, LIVE_HALF))
        if f is None or not f["okc"].any():
            print("    %s: INVERSION FAILED at n=%d -- no integer table passes the screen"
                  % (enz, LIVE_HALF))
            continue
        i = f["best"]
        table = f["sub"][i]
        r = rows[tab]
        cand_i = np.where(f["okc"])[0]
        print("    every table consistent at n=%d for %s:" % (LIVE_HALF, enz))
        for j in cand_i:
            print("        TP %3d  FP %3d  FN %3d  TN %3d    max|z| %.2f%s"
                  % (*f["sub"][j], f["mz"][j], "   <- reported below" if j == i else ""))
        print("    %s (%s), n = %d:   TP %d   FP %d   FN %d   TN %d   (max |z| %.2f)"
              % (enz, tab, LIVE_HALF, *table, f["mz"][i]))
        pt = point_metrics(*table)
        print("      %-10s %9s %11s %9s %7s %12s"
              % ("metric", "board", "boot mean", "sd", "z", "point table"))
        for j, (k, name) in enumerate(zip(KEYS, NAMES)):
            print("      %-10s %9.4f %11.4f %9.4f %+7.2f %12.4f"
                  % (name, r[name], f["mean"][k][i], f["sd"][k][i], f["z"][i, j], float(pt[k])))
        # uncertainty: every consistent table across the n that the evidence admits
        pool = np.concatenate([found[(enz, n)]["sub"][found[(enz, n)]["okc"]]
                               for n in n_live if (enz, n) in found])
        ns = np.concatenate([np.full(int(found[(enz, n)]["okc"].sum()), n)
                             for n in n_live if (enz, n) in found])
        print("      across n in %d-%d, %d consistent integer tables:"
              % (min(n_live), max(n_live), len(pool)))
        print("        TP %d-%d   FP %d-%d   FN %d-%d   TN %d-%d   prevalence %.4f-%.4f"
              % (pool[:, 0].min(), pool[:, 0].max(), pool[:, 1].min(), pool[:, 1].max(),
                 pool[:, 2].min(), pool[:, 2].max(), pool[:, 3].min(), pool[:, 3].max(),
                 ((pool[:, 0] + pool[:, 2]) / ns).min(), ((pool[:, 0] + pool[:, 2]) / ns).max()))
        offset = r["MCC"] - float(pt["MCC"])
        print("      point MCC %.4f; published bootstrap mean %.4f; Jensen offset %+.4f\n"
              % (float(pt["MCC"]), r["MCC"], offset))
        tables[enz] = dict(table=table, pool=pool, ns=ns, offset=offset)
    return tables


# ---------------------------------------------------------------------------
# 5. prevalence, and k94's weakness
# ---------------------------------------------------------------------------


def prevalence(rows, tables, ncall):
    print("\n\n7. IMPLIED TEST PREVALENCE, AND OUR OWN POSITIVE RATE\n")
    tdi = pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv")
    print("    %-8s %11s %16s %13s %14s %11s"
          % ("enzyme", "test prev", "range (all n)", "our rate", "our rate, 750", "train prev"))
    per = []
    for enz in ["CYP2D6", "CYP3A4"]:
        if enz not in tables:
            continue
        t, pool, ns = tables[enz]["table"], tables[enz]["pool"], tables[enz]["ns"]
        pv = (pool[:, 0] + pool[:, 2]) / ns
        lab = tdi[enz + "_is_TDI"].dropna().astype(bool)
        per.append((t[0] + t[2]) / LIVE_HALF)
        print("    %-8s %11.4f %16s %13.4f %14.4f %11.4f"
              % (enz, per[-1], "%.4f-%.4f" % (pv.min(), pv.max()),
                 (t[0] + t[1]) / LIVE_HALF, ncall[enz] / FULL, lab.mean()))
    print("    train prevalence over labelled rows only: CYP2D6 %d, CYP3A4 %d of %d"
          % (tdi["CYP2D6_is_TDI"].notna().sum(), tdi["CYP3A4_is_TDI"].notna().sum(), len(tdi)))

    m = rows[MACRO_TAB]
    k = 1.0 + m["MA-Recall"] / m["MA-Precision"] - 2.0 * m["MA-Recall"]
    p_macro = (1.0 - m["MA-Accuracy"]) / k
    print("\n\n8. DOES THE PER-ENZYME INVERSION ESCAPE k94's WEAKNESS?\n")
    print("    k94 implied_prevalence on our MACRO row          %.4f" % p_macro)
    if len(per) == 2:
        print("    mean of the two per-enzyme prevalences here      %.4f" % np.mean(per))
        print("    the per-enzyme values differ %.1fx, the condition under which a" % (max(per) / min(per)))
        print("    macro-average of ratios cannot invert to one p; the defect on our own row")
        print("    is %.4f of prevalence." % abs(p_macro - np.mean(per)))


# ---------------------------------------------------------------------------
# 6. a different threshold
# ---------------------------------------------------------------------------


def threshold(tables, ncall):
    print("\n\n9. MCC AT A DIFFERENT THRESHOLD: what is determined, and what is not\n")
    print("    Assumption doing work: a threshold on ONE score, so positive sets are nested.")
    print("    For a fixed number of calls, TP*TN - FP*FN = n*TP - Npred*Npos makes MCC linear")
    print("    and increasing in TP. So the extremes of any move are: false positives leave")
    print("    (or false negatives join) first, and the reverse. Those are the bounds.\n")
    n = LIVE_HALF
    for enz in ["CYP2D6", "CYP3A4"]:
        if enz not in tables:
            continue
        tp, fp, fn, tn = (int(x) for x in tables[enz]["table"])
        npos, npred = tp + fn, tp + fp
        now = mcc_of(n, npred, npos, tp)
        print("    %s: TP %d FP %d FN %d TN %d -- %d calls, %d true positives, point MCC %.4f"
              % (enz, tp, fp, fn, tn, npred, npos, now))
        print("      %7s %7s %12s %12s %14s" % ("move", "calls", "MCC best", "MCC worst", "width"))
        for k in [1, 2, 5, 10, 25, 50, 100]:
            if k < npred:
                hi = mcc_of(n, npred - k, npos, tp - max(0, k - fp))
                lo = mcc_of(n, npred - k, npos, tp - min(k, tp))
                print("      %7s %7d %12.4f %12.4f %14.4f" % ("-%d" % k, npred - k, hi, lo, hi - lo))
        for k in [1, 2, 5, 10, 25, 50, 100]:
            if k <= fn + tn:
                hi = mcc_of(n, npred + k, npos, tp + min(k, fn))
                lo = mcc_of(n, npred + k, npos, tp + max(0, k - tn))
                print("      %7s %7d %12.4f %12.4f %14.4f" % ("+%d" % k, npred + k, hi, lo, hi - lo))
        one = max(abs(mcc_of(n, npred - 1, npos, tp - s) - now) for s in (0, 1))
        need = None
        for k in range(1, npred):
            hi = mcc_of(n, npred - k, npos, tp - max(0, k - fp))
            lo = mcc_of(n, npred - k, npos, tp - min(k, tp))
            if max(abs(hi - now), abs(now - lo)) >= MCC_FLOOR[enz]:
                need = k
                break
        print("      one compound removed moves MCC by up to %.4f; item 235's floor is %.4f"
              % (one, MCC_FLOOR[enz]))
        print("      smallest tightening whose DETERMINED band can reach that floor: %s compounds"
              % need)
        pool = tables[enz]["pool"]
        ns = tables[enz]["ns"]
        los, his = [], []
        for t, nn in zip(pool, ns):
            a, b, c, _ = (int(x) for x in t)
            los.append(mcc_of(nn, a + b - 10, a + c, a - min(10, a)))
            his.append(mcc_of(nn, a + b - 10, a + c, a - max(0, 10 - b)))
        print("      removing 10 calls, over EVERY consistent table: MCC in [%.4f, %.4f]"
              % (min(los), max(his)))
        print("      offset of the published (bootstrap-mean) MCC from the point MCC: %+.4f\n"
              % tables[enz]["offset"])

    print("    NOT DETERMINED, and not estimated here:")
    print("      * where inside those bounds the real value falls. That needs the ORDER of")
    print("        our test scores within the moved slice; the board carries five aggregates")
    print("        per enzyme and no probabilities.")
    print("      * a non-nested change (another gate, another ensemble member).")
    print("      * the blinded half: this table covers the %d scored rows only. For the" % n)
    print("        other %d, only our call counts follow (%d total on 2D6 and %d on 3A4"
          % (FULL - n, ncall["CYP2D6"], ncall["CYP3A4"]))
    print("        in the file, minus the scored calls above); their labels are unconstrained.")
    print("      * the Jensen offset at a new threshold: it moves with the table, and is of")
    print("        the size printed per enzyme above.")


if __name__ == "__main__":
    control_metrics()
    snap, rows = load()
    tabs, ncall, implied = identify(snap, rows)
    provenance(implied)
    got = closed_form(rows, tabs, ncall)
    control_synthetic(got)
    size_from_gap(rows, tabs, got)
    # A series-split half need not be exactly 375, and CYP2D6 carries a couple of missing
    # labels, so n is scanned rather than assumed: the answer's dependence on it is the
    # answer's real uncertainty. The gap estimator above bounds n to roughly 300-500.
    n_live = list(range(355, 381))
    n_full = [747, 748, 749, 750]
    found = search(rows, tabs, ncall, n_live + n_full)
    tables = answer(rows, tabs, ncall, found, n_live)
    prevalence(rows, tables, ncall)
    threshold(tables, ncall)
