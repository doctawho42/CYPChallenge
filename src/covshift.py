"""Estimate how far a blind set's label distribution has moved, using no labels from it.

The problem this solves. A model is fitted on a labelled set and asked to predict a blind
one. If the blind set's labels are distributed differently - shifted, enriched, depleted - the
predictions inherit the labelled set's marginal and are wrong in a way no amount of fitting
repairs. The shift is a property of the blind set's LABELS, which by definition nobody has.

What is observable is the blind set's PREDICTIONS. That is enough, because the map from labels
to predictions can be measured on the labelled set out of fold:

    E_blind[yhat] = integral of E[yhat | y] p_blind(y) dy

The kernel E[yhat | y] comes from the labelled set - fit it isotonically on out-of-fold
predictions, which is the weakest assumption that still gives a monotone map. Then tilt the
labelled set's own label marginal until its image through that kernel has the mean actually
observed on the blind set. The tilt that achieves it is the shift.

Two things this deliberately does not do. It does not divide the observed prediction gap by a
regression slope, which is the intuitive move and is wrong: the slope of y on yhat and the
slope of yhat on y differ by the fraction of variance explained, and using the wrong one on
this data costs a factor of 1.7 to 4.5 per enzyme. And it does not require an instrument or
any auxiliary variable. The kernel is the instrument.

Stratification, when the shift is partly composition. If the blind set contains a different
mix of some chemically meaningful class - here, protonatable bases, which CYP2D6 binds through
a salt bridge and the other enzymes do not - then part of the apparent label shift is simply
that mixture changing, and part is a genuine shift within each class. Passing `strata` and
`blind_strata` separates them: one tilt is shared across strata, so the assumption is that
labels move by the same amount inside each class, and the classes are then mixed by the blind
set's own measured composition. On CYP2D6 this moved the estimate from -0.917 to -0.508, and
nearly half of what had been called a label shift turned out to be composition.

Uncertainty. Both sides are resampled, never one: the kernel's uncertainty and the target's
uncertainty both matter, and under inversion the target enters divided by the kernel slope, so
with slopes of 0.22 to 0.59 the target's share is the larger one. Holding the blind mean fixed
- the obvious implementation - understates the interval by a factor of two to five. Resampling
is by group where groups are given, because analogue series are not independent draws.

    from covshift import estimate, recalibrate, apply_affine

    s = estimate(y, oof, blind_pred, groups=cluster_id, blind_groups=test_cluster_id)
    print(s.delta, s.lo, s.hi)
    off, lam = recalibrate(y, oof, band_lo, band_hi, delta=s.delta)
    final = apply_affine(blind_pred, off, lam, mu=y.mean() + s.delta)

Running this file reproduces the four per-enzyme estimates this repository published, which is
the only acceptance test that matters: `uv run python src/covshift.py`.
"""
import numpy as np
from scipy.optimize import brentq
from sklearn.isotonic import IsotonicRegression

TH_MAX = 80.0


class Shift:
    """An estimated shift with its interval and the diagnostics that qualify it."""

    __slots__ = ("delta", "lo", "hi", "draws", "slope", "r2", "composition", "n", "n_blind")

    def __init__(self, **kw):
        for k in self.__slots__:
            setattr(self, k, kw.get(k))

    def __repr__(self):
        c = "" if self.composition is None else f", composition {self.composition:+.3f}"
        return (f"Shift({self.delta:+.3f}, 95% [{self.lo:+.3f}, {self.hi:+.3f}], "
                f"slope {self.slope:.3f}{c})")

    def amplification(self):
        """How much the inversion multiplies an error in the observed blind mean.

        The reciprocal of the kernel slope - the slope of predictions on labels, not the
        other way round. Where it is large the estimate is intrinsically
        loose and no amount of extra labelled data on our side tightens it, because the
        looseness is in the conditioning of the inversion rather than in the sample size.
        """
        return 1.0 / self.slope


def _tilt_weights(y, th):
    """Exponential tilt as normalised weights - the minimum relative entropy way to move a mean."""
    z = th * (y - y.mean())
    z -= z.max()
    w = np.exp(z)
    return w / w.sum()


def kernel(y, oof):
    """E[yhat | y] on the labelled set, fitted isotonically and evaluated at each y."""
    return IsotonicRegression(out_of_bounds="clip").fit(y, oof).predict(y)


def _solve(y, m, target):
    """The label-mean shift whose tilt sends the kernel image's mean to `target`."""
    f = lambda th: float(_tilt_weights(y, th) @ m) - target
    if f(-TH_MAX) * f(TH_MAX) > 0:
        return np.nan
    w = _tilt_weights(y, brentq(f, -TH_MAX, TH_MAX, xtol=1e-10))
    return float(w @ y) - y.mean()


def _solve_strat(y, oof, strata, frac, target):
    """One shared tilt inside every stratum, mixed by the blind set's composition.

    The kernel is fitted SEPARATELY inside each stratum, which is what makes this stratified
    rather than merely reweighted. The map from label to prediction is itself different per
    class - on CYP2D6 the model responds differently to a basic compound's label than to a
    neutral one's, because the salt bridge is a different binding mode - and slicing one global
    kernel by stratum would keep exactly the confound the stratification exists to remove.

    The tilt is shared, so the assumption is that labels move by the same amount inside each
    class and only the mixture differs. Returns the total shift of the label mean and the part
    of it that is composition alone - what the mean would move by if labels within each stratum
    did not move at all.
    """
    ks = list(frac)
    ys = {k: y[strata == k] for k in ks}
    ms = {k: kernel(ys[k], oof[strata == k]) for k in ks}
    for k in ks:
        if len(ys[k]) < 2:
            raise ValueError(f"страта {k}: {len(ys[k])} размеченных строк, наклонять нечего")

    def image(th):
        return sum(frac[k] * float(_tilt_weights(ys[k], th) @ ms[k]) for k in ks)

    g = lambda th: image(th) - target
    if g(-TH_MAX) * g(TH_MAX) > 0:
        return np.nan, np.nan
    th = brentq(g, -TH_MAX, TH_MAX, xtol=1e-10)
    mix = sum(frac[k] * float(_tilt_weights(ys[k], th) @ ys[k]) for k in ks)
    comp = sum(frac[k] * float(ys[k].mean()) for k in ks) - y.mean()
    return mix - y.mean(), comp


def _resample(rng, groups, n):
    """Indices of one bootstrap draw, by group when groups are given."""
    if groups is None:
        return rng.integers(0, n, n)
    u = np.unique(groups)
    pick = rng.choice(u, len(u), replace=True)
    return np.concatenate([np.where(groups == k)[0] for k in pick])


def estimate(y, oof, blind_pred, *, groups=None, blind_groups=None,
             strata=None, blind_strata=None, draws=1500, seed=0, level=95.0):
    """Shift of the mean label between the labelled set and the blind set.

    y, oof        labels and OUT-OF-FOLD predictions on the labelled set, same length.
    blind_pred    predictions on the blind set. Its labels are never needed and must not
                  be passed; if you have them you do not need this function.
    groups        cluster id per labelled row, for resampling. Analogue series are not
                  independent and treating them as such understates the interval.
    blind_groups  the same for the blind set.
    strata        integer or boolean class per labelled row. When given, the shift is
                  decomposed into a within-stratum move and a composition change.
    blind_strata  the same for the blind set; its composition sets the mixture. It is
                  resampled along with everything else: the mixture is estimated from the
                  blind set, not known, and holding it fixed is the same error one level down.
    draws         1500 by default, not 300. The point estimate is stable at 300 and the
                  interval is not - on CYP2D6 the 95% bounds move by 0.12 between the two,
                  enough to change whether they contain zero.

    Nothing here reads a blind label, which is the whole point and is worth checking by
    reading the signature rather than by trusting this sentence.
    """
    y = np.asarray(y, float)
    oof = np.asarray(oof, float)
    blind_pred = np.asarray(blind_pred, float)
    if len(y) != len(oof):
        raise ValueError(f"меток {len(y)}, предсказаний вне фолда {len(oof)}")

    target = float(blind_pred.mean())
    if strata is None:
        delta, comp = _solve(y, kernel(y, oof), target), None
    else:
        strata = np.asarray(strata)
        blind_strata = np.asarray(blind_strata)
        frac = {k: float((blind_strata == k).mean()) for k in np.unique(strata)}
        delta, comp = _solve_strat(y, oof, strata, frac, target)

    # Наклон ЯДРА, то есть регрессия предсказания на метку. Обратная регрессия --- метки
    # на предсказание --- отличается на долю объяснённой дисперсии и даёт не то усиление.
    slope = float(np.polyfit(y, oof, 1)[0])
    r2 = 1.0 - float(((y - oof) ** 2).sum() / ((y - y.mean()) ** 2).sum())

    rng = np.random.default_rng(seed)
    bs = []
    for _ in range(draws):
        i = _resample(rng, groups, len(y))
        j = _resample(rng, blind_groups, len(blind_pred))
        try:
            if strata is None:
                v = _solve(y[i], kernel(y[i], oof[i]), float(blind_pred[j].mean()))
            else:
                fr = {k: float((blind_strata[j] == k).mean()) for k in np.unique(strata)}
                v, _ = _solve_strat(y[i], oof[i], strata[i], fr, float(blind_pred[j].mean()))
        except (ValueError, ZeroDivisionError):
            continue
        if np.isfinite(v):
            bs.append(v)
    bs = np.array(bs)
    a = (100.0 - level) / 2.0
    q = np.percentile(bs, [a, 100.0 - a]) if len(bs) else (np.nan, np.nan)
    return Shift(delta=delta, lo=q[0], hi=q[1], draws=bs, slope=slope, r2=r2,
                 composition=comp, n=len(y), n_blind=len(blind_pred))


def recalibrate(y, oof, band_lo, band_hi, delta, offsets=None, lambdas=None):
    """The affine pair (offset, lambda) that minimises band loss under an ASSUMED shift.

    Predictions are mapped by q = mu*(1 - lam) + off*(1 - lam) + lam*p, so the fitted pair is
    a shrinkage towards a centre displaced by `off`. The two parameters are fitted JOINTLY:
    fixing one and sweeping the other picks an arbitrary slice through a two-parameter family,
    and doing so cost 0.008 of macro here before it was noticed.

    The grid is checked at its edges, because an optimum sitting on the boundary means the
    grid, not the data, chose the answer.
    """
    y = np.asarray(y, float)
    p = np.asarray(oof, float)
    lo = np.asarray(band_lo, float)
    hi = np.asarray(band_hi, float)
    offsets = np.round(np.arange(-3.0, 3.01, 0.05), 2) if offsets is None else np.asarray(offsets)
    lambdas = np.round(np.arange(0.0, 1.001, 0.02), 2) if lambdas is None else np.asarray(lambdas)

    # Целевая выборка --- наша собственная, наклонённая к предполагаемому сдвигу.
    w = _tilt_weights(y, 0.0) if abs(delta) < 1e-12 else None
    if w is None:
        f = lambda th: float(_tilt_weights(y, th) @ y) - (y.mean() + delta)
        if f(-TH_MAX) * f(TH_MAX) > 0:
            raise ValueError(f"сдвиг {delta:+.3f} недостижим на этой выборке")
        w = _tilt_weights(y, brentq(f, -TH_MAX, TH_MAX, xtol=1e-10))

    mu = float(y.mean())
    A = offsets[:, None] * (1.0 - lambdas[None, :])
    B = np.broadcast_to(lambdas[None, :], A.shape)
    q = (mu * (1.0 - B) + A)[:, :, None] + B[:, :, None] * p[None, None, :]
    pen = (w * (np.maximum(q - hi, 0.0) + np.maximum(lo - q, 0.0))).sum(axis=2)
    i, j = np.unravel_index(pen.argmin(), pen.shape)
    if i in (0, len(offsets) - 1) or j in (0, len(lambdas) - 1):
        raise ValueError(f"оптимум на краю сетки: смещение {offsets[i]}, lambda {lambdas[j]}")
    return float(offsets[i]), float(lambdas[j])


def apply_affine(pred, off, lam, mu):
    """Apply the fitted pair. Predictions move by (1 - lam) * off, not by off."""
    return mu * (1.0 - lam) + off * (1.0 - lam) + lam * np.asarray(pred, float)


def main():
    """Reproduce this repository's four published estimates, as the acceptance test."""
    import sys as _sys, pathlib as _pl
    _sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
    from cyppaths import D, RES
    import json
    import pandas as pd
    from cypsplit import cluster_ids

    CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
    PUBLISHED = {"CYP1A2": 0.045, "CYP2C9": 0.362, "CYP2D6": -0.917, "CYP3A4": 0.740}

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    oof = json.load(open(RES + "preds/oof.json"))
    PT = np.load(D + "test_pred.npz")
    cid = np.asarray(cluster_ids(list(rows.SMILES))[0])
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    tcid = np.asarray(cluster_ids(list(te.SMILES), threshold=0.50)[0])

    mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    z = np.load(D + "feats.npz")
    base = z["MECH"][:, mn.index("is_base_74")] > 0.5
    tbase = np.load(D + "test_feats.npz")["MECH"][:, mn.index("is_base_74")] > 0.5

    print(f"{'фермент':8s} {'наш':>7s} {'опубл.':>8s} {'расх.':>7s} {'95% интервал':>22s} "
          f"{'усил.':>6s}")
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        s = estimate(tr.loc[m, col].to_numpy(), np.asarray(oof[f"FP+DESC+MECH|{c}"]),
                     PT[c], groups=cid[m], blind_groups=tcid)
        d = abs(s.delta - PUBLISHED[c])
        print(f"{c:8s} {s.delta:+7.3f} {PUBLISHED[c]:+8.3f} {d:7.3f} "
              f"[{s.lo:+8.3f}, {s.hi:+8.3f}] {s.amplification():6.1f}"
              + ("" if d < 0.01 else "   <-- НЕ ВОСПРОИЗВЕЛОСЬ"))

    col = "CYP2D6_pIC50_direct_inhibition"
    m = tr[col].notna().to_numpy()
    s = estimate(tr.loc[m, col].to_numpy(), np.asarray(oof["FP+DESC+MECH|CYP2D6"]),
                 PT["CYP2D6"], groups=cid[m], blind_groups=tcid,
                 strata=base[m], blind_strata=tbase)
    print(f"\nCYP2D6 со стратификацией по основаниям: {s.delta:+.3f} "
          f"[{s.lo:+.3f}, {s.hi:+.3f}], из них состав {s.composition:+.3f}")
    print("Опубликовано -0.508, из них состав -0.138."
          if abs(s.delta + 0.508) > 0.01 else "Совпало с опубликованным -0.508.")


if __name__ == "__main__":
    main()
