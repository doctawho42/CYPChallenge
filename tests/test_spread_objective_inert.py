"""The label SPREAD must be expressible in the shrinkage objective, and adding it must change nothing.

`src/submit.py:fit_shrinkage` picks the per-enzyme (shift, lambda) by minimising the soft error
against OUR OWN labels, reweighted by `reweight.tilt` so their MEAN sits `delta` higher. That was
the only distribution shift the objective could express. The leaderboard's MAE and R2 columns bound
the test's label SPREAD -- sd(y_test) >= MAE/sqrt(1-R2) by Jensen, agreeing across 273 entrant rows
inside each enzyme -- and it is 1.02 to 1.48 times our training spread, widest on CYP2D6. A lambda
fitted against the narrow marginal keeps a compression the wider test punishes.

So `reweight.widen` is the exact analogue of `tilt` for the second moment, and `fit_shrinkage` grew
a `spread` argument. This pins the two things that have to be true:

  * `spread = 1.0` returns ones EXACTLY, so every number the submission has published is untouched;
  * lambda is monotone in the spread, which is the mechanism the whole claim rests on. Asserting it
    on synthetic data rather than on the real tables keeps this test in milliseconds and makes it
    run on a checkout with no data fetched.

Non-vacuity is asserted, not assumed: the compressed-prediction fixture must give a lambda that
actually moves, or a broken `widen` that returned ones for every ratio would pass silently.
"""
import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

from reweight import tilt, widen  # noqa: E402
import submit  # noqa: E402


def _fixture(n=600, seed=0):
    """Labels with bands, and a prediction mildly compressed relative to their spread.

    The compression is 0.75 rather than the 0.32-0.61 the real predictions carry, and the band
    is wide: both deliberately, so that the optimal lambda sits INSIDE `submit.GRID` at
    spread 1.0 and has room to rise. With the real degree of compression the optimum pins at
    the grid's 2.0 edge for every spread, and a saturated parameter cannot show a trend --
    which is the first version of this test and is how the fixture got tuned.
    """
    rng = np.random.default_rng(seed)
    y = rng.normal(5.0, 1.0, n)
    p = y.mean() + 0.75 * (y - y.mean()) + rng.normal(0, 0.35, n)
    half = 0.35
    return y, p, y - half, y + half


def test_widen_at_one_is_exactly_ones():
    y, _p, _lo, _hi = _fixture()
    assert np.array_equal(widen(y, 1.0), np.ones_like(y))
    assert np.array_equal(widen(y, 1.0 + 1e-15), np.ones_like(y))


@pytest.mark.parametrize("ratio", [1.1, 1.25, 1.5])
def test_widen_hits_the_spread_and_holds_the_mean(ratio):
    y, _p, _lo, _hi = _fixture()
    w = widen(y, ratio)
    m = float(np.average(y, weights=w))
    sd = float(np.sqrt(np.average((y - m) ** 2, weights=w)))
    assert abs(sd - ratio * y.std()) < 1e-6 * max(1.0, ratio * y.std())
    assert abs(m - y.mean()) < 1e-8
    # and it must differ from tilt, which moves the mean instead
    assert not np.allclose(w, tilt(y, 0.3))


def test_widen_refuses_a_spread_the_labels_cannot_express():
    y, _p, _lo, _hi = _fixture()
    with pytest.raises(ValueError):
        widen(y, 50.0)


def _fit(spread, delta=(0.0, 0.0, 0.0, 0.0)):
    y, p, lo, hi = _fixture()
    Y = np.stack([y] * 4, 1)
    submit.LO = np.stack([lo] * 4, 1)
    submit.HI = np.stack([hi] * 4, 1)
    M = np.ones_like(Y, dtype=bool)
    return submit.fit_shrinkage([p] * 4, Y, M, delta=delta, spread=spread)


def test_spread_one_is_inert():
    """The default path must be bit-for-bit what it was before `spread` existed."""
    base = _fit((1.0, 1.0, 1.0, 1.0))
    # the pre-change call signature, with spread left at its default
    y, p, lo, hi = _fixture()
    Y = np.stack([y] * 4, 1)
    submit.LO = np.stack([lo] * 4, 1)
    submit.HI = np.stack([hi] * 4, 1)
    M = np.ones_like(Y, dtype=bool)
    default = submit.fit_shrinkage([p] * 4, Y, M)
    assert base == default
    for d in (0.0, 0.4):
        assert _fit((1.0,) * 4, delta=(d,) * 4) == _fit((1.0,) * 4, delta=(d,) * 4)


def test_lambda_is_monotone_in_the_spread():
    lams = [_fit((r,) * 4)[0][0] for r in (1.0, 1.15, 1.3, 1.45)]
    assert all(b >= a for a, b in zip(lams, lams[1:])), lams
    # non-vacuity: it must actually move, or a widen() that returned ones would pass
    assert lams[-1] > lams[0] + 1e-9, lams
    # and it must not be reading the grid's edge, where no trend could be visible
    assert lams[-1] < submit.GRID[-1] - 1e-9, lams


def test_spread_and_delta_compose():
    """Both reweightings must apply, so the fit differs from either alone."""
    a = _fit((1.3,) * 4, delta=(0.0,) * 4)
    b = _fit((1.0,) * 4, delta=(0.4,) * 4)
    c = _fit((1.3,) * 4, delta=(0.4,) * 4)
    assert c != a and c != b
