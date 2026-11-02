"""Several frozen-encoder members must be expressible, and adding the capacity must change nothing.

WHY THIS EXISTS. Item 346 adopted a composition `src/submit.py` could not build: six members plus
`chemeleon` plus the gbm head, +0.0047 of macro rank and -0.0085 of pair at sign 4/4, called there
"the first adoptable gain since item 324". `_combine` held exactly one `extra`, so the adopted
composition had no route to a submission file at all, and the board has never seen it -- item 344
tabulates nine uploads and the ninth predates items 345 to 351.

The capacity is added as a SECOND parameter, `extras`, rather than by widening `extra`, and this
pins both halves of what that has to mean:

  * nothing moves. `extras=None` and `extra=<one member>` must return what they returned before
    the parameter existed, element for element, because the file on the board was built by those
    paths and `tests/test_probe_member_inert.py` compares the no-extra path against the pre-probe
    module recovered from git.
  * the new members must reach EVERY enzyme. That is the same hazard item 324's test was written
    for, and it is not the arithmetic of a mean: `_keep` reads SOLO, which drops members on three
    enzymes of four. Members routed through that filter would vanish on CYP1A2, CYP2C9 and CYP3A4
    and survive only on CYP2D6 -- the enzyme with ~zero measured gain -- and the run would finish
    cleanly, write a file and mean nothing.

The synthetic members are distinct constants per member, so an implementation that returned any
single member, or that silently dropped one, fails. Non-vacuity is asserted rather than assumed:
the extras path must differ from the no-extras path on all four enzymes before the equalities
below are worth anything.

WHAT THIS DOES NOT COVER, said here so nobody reads a green test as a shipped composition:
`main()` is NOT wired to build more than one extra member, and the gbm head of item 335 has no
device in this file at all. What is pinned here is the combination and the loader, which is where
the hazard lives; the wiring is a separate change and `docs/board_gap.md` says so.
"""
import inspect
import pathlib
import sys

import numpy as np

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import submit  # noqa: E402

N = 7
MEMBERS = ("поферментно", "пул", "GP", "гребневая", "ствол")


def _parts():
    """Five members, each a distinct constant, so a dropped or duplicated one is visible."""
    return [(k, [np.full(N, 10.0 * (i + 1) + e, float) for e in range(len(submit.CYPS))])
            for i, k in enumerate(MEMBERS)]


def _extra(v):
    return [np.full(N, v + e, float) for e in range(len(submit.CYPS))]


def _kept_mean(parts, e):
    return np.mean([P[e] for k, P in parts if submit._keep(e, k)], axis=0)


def test_no_extras_is_the_shipped_expression():
    parts = _parts()
    for e in range(len(submit.CYPS)):
        assert np.array_equal(submit._combine(parts, e), _kept_mean(parts, e))
        assert np.array_equal(submit._combine(parts, e, extras=None), _kept_mean(parts, e))
        assert np.array_equal(submit._combine(parts, e, extras=[]), _kept_mean(parts, e))


def test_single_extra_path_is_unchanged():
    parts, X = _parts(), _extra(100.0)
    for e in range(len(submit.CYPS)):
        want = np.mean([P[e] for k, P in parts if submit._keep(e, k)] + [X[e]], axis=0)
        assert np.array_equal(submit._combine(parts, e, X), want)
        assert np.array_equal(submit._combine(parts, e, extra=X, extras=None), want)


def test_extras_reach_every_enzyme_and_are_not_filtered_by_solo():
    parts, X, Y = _parts(), _extra(100.0), _extra(200.0)
    for e in range(len(submit.CYPS)):
        kept = [P[e] for k, P in parts if submit._keep(e, k)]
        want = np.mean(kept + [X[e], Y[e]], axis=0)
        got = submit._combine(parts, e, extras=[X, Y])
        assert np.array_equal(got, want), submit.CYPS[e]
        # non-vacuity: the control must be able to fire on every enzyme
        assert not np.array_equal(got, _kept_mean(parts, e)), submit.CYPS[e]
        # and the DENOMINATOR must be len(kept) + 2. Asserting `len(kept)+2 == len(kept)+2`
        # would be a control that cannot fail, which is item 319's named defect; so the two
        # wrong denominators are computed and required to differ from what came back.
        total = sum(kept) + X[e] + Y[e]
        assert not np.array_equal(got, total / (len(kept) + 1)), submit.CYPS[e]
        assert not np.array_equal(got, total / (len(kept) + 3)), submit.CYPS[e]
        assert np.allclose(got, total / (len(kept) + 2)), submit.CYPS[e]


def test_extra_and_extras_compose_with_extra_first():
    parts, X, Y = _parts(), _extra(100.0), _extra(200.0)
    for e in range(len(submit.CYPS)):
        assert np.array_equal(submit._combine(parts, e, extra=X, extras=[Y]),
                              submit._combine(parts, e, extras=[X, Y]))


def test_member_count_scales_with_extras():
    """Three extras must average over kept+3, not over kept+1 or kept."""
    parts = _parts()
    ex = [_extra(100.0), _extra(200.0), _extra(300.0)]
    for e in range(len(submit.CYPS)):
        kept = [P[e] for k, P in parts if submit._keep(e, k)]
        want = np.mean(kept + [x[e] for x in ex], axis=0)
        assert np.array_equal(submit._combine(parts, e, extras=ex), want)


def test_probe_loader_defaults_still_point_at_the_shipped_pair():
    """Parameterising the loader must not change WHICH embedding ships by default."""
    sig = inspect.signature(submit._probe_load)
    for name in ("emb", "emb_te"):
        assert sig.parameters[name].default is None, name
    assert submit.PROBE_EMB.endswith("emb_chemprop_medium.npz")
    assert submit.PROBE_EMB_TEST.endswith("emb_chemprop_medium_test.npz")
    assert sig.parameters["tag"].default == submit.PROBE_KIND


def test_oof_predictions_threads_extras():
    assert "extras" in inspect.signature(submit.oof_predictions).parameters
