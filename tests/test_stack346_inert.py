"""Item 346's stack must be expressible in `src/submit.py`, and adding it must change nothing.

Item 346 adopted a composition this file could not build: the shipped six, plus a ridge head on a
second frozen checkpoint, plus a gbm head on the first probe's columns --- +0.0047 of macro rank
and -0.0085 of pair at sign 4/4, both above their floors, and the only measured gain over what
ships. Item 351 records why it never reached the board: `main()` threaded ONE extra member and ran
ONE head, so the stack could be measured and not submitted.

`_combine` already took `extras` as a plural list. This pins the three things that have to be true
of the threading:

  * an EMPTY `extras` is arithmetically the old call, so every number the submission has published
    is untouched --- and non-vacuity is asserted, because a `_combine` that ignored `extras`
    entirely would pass that on its own;
  * the new flags are OFF by default, and a half-built stack is refused rather than averaged: its
    gain was measured as a whole, and item 361 priced a member SLOT at about 0.006 of rank, so
    members seven and eight are not interchangeable with member six;
  * the gbm head carries `verify/k114_gate93.py`'s hyper-parameters and NOT `gbm_reg()`'s. Those
    differ in exactly one keyword --- `l2_regularization`, 0.0 against 1.0 --- and items 335 and 346
    measured the member on 0.0. Reusing `gbm_reg()` would ship a member nobody measured, and the
    diff would read like a simplification.
"""
import pathlib
import re
import subprocess
import sys

import numpy as np

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import submit as S  # noqa: E402


def _parts(n=9, seed=0):
    """Member predictions shaped as `_combine` wants them: (kind, per-enzyme list)."""
    rng = np.random.default_rng(seed)
    kinds = ["поферментно", "пул", "GP", "гребневая", "ствол"]
    return [(k, [rng.normal(size=n) for _ in range(len(S.CYPS))]) for k in kinds]


def test_empty_extras_is_exactly_the_old_call():
    parts = _parts()
    extra = [np.arange(9.0) for _ in range(len(S.CYPS))]
    for e in range(len(S.CYPS)):
        a = S._combine(parts, e, extra)
        b = S._combine(parts, e, extra, extras=[])
        assert np.array_equal(a, b), f"{S.CYPS[e]}: пустой extras сдвинул среднее"
        c = S._combine(parts, e, extra, extras=None)
        assert np.array_equal(a, c), f"{S.CYPS[e]}: extras=None сдвинул среднее"


def test_extras_actually_enter_the_mean():
    """Non-vacuity: a `_combine` that ignored `extras` would pass the test above by itself."""
    parts = _parts()
    extra = [np.zeros(9) for _ in range(len(S.CYPS))]
    more = [[np.full(9, 10.0) for _ in range(len(S.CYPS))]]
    for e in range(len(S.CYPS)):
        a = S._combine(parts, e, extra)
        b = S._combine(parts, e, extra, extras=more)
        assert not np.allclose(a, b), f"{S.CYPS[e]}: extras не вошёл в среднее"
        # The count is NOT 5 + 1: `_keep` admits three members on CYP1A2 and two on CYP2C9
        # and CYP3A4 (items 282--285), so it is taken from `_keep` rather than assumed. That
        # is also the interaction `_combine`'s docstring warns about -- a new member must NOT
        # be filtered through SOLO, or it would vanish on three enzymes of four.
        n = sum(1 for k, _P in parts if S._keep(e, k)) + 1      # + the `extra`
        assert np.allclose(b, (a * n + 10.0) / (n + 1)), \
            f"{S.CYPS[e]}: вес нового члена не равный при {n} уже входящих"
        assert n < len(parts) + 1 or e == 2, \
            f"{S.CYPS[e]}: _keep пустил все пять — проверьте SOLO, а не этот тест"


def test_gbm_head_is_k114s_and_not_ablates():
    """The one keyword that decides whether the shipped member is the measured one."""
    assert set(S.GBM_KW) == {"max_iter", "learning_rate", "random_state"}, S.GBM_KW
    assert S.GBM_KW["max_iter"] == 300 and S.GBM_KW["learning_rate"] == 0.06
    # sklearn's default, which is what k114 ran on. Naming it here is the point.
    assert "l2_regularization" not in S.GBM_KW
    got = S.HistGradientBoostingRegressor(**S.GBM_KW).l2_regularization
    assert got == 0.0, f"голова получила l2={got}, а мерили на 0.0"
    assert S.gbm_reg().l2_regularization == 1.0, "gbm_reg() перестал отличаться — проверьте оба"


def test_gbm_params_match_k114_source():
    """Cross-file: if either side's numbers move, this fails loudly rather than drifting."""
    src = (ROOT / "verify" / "k114_gate93.py").read_text(encoding="utf-8")
    m = re.search(r'elif kind == "gbm":(.{0,400})', src, re.S)
    assert m, "ветка gbm в k114 не найдена — переименовали?"
    branch = m.group(1)
    assert "HistGradientBoostingRegressor" in branch
    assert "max_iter=300" in branch and "learning_rate=0.06" in branch
    assert "l2_regularization" not in branch, "k114 завёл l2 — submit.py надо сверить"


def _run(args):
    return subprocess.run([sys.executable, str(ROOT / "src" / "submit.py")] + args,
                          capture_output=True, text=True, timeout=180, cwd=str(ROOT))


def test_flags_are_off_by_default():
    h = _run(["--help"])
    assert h.returncode == 0, h.stderr[-400:]
    for flag in ("--probe2", "--gbm", "--stack346"):
        assert flag in h.stdout, f"{flag} нет в справке"
    # They are store_true, so absence from the shipped command line means off. The shipped
    # provenance records "все флаги по умолчанию", so this is the guard on that sentence.
    assert "ВЫКЛЮЧЕНО по умолчанию" in h.stdout


def test_half_built_stack_is_refused():
    for args in (["--probe2"], ["--gbm"], ["--probe2", "--gbm"]):
        r = _run(args)
        assert r.returncode != 0, f"{args} не был отвергнут"
        assert "требуют --probe" in (r.stdout + r.stderr), (r.stdout + r.stderr)[-400:]
