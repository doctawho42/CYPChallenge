"""The sixth member must be INERT unless --probe is given, or every number the submission has
published moves.

`src/submit.py` grew a RidgeCV head on a frozen chemprop encoder (item 324). The file on the
board was built by the five-member path, so the condition for that change to be safe is exact:
with `extra=None` the combination `oof_predictions` performs must be the one the pre-probe file
performed, element for element, and the flag must default to off.

Why this needs a test rather than a one-off check. The hazard is not the arithmetic of a mean,
it is the FILTER in front of it. `_keep` reads SOLO, which names which of our own five members
each enzyme keeps -- and three enzymes of four drop members there. A sixth member routed through
that filter would vanish on CYP1A2, CYP2C9 and CYP3A4 and survive only on CYP2D6, which is the
enzyme with ~zero measured gain; the run would finish cleanly, write a file, and mean nothing.
So this pins both halves: extra=None changes nothing, and extra=<probe> reaches every enzyme.

Two things make it able to fail rather than pass vacuously. It compares against the PRE-PROBE
submit.py recovered from git by content hash, not against itself, and it asserts first that the
two modules genuinely differ. The synthetic members are drawn from a fixed generator and differ
from each other, so an implementation that returned any single member would fail.

If this fails after an intentional change to the composition, do not relax it -- re-derive the
submitted file in the same commit, as tests/test_split.py says for the fold digest.
"""
import importlib.util
import pathlib
import subprocess
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

# The pre-probe submit.py, pinned BY CONTENT. A blob hash is content-addressed and cannot
# drift the way `git show HEAD:src/submit.py` did in the trunk test's first version, where the
# baseline became the file under test the moment it was committed.
PRISTINE_PREPROBE_BLOB = "cedcd8d4bc825b1bbf5765d6c2e24ae370298b7c"

N = 37          # rows per enzyme in the synthetic members; small and arbitrary
KINDS = ("поферментно", "пул", "GP", "гребневая", "ствол")


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def modules(tmp_path_factory):
    """The current submit beside the pre-probe one, recovered from git."""
    r = subprocess.run(["git", "cat-file", "-p", PRISTINE_PREPROBE_BLOB],
                       capture_output=True, text=True, cwd=ROOT)
    if r.returncode != 0:
        pytest.skip(f"pinned pre-probe blob unreachable: {r.stderr.strip()[:80]}")
    old_path = tmp_path_factory.mktemp("submit") / "pristine_submit.py"
    old_path.write_text(r.stdout, encoding="utf-8")
    try:
        new = _load("submit_probe_new", ROOT / "src" / "submit.py")
        old = _load("submit_probe_old", old_path)
    except (ImportError, SystemExit) as exc:      # organisers' submodule absent
        pytest.skip(f"src/submit.py not importable here: {exc}")
    return new, old


@pytest.fixture(scope="module")
def parts():
    """Five synthetic members, one array per enzyme each, all different from one another."""
    g = np.random.default_rng(20260925)
    return [(k, [g.normal(5.0, 1.0, N) for _ in range(4)]) for k in KINDS]


@pytest.fixture(scope="module")
def probe():
    g = np.random.default_rng(11)
    return [g.normal(5.0, 1.0, N) for _ in range(4)]


def _pristine_combine(old, parts, e):
    """The pre-probe expression, executed from the pre-probe module rather than retyped.

    `oof_predictions` is called with its own `oof_members` and `dz_pass` replaced, so the
    expression under test is the one that file actually carries, not a copy of it here.
    """
    old.oof_members = lambda *a, **k: parts
    old.dz_pass = lambda p, *a, **k: p
    return old.oof_predictions(None, None, None, None, "ансамбль5")[e]


def test_modules_actually_differ(modules):
    """Non-vacuity: 'identical' must not be able to mean 'the same file twice'.

    This FAILS rather than skips. A skip is what silently disarmed the trunk test's first
    version, so the condition that would once have skipped is the failure here.
    """
    new, old = modules
    assert hasattr(new, "_combine"), "src/submit.py has no _combine to test"
    assert not hasattr(old, "_combine"), (
        "the pinned baseline already carries the probe refactor, so the comparison below "
        "would be the file against itself. Re-pin PRISTINE_PREPROBE_BLOB."
    )
    assert not hasattr(old, "PROBE_ALPHAS"), "the pinned baseline already carries the probe"


def test_probe_is_off_by_default(modules):
    """The flag exists, is a switch, and is off unless asked for."""
    new, _ = modules
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--probe", action="store_true")
    assert ap.parse_args([]).probe is False
    assert ap.parse_args(["--probe"]).probe is True
    src = (ROOT / "src" / "submit.py").read_text(encoding="utf-8")
    assert '"--probe", action="store_true"' in src, "--probe is no longer a bare switch"


def test_extra_none_is_bit_identical(modules, parts):
    """extra=None must reproduce the pre-probe combination exactly, on every enzyme."""
    new, old = modules
    new.oof_members = lambda *a, **k: parts
    new.dz_pass = lambda p, *a, **k: p
    got = new.oof_predictions(None, None, None, None, "ансамбль5")
    for e in range(4):
        want = _pristine_combine(old, parts, e)
        diff = float(np.abs(np.asarray(got[e]) - np.asarray(want)).max())
        assert diff == 0.0, (
            f"enzyme {e}: the sixth member is NOT inert, max|diff| = {diff:.3e}. "
            "The submitted numbers have moved; see this file's docstring."
        )


def test_combination_can_change_at_all(modules, parts, probe):
    """The comparison above must be capable of failing: with a probe it does change."""
    new, old = modules
    new.oof_members = lambda *a, **k: parts
    new.dz_pass = lambda p, *a, **k: p
    with_probe = new.oof_predictions(None, None, None, None, "ансамбль5", extra=probe)
    for e in range(4):
        want = _pristine_combine(old, parts, e)
        assert float(np.abs(np.asarray(with_probe[e]) - np.asarray(want)).max()) > 0.0


def test_probe_reaches_every_enzyme(modules, parts, probe):
    """The trap this file exists for: SOLO must not filter the sixth member away.

    Checked by arithmetic rather than by inspection -- the mean with the probe must equal the
    mean of (kept members + probe), which on an enzyme where the probe had been dropped would
    instead equal the mean of the kept members alone.
    """
    new, _ = modules
    for e, c in enumerate(new.CYPS):
        kept = [P[e] for k, P in parts if new._keep(e, k)]
        assert len(kept) == len(new.SOLO.get(c, KINDS)), f"{c}: SOLO fixture out of step"
        want = np.mean(kept + [probe[e]], axis=0)
        got = new._combine(parts, e, probe)
        assert float(np.abs(got - want).max()) == 0.0, f"{c}: probe not in the mean"
        assert float(np.abs(got - np.mean(kept, axis=0)).max()) > 0.0, (
            f"{c}: the probe was filtered out by SOLO -- it would be inert on this enzyme"
        )
