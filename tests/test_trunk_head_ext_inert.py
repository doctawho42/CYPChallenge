"""The third head must be INERT, or every number src/trunk.py has published moves.

`Net` grew a `head_ext` so that the primary, the screen and an external source can each keep
their own output scale (item 318). The trunk's published numbers were taken with two heads, so
the condition for that change to be safe is exact: with `ext=None` and `lam_ext=0`, `run_fold`
must return predictions BIT FOR BIT identical to the two-head file's.

Why this needs a test rather than a one-off check. The hazard is not in the arithmetic, it is in
the random number stream. `nn.Linear` draws from the CPU generator, and on `device="cpu"` so does
every dropout mask in the training loop — so a third head that is merely constructed LAST still
leaves the weights alone and moves every prediction, through the masks. `Net.__init__` hands the
generator back after building `head_ext` for exactly this reason. A future edit that reorders the
constructor, drops the wind-back, or adds a fourth module would break it silently, and silently is
the point: nothing would fail, the numbers would simply stop matching the document.

Two things make this test able to fail rather than pass vacuously. It compares against the
PRISTINE two-head file recovered from git rather than against itself, and it asserts first that
the two modules genuinely differ. It runs on cpu deliberately: on mps the dropout masks come from
the MPS generator, which the CPU-side initialisation never touches, so mps would pass either way.

If this fails after an intentional architecture change, do not relax it — re-derive the published
numbers in the same commit, as `tests/test_split.py` says for the fold digest.
"""
import importlib.util
import pathlib
import subprocess
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from cyppaths import D  # noqa: E402

ROOT = pathlib.Path(__file__).resolve().parents[1]
EPOCHS_FAST = 5
FOLD = 1


def _load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def modules(tmp_path_factory):
    """The current trunk beside the two-head one, recovered from git."""
    pytest.importorskip("torch")
    r = subprocess.run(["git", "show", "HEAD:src/trunk.py"],
                       capture_output=True, text=True, cwd=ROOT)
    if r.returncode != 0:
        pytest.skip(f"cannot recover the committed trunk: {r.stderr.strip()[:80]}")
    old_path = tmp_path_factory.mktemp("trunk") / "pristine_trunk.py"
    old_path.write_text(r.stdout)
    new = _load("trunk_head_ext_new", ROOT / "src" / "trunk.py")
    old = _load("trunk_head_ext_old", old_path)
    return new, old


def test_modules_actually_differ(modules):
    """Non-vacuity: 'identical' must not be able to mean 'the same file twice'.

    If HEAD already carries head_ext the comparison below is meaningless, so say so loudly.
    """
    new, old = modules
    assert hasattr(new.Net(8), "head_ext"), "src/trunk.py has no head_ext to test"
    if hasattr(old.Net(8), "head_ext"):
        pytest.skip("HEAD already has head_ext; the pristine two-head baseline is gone")


def test_head_ext_is_bit_identical_when_inert(modules):
    """ext=None, lam_ext=0 must reproduce the two-head predictions exactly, on cpu."""
    new, old = modules
    rows_csv = pathlib.Path(D + "rows.csv")
    if not rows_csv.exists():
        pytest.skip("data/rows.csv absent; run src/feats.py first")

    from cypsplit import butina_folds

    new.EPOCHS = old.EPOCHS = EPOCHS_FAST
    X, y, scr = old.load("DESC+MECH")[:3]
    fold, _ = butina_folds(list(pd.read_csv(rows_csv).SMILES))

    for lam in (0.0, 3.0):
        a = np.asarray(old.run_fold(X, y, scr, fold, FOLD, lam, 0, "cpu", mode="twohead"))
        b = np.asarray(new.run_fold(X, y, scr, fold, FOLD, lam, 0, "cpu", mode="twohead",
                                    ext=None, lam_ext=0.0))
        assert a.shape == b.shape
        diff = float(np.abs(a - b).max())
        assert diff == 0.0, (
            f"lam={lam}: the third head is NOT inert, max|diff| = {diff:.3e}. "
            "Published trunk numbers have moved; see this file's docstring."
        )


def test_lam_ext_without_a_block_refuses(modules):
    """A lambda with no target block would be a silent zero, so it must raise instead.

    masked_mse over an empty mask returns 0, so the arm would run, cost nothing and mean
    nothing — the shape CLAUDE.md calls a query that could not have succeeded.
    """
    new, _ = modules
    rows_csv = pathlib.Path(D + "rows.csv")
    if not rows_csv.exists():
        pytest.skip("data/rows.csv absent; run src/feats.py first")

    from cypsplit import butina_folds

    new.EPOCHS = EPOCHS_FAST
    X, y, scr = new.load("DESC+MECH")[:3]
    fold, _ = butina_folds(list(pd.read_csv(rows_csv).SMILES))
    with pytest.raises(SystemExit):
        new.run_fold(X, y, scr, fold, FOLD, 0.0, 0, "cpu", mode="twohead",
                     ext=None, lam_ext=3.0)
