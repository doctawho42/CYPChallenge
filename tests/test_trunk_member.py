"""Guards on the fifth ensemble member, which is the only one read from a file.

Item 121 added three guards to `src/submit.py` around the trunk member, and each of them
protects against something that would otherwise fail silently rather than loudly. They were
verified once by hand when they were written. That is not the same as verified.

Not covered here: the identity of the test-time path with `run_fold`, which is checked by
`src/trunk.py --check-test-path` and needs a torch training run — too slow for a test suite
that is meant to run on every push. What is covered is everything that costs milliseconds.
"""
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

from cyppaths import D, RES  # noqa: E402

ROWS = pathlib.Path(D) / "rows.csv"
FEATS = pathlib.Path(D) / "feats.npz"
TRUNK = pathlib.Path(RES) / "preds" / "trunk_twohead.json"

pytestmark = pytest.mark.skipif(
    not (ROWS.exists() and FEATS.exists() and TRUNK.exists()),
    reason="needs data/rows.csv, data/feats.npz and results/preds/trunk_twohead.json",
)


@pytest.fixture(scope="module")
def S():
    import submit
    return submit


@pytest.fixture(scope="module")
def data():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    from cypsplit import butina_folds
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float)
                  for c in ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]], 1)
    fold, _ = butina_folds(list(rows.SMILES))
    return y, ~np.isnan(y), fold


def test_blocks_round_trip(S):
    """The trunk takes FP/DESC/MECH apart because it log1p's the fingerprint. If the widths
    are ever guessed instead of read, the whole matrix shifts by a column and nothing errors."""
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    fp, de, me = S._trunk_blocks(X)
    assert np.array_equal(fp, z["FP"])
    assert np.array_equal(de, z["DESC"])
    assert np.array_equal(me, z["MECH"])


def test_blocks_reject_wrong_width(S):
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"]])
    with pytest.raises(SystemExit):
        S._trunk_blocks(X)


def test_oof_trunk_lengths(S, data):
    y, mask, fold = data
    P = S._oof_trunk(y, mask, fold)
    assert [len(p) for p in P] == [int(mask[:, e].sum()) for e in range(4)]


def test_oof_trunk_rejects_moved_split(S, data):
    """The guard that matters most. Every other member is recomputed against whatever split
    is current; this one is read from a file computed on seed 0. A moved split would leave it
    silently in-sample, and nothing else in the pipeline would notice."""
    y, mask, fold = data
    bad = fold.copy()
    bad[0] = (bad[0] + 1) % 5
    with pytest.raises(SystemExit):
        S._oof_trunk(y, mask, bad)


def test_clip_is_applied(S, data):
    """Item 79: without the clip, seed 0 contains one compound the trunk predicts at -360,
    and the affine pair does not absorb it. The measurement would be of that molecule."""
    y, mask, fold = data
    P = S._oof_trunk(y, mask, fold)
    for e in range(4):
        yy = y[mask[:, e], e]
        assert P[e].min() >= yy.min() - 2.0 - 1e-9
        assert P[e].max() <= yy.max() + 2.0 + 1e-9


def test_default_mode_has_no_trunk(S, data):
    """`ансамбль` must stay a four-member average. The fifth member is behind a flag because
    the choice of what to submit belongs to the team, and a default that drifts is not a
    choice anyone made."""
    y, mask, fold = data
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    called = []
    # Только проверка ветвления: тяжёлые члены заглушены, ловится сам факт вызова.
    saved = {n: getattr(S, n) for n in ("_oof_trunk", "_oof_one", "_oof_gp", "_oof_ridge")}
    stub = lambda *a, **k: [np.zeros(int(mask[:, e].sum())) for e in range(4)]
    try:
        S._oof_trunk = lambda *a, **k: (called.append(1), stub())[1]
        for n in ("_oof_one", "_oof_gp", "_oof_ridge"):
            setattr(S, n, stub)
        S.oof_predictions(X, y, mask, fold, "ансамбль")
        assert called == [], "режим «ансамбль» тронул пятый член"
        S.oof_predictions(X, y, mask, fold, "ансамбль5")
        assert called == [1], "режим «ансамбль5» не тронул пятый член"
    finally:
        for n, f in saved.items():
            setattr(S, n, f)
