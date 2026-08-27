"""Golden-value guard on the cross-validation split.

Every reported number depends on the fold assignment, and numbers printed by different
scripts are comparable only while their folds agree. These digests were captured from
the original inline recipe before it was factored into cypsplit.py, so a failure here
means the split moved and the tables in docs/ no longer describe the current code.

If a change to the split is intentional: update the digests in the same commit, say so
in the message, and re-run the pipeline. Never adjust a digest to make CI green.
"""
import pathlib
import sys

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))

from cyppaths import D  # noqa: E402
from cypsplit import butina_folds, fold_digest  # noqa: E402

ROWS = pathlib.Path(D) / "rows.csv"

pytestmark = pytest.mark.skipif(
    not ROWS.exists(),
    reason="needs data/rows.csv - run `bash data/fetch.sh && python src/feats.py` first",
)

# The canonical split: fingerprints at 2048 bits, threshold 0.35, seed 0.
# Used by ablate, base, decision, range, rangectl, tdibase, tdiprob, f12_cvhard.
CANONICAL = "2d93c19815e14261"
N_CLUSTERS = 4703
FOLD_SIZES = [1016, 938, 992, 933, 1026]

# mech2/mech3/mech4 cluster at 1024 bits, so they have always had their own folds.
# Pinned too, so that difference stays deliberate rather than drifting further.
MECH_1024 = "d7b7ca15ecafe232"

# verify/f3_seeds.py measures how much the split alone moves the answer.
SEEDS = {1: "b26e229cfaace213", 2: "14f485f315fd992c", 3: "1193b75ee907b239"}


@pytest.fixture(scope="module")
def smiles():
    return list(pd.read_csv(ROWS).SMILES)


def test_canonical_split_is_unchanged(smiles):
    fold, n_clusters = butina_folds(smiles)
    assert n_clusters == N_CLUSTERS
    assert np.bincount(fold).tolist() == FOLD_SIZES
    assert fold_digest(fold) == CANONICAL


def test_mech_scripts_keep_their_own_1024_bit_split(smiles):
    fold, _ = butina_folds(smiles, fp_size=1024)
    assert fold_digest(fold) == MECH_1024
    assert fold_digest(fold) != CANONICAL, "mech folds must stay distinct from canonical"


@pytest.mark.parametrize("seed,digest", sorted(SEEDS.items()))
def test_seed_sweep_is_unchanged(smiles, seed, digest):
    fold, _ = butina_folds(smiles, seed=seed)
    assert fold_digest(fold) == digest


def test_clusters_never_straddle_folds(smiles):
    """The whole point of the split: a Butina cluster lands in exactly one fold."""
    from cypsplit import cluster_ids

    cid, _ = cluster_ids(smiles)
    fold, _ = butina_folds(smiles)
    df = pd.DataFrame({"cid": cid, "fold": fold})
    assert (df.groupby("cid").fold.nunique() == 1).all()
