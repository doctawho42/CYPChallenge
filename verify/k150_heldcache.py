"""Fill the chemeleon-probe cache for the held-out seeds, with k122's own function.

Reusing `k122_stack.emb_probe` rather than reimplementing it is the whole point: the cache for
seeds 0-3 was written by that function, so the held-out seeds are produced by the same device and
the confirmation compares like with like.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd
from cypsplit import butina_folds, fold_digest
from submit import CYPS
from k122_stack import emb_probe, EMB
rows = pd.read_csv(D + "rows.csv")
tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
        .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS],1)
M = ~np.isnan(Y)
# emb_probe takes a PATH, not an array: it calls load_any, which carries the alignment
# control for nameless files and verifies molecule_name for the ones that have it.
E = EMB["chemeleon"]
for s in (4,5,6,7):
    fold, _ = butina_folds(list(rows.SMILES), seed=s); dg = fold_digest(fold)
    P = emb_probe("chemeleon", E, rows, Y, M, fold, s, dg)
    print("  сид %d дайджест %s готов" % (s, dg), flush=True)
print("готово")
