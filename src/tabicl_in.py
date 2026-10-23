"""Build the second TabICL input bundle: E (frozen 600) for the PCA grid, S (surrogate_adme 2048).

S is the instrument, not a candidate.  The 1/n law measured today says no ninth member can clear
the floor at n=8, so the remaining question is whether the law's missing interaction term stays
missing for a DECORRELATED member.  Every other embedding on disk already has a member standing on
it inside the composition (chemeleon's probe, chemprop_medium's probe); surrogate_adme has none.

S is also the one embedding file with NO molecule_name array, so its row order is an assumption
rather than a check.  align_control is the control that must succeed: a probe on the true order has
to beat a probe on a rolled copy of the same matrix by more than 0.05 of Spearman.  A roll leaves
the matrix, its scale and its rank untouched and breaks only the correspondence to the labels, so
it exercises exactly the suspect clause.
"""
import json
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd
from cypsplit import butina_folds, fold_digest
from submit import CYPS
from k121_chemeleon import align_control
rows = pd.read_csv(D + "rows.csv")
tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
        .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS],1)
M = ~np.isnan(Y)
fold0, _ = butina_folds(list(rows.SMILES), seed=0)

E = np.load(D+"emb_chemprop_medium.npz", allow_pickle=True)
assert list(E["molecule_name"]) == list(rows.Molecule_Name), "medium: порядок не совпал"
print("  E  (%s) имена сверены со rows.csv поэлементно" % (E["train"].shape,))

Sz = np.load(D+"emb_chemprop_surrogate_adme.npz", allow_pickle=True)
S = Sz["train"].astype(np.float64)
assert S.shape[0] == len(rows), S.shape
print("  S  (%s) имён в файле нет, ключи: %s" % (S.shape, list(Sz.files)))
ok = align_control(S, Y, M, fold0, e=0)
print("  S  КОНТРОЛЬ выравнивания (истинный порядок против свёрнутого): %s" % ok)
if not ok:
    raise SystemExit("ОТКАЗ: порядок строк surrogate_adme не подтверждён")

d = {"E": E["train"].astype(np.float32), "S": S.astype(np.float32),
     "Y": Y, "cyps": np.array(CYPS)}
digs = []
for s in (0,1,2,3):
    f, _ = butina_folds(list(rows.SMILES), seed=s)
    d["fold%d" % s] = f; digs.append(fold_digest(f))
assert digs[0] == "2d93c19815e14261", digs[0]
d["digests"] = np.array(digs)
np.savez_compressed(D+"tabicl_in2.npz", **d)
import os
print("\n  записано %stabicl_in2.npz, %.1f МБ, дайджесты %s"
      % (D, os.path.getsize(D+"tabicl_in2.npz")/1e6, ",".join(x[:8] for x in digs)))
