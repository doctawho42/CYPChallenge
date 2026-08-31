"""Features for the NCGC panel, built by the same code and pinned to the same columns.

Why a separate script rather than a line inside the ablation. `src/feats.py` selects its
descriptor columns by `isna().mean() < 0.05`, a filter fitted on our 4905 training molecules.
Recomputing that filter on 13126 foreign structures would silently produce a different column
set, and a check on the column *count* would not catch a permutation either -- which is why
`build()` takes `desc_names` and `mech_names` and reindexes by name. The test set already goes
through that path; the panel must too, or its rows would be describing different quantities in
the same positions.

The panel reaches further into chemical space than anything here has: 13126 structures against
our 4905, overlapping in 183 to 319. So `_guard` is expected to fire -- Ipc grows exponentially
with molecule size and overflows float32 on public compounds while our own set tops out at
5.1e14. Its warnings are the point of running this separately and reading the output.

Reads data/ncgc/panel.csv (written by src/ncgcmerge.py). Writes data/ncgc/feats.npz with the
canonical SMILES in the row order of the matrices. ~5 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D

import time

import numpy as np
import pandas as pd

from feats import build

t0 = time.time()
panel = pd.read_csv(D + "ncgc/panel.csv", low_memory=False)
smiles = sorted(panel.smiles.dropna().unique())
print(f"уникальных структур в панели: {len(smiles)}")

dn = [l.strip() for l in open(D + "desc_names.csv")]
mn = [l.strip() for l in open(D + "mech_names.csv")]
print(f"колонки пинуются по именам: DESC {len(dn)}, MECH {len(mn)}\n")

FP, dsc, M, ok = build(smiles, desc_names=dn, mech_names=mn)
kept = [smiles[i] for i in ok]
print(f"\nRDKit разобрал {len(kept)} из {len(smiles)}")

np.savez_compressed(D + "ncgc/feats.npz", FP=FP, DESC=dsc.to_numpy(np.float32),
                    MECH=M.to_numpy(np.float32), SMILES=np.array(kept, dtype=object))
print(f"FP {FP.shape} DESC {dsc.shape} MECH {M.shape}   ({time.time()-t0:.0f} с)")
print(f"сохранено: {D}ncgc/feats.npz")

z = np.load(D + "feats.npz")
print(f"\nсверка ширины с обучающей матрицей: FP {z['FP'].shape[1]}=={FP.shape[1]}, "
      f"DESC {z['DESC'].shape[1]}=={dsc.shape[1]}, MECH {z['MECH'].shape[1]}=={M.shape[1]}")
assert (z["FP"].shape[1], z["DESC"].shape[1], z["MECH"].shape[1]) == \
       (FP.shape[1], dsc.shape[1], M.shape[1]), "ширина не совпала --- матрицы несовместимы"
print("совпало")
