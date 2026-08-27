# Data

Not stored in git: these are the published challenge files, they are freely available and
together come to about seven megabytes.

```bash
bash data/fetch.sh
```

Five files should then be here.

| File | Rows | What it is |
|---|---|---|
| `cyp-challenge-TRAIN_inhibition.csv` | 4905 | direct-inhibition pIC50, interval bounds, standard error |
| `cyp-challenge-TRAIN_TDI.csv` | 6145 | the same for the pre-incubation condition, plus TDI labels |
| `cyp-challenge-TRAIN_Emax.csv` | 6146 | limiting depth of suppression, both arms |
| `cyp-challenge-single-concentration-TRAIN.csv` | 17504 | single-point screen at 49.5 uM |
| `cyp-challenge-TEST-BLINDED.csv` | 750 | the test set, name and SMILES only |

Derived files are built by `src/feats.py` and are not committed either:

- `feats.npz` — the feature matrices `FP` (2048), `DESC` (217), `MECH` (30);
- `rows.csv` — **the row order of the feature matrix**. `feats.py` drops molecules RDKit
  cannot parse, so every consumer realigns its table against this file;
- `desc_names.csv`, `mech_names.csv` — feature names, for reading the tables;
- `clusters_*.npz` — the Butina clustering cache written by `cypsplit.py`. The filename
  carries a hash of the SMILES and every clustering parameter, so a stale hit is not
  possible; delete them freely.

Saved out-of-fold predictions live in `results/preds/` and **are** committed: recomputing
them takes about an hour.
