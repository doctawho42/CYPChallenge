"""Feature-invariant fingerprints: item 90's finding carried over to the feature side.

The argument. Item 90 measured pharmacophore similarity beating substructure similarity as a
METRIC, 0.187 against 0.103. Item 96 then showed that does not transfer to features - Morgan bits
are strong features even though Morgan distance is a weak metric, because a tree asks "is this
fragment present" rather than "how similar are these vectors". So the question is open on the
feature side and has to be asked separately: do bits that encode the ROLE of an atom beat bits
that encode its element?

FCFP is the same Morgan algorithm with feature invariants - atoms typed as donor, acceptor,
aromatic, halogen, basic or acidic rather than by element. Two molecules with different skeletons
and the same arrangement of functional groups share bits here and share none in ECFP. For a
binding site that is the relevant sense of "the same".

One argument to the generator, the same vector length, the same learner, the same folds. Three
arms: the current ECFP set, the same with ECFP swapped for FCFP, and both together, since they
are not mutually exclusive and the union may beat either.

Scored as item 80 requires. The prediction: if the metric result carries over, FCFP beats ECFP on
rank; if item 96's distinction is the whole story, they are equal and the union wins slightly by
having more columns.

Reads data/feats.npz and data/rows.csv. Writes results/preds/oof_fcfp.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/oof_fcfp.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    DM = np.hstack([z["DESC"], z["MECH"]])

    cache = _pl.Path(D + "fcfp.npz")
    if cache.exists():
        FC = np.load(cache)["FC"]
    else:
        print("считаю FCFP", flush=True)
        t0 = time.time()
        gen = rdFingerprintGenerator.GetMorganGenerator(
            radius=2, fpSize=2048,
            atomInvariantsGenerator=rdFingerprintGenerator.GetMorganFeatureAtomInvGen())
        FC = np.array([gen.GetCountFingerprintAsNumPy(Chem.MolFromSmiles(s))
                       for s in rows.SMILES], dtype=np.float32)
        np.savez_compressed(cache, FC=FC)
        print(f"  {FC.shape} за {time.time()-t0:.0f} с", flush=True)

    SETS = {"ECFP (как сейчас)": np.hstack([z["FP"], DM]),
            "FCFP вместо ECFP": np.hstack([FC, DM]),
            "оба": np.hstack([z["FP"], FC, DM])}
    print("наборы: " + ", ".join(f"{k} ({v.shape[1]})" for k, v in SETS.items()) + "\n")

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for nm, X in SETS.items():
            t0 = time.time()
            r = {"seed": seed, "набор": nm}
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                Xi, fi = X[m], fold[m]
                p = np.zeros_like(y)
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    p[te] = HistGradientBoostingRegressor(**KW).fit(Xi[trn], y[trn]).predict(Xi[te])
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{nm}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            r["MACRO пара"] = round(float(np.mean([r[f"{c} пара"] for c in CYPS])), 4)
            r["MACRO rho"] = round(float(np.mean([r[f"{c} rho"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {nm:18s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("набор")[["MACRO пара", "MACRO rho"] +
                              [f"{c} rho" for c in CYPS]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")


if __name__ == "__main__":
    main()
