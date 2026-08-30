"""Does three-dimensional shape add anything the two-dimensional features do not?

The gap it addresses. Every feature in the matrix is 2D, and the mechanistic block is about
basicity, which is to say about CYP2D6. Shape is the other half of isoform selectivity and is
encoded nowhere: CYP1A2 is a narrow planar slot, CYP3A4 a large flexible cavity, and nothing in
2295 columns separates a flat polyaromatic from a globular molecule of the same weight and logP.
src/shape3d.py adds sixteen columns from one ETKDG conformer.

The prediction, recorded before the run and differing by enzyme, which is what makes it useful:
planarity should matter on CYP1A2, volume and flexibility on CYP3A4, and neither much on CYP2D6,
where charge decides and where item 82 already showed the working features are the geometric ones
in the mechanistic block.

One of the sixteen is not a shape index. The mechanistic block measures the distance from a basic
nitrogen to an aromatic ring in BONDS while the CYP2D6 pharmacophore is defined in ANGSTROMS,
5 to 7, so the conformer supplies the real distance. If CYP2D6 moves at all, that column is the
first place to look.

Scored as item 80 requires: raw, after the affine pair, and as a change in rank correlation. The
first is reported only to show how much it misleads.

Reads data/feats.npz and data/shape3d.npz. Writes results/preds/oof_shape.json.
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
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/oof_shape.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    S = np.load(D + "shape3d.npz")["train"]
    if len(S) != len(rows):
        raise SystemExit(f"форма {len(S)} строк против {len(rows)}")
    base = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    SETS = {"без формы": base, "с формой": np.hstack([base, S])}
    print(f"база {base.shape[1]} колонок, с формой {SETS['с формой'].shape[1]}\n")

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for name, X in SETS.items():
            t0 = time.time()
            r = {"seed": seed, "набор": name}
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
                out[f"{seed}|{name}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            r["MACRO пара"] = round(float(np.mean([r[f"{c} пара"] for c in CYPS])), 4)
            r["MACRO rho"] = round(float(np.mean([r[f"{c} rho"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {name:10s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("набор")[["MACRO пара", "MACRO rho"] +
                              [f"{c} rho" for c in CYPS]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Смотреть на прирост ранга по ферментам, а не на макро: предсказание было
поферментным и разнонаправленным, и макро его усреднит в кашу. Планарность на 1A2, объём и
гибкость на 3A4, около нуля на 2D6 --- если так, предсказание сбылось; если прирост окажется
равномерным, значит работает не форма как механизм, а просто шестнадцать лишних колонок.""")


if __name__ == "__main__":
    main()
