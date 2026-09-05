"""What is the electrostatic channel worth, when a rough version of it is already in the matrix?

An outside proposal is to compute DFT-quality ESP/RESP charges as descriptors. The premise this
file tests is that the proposal is not "add electrostatics" but "refine electrostatics from
Gasteiger to DFT", because a coarse version is already present:

    20 заряд-производных в DESC   MaxPartialCharge, MinPartialCharge, MaxAbsPartialCharge,
                                  MinAbsPartialCharge, BCUT2D_CHGHI/CHGLO, PEOE_VSA1..14
                                  -- PEOE_VSA is molecular surface area binned by Gasteiger
                                  partial charge, i.e. a coarse ESP field folded into a histogram
    21 EState в DESC              EState_VSA1..11, VSA_EState1..10
    10 признаков мостика в MECH   pka_max_basic, frac_prot_74, is_base_74, ph74_net_charge,
                                  ph74_n_cation/anion, topo_bN_to_arom_min/mean, pharm_2d6,
                                  pharm_2d6_x_prot
     3 координации гема в MECH    n_pyridineN, n_imidazoleN, n_tetrazole

Permutation of a block at PREDICTION time bounds what the block currently delivers, and therefore
bounds what refining the same quantity can deliver: a refinement can at most replace the coarse
version with a perfect one, and the coarse version's whole contribution is its permutation cost.
Item 128's procedure, item 227's ladder -- measure the derivative on the cheap version first.

**The control this needs and an outside run of it did not have.** A 20-column block compared
against a 216-column block is not a comparison: bigger blocks cost more because they are bigger.
So every named block is measured against a SIZE-MATCHED random block of the same width drawn from
the same source, ten draws, and what is reported is the named block's cost ABOVE its null. Without
that line the charge block's 0.022 on CYP2D6 cannot be told from any twenty columns.

Spearman is the primary metric: it is the one the affine pair preserves (item 77), so a permutation
cost in rank is a cost the submission would actually see. RMSE is printed beside it and is NOT the
submission's metric.

Reads data/feats.npz, data/rows.csv, data/desc_names.csv, data/mech_names.csv.
Five folds, four seeds, one HistGB fit per fold -- the arms are permutations of its input, not
refits, so the whole sweep costs 80 fits.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES

import argparse
import json
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor

from cypsplit import butina_folds

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
REPS = 5           # повторов перестановки на плечо
NULLS = 10         # случайных блоков той же ширины


def blocks(desc_names, mech_names, n_fp):
    d = {n: n_fp + i for i, n in enumerate(desc_names)}
    m = {n: n_fp + len(desc_names) + i for i, n in enumerate(mech_names)}
    chg = [d[n] for n in desc_names
           if any(k in n for k in ("PartialCharge", "BCUT2D_CHG", "PEOE_VSA"))]
    est = [d[n] for n in desc_names
           if n.startswith("EState_VSA") or n.startswith("VSA_EState")]
    salt = [m[n] for n in mech_names
            if any(k in n for k in ("pka", "prot_74", "is_base", "ph74", "topo_bN", "pharm_2d6"))]
    heme = [m[n] for n in mech_names if n in ("n_pyridineN", "n_imidazoleN", "n_tetrazole")]
    return {
        "заряд (20)": chg,
        "EState (21)": est,
        "мостик MECH (10)": salt,
        "гем MECH (3)": heme,
        "весь DESC (217)": list(d.values()),
        "весь MECH (30)": list(m.values()),
        "весь FP (2048)": list(range(n_fp)),
    }, list(d.values()), list(m.values())


def permuted_pred(model, Xte, cols, rng, reps=REPS):
    out = []
    for _ in range(reps):
        Z = Xte.copy()
        Z[:, cols] = Z[rng.permutation(len(Z))][:, cols]
        out.append(model.predict(Z))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/chgperm.json")
    a = ap.parse_args()

    z = np.load(D + "feats.npz")
    FP, DESC, MECH = z["FP"], z["DESC"], z["MECH"]
    X = np.hstack([FP, DESC, MECH])
    rows = pd.read_csv(D + "rows.csv")
    desc_names = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
    mech_names = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    BL, desc_cols, mech_cols = blocks(desc_names, mech_names, FP.shape[1])
    print("блоки: " + ", ".join(f"{k}={len(v)}" for k, v in BL.items()), flush=True)

    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)

    recs = []
    for seed in [int(s) for s in a.seeds.split(",")]:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        print(f"\n=== сид {seed} ===", flush=True)
        for e, c in enumerate(CYPS):
            obs = np.isfinite(Y[:, e])
            y = Y[obs, e]
            Xi, fi = X[obs], fold[obs]
            base = np.zeros(len(y))
            arms = {k: [np.zeros(len(y)) for _ in range(REPS)] for k in BL}
            nulls = {"нуль 20 из DESC": (20, desc_cols), "нуль 10 из MECH": (10, mech_cols),
                     "нуль 3 из MECH": (3, mech_cols), "нуль 21 из DESC": (21, desc_cols)}
            nl = {k: [np.zeros(len(y)) for _ in range(NULLS)] for k in nulls}
            t0 = time.time()
            for f in range(5):
                trn, te = fi != f, fi == f
                if te.sum() == 0:
                    continue
                mdl = HistGradientBoostingRegressor(**KW).fit(Xi[trn], y[trn])
                base[te] = mdl.predict(Xi[te])
                rng = np.random.default_rng(1000 * seed + f)
                for k, cols in BL.items():
                    for r, p in enumerate(permuted_pred(mdl, Xi[te], cols, rng)):
                        arms[k][r][te] = p
                for k, (w, src) in nulls.items():
                    for r in range(NULLS):
                        cols = list(rng.choice(src, w, replace=False))
                        nl[k][r][te] = permuted_pred(mdl, Xi[te], cols, rng, reps=1)[0]
            r0 = spearmanr(y, base).statistic
            e0 = float(np.sqrt(np.mean((base - y) ** 2)))
            print(f"  {c}: n {len(y)}, {time.time()-t0:.0f} с | база RMSE {e0:.4f} rho {r0:.4f}",
                  flush=True)
            recs.append({"seed": seed, "cyp": c, "arm": "база", "rho": float(r0),
                         "rmse": e0, "drho": 0.0, "n": int(len(y))})
            for k in list(BL) + list(nulls):
                ps = arms[k] if k in BL else nl[k]
                rs = [spearmanr(y, p).statistic for p in ps]
                es = [float(np.sqrt(np.mean((p - y) ** 2))) for p in ps]
                dr = r0 - float(np.mean(rs))
                print(f"      {k:18s} RMSE {np.mean(es):.4f}  rho {np.mean(rs):.4f}  "
                      f"цена ранга {dr:+.4f}  (sd по повторам {np.std(rs):.4f})", flush=True)
                recs.append({"seed": seed, "cyp": c, "arm": k, "rho": float(np.mean(rs)),
                             "rmse": float(np.mean(es)), "drho": float(dr), "n": int(len(y))})
        json.dump(recs, open(a.out, "w"))

    df = pd.DataFrame(recs)
    p = df.pivot_table(index="arm", columns="cyp", values="drho")
    print("\n\nЦЕНА РАНГА (rho базы минус rho после перестановки), среднее по 4 сидам")
    print(p.round(4).to_string())
    print("\nНАЗВАННЫЙ БЛОК СВЕРХ СВОЕГО РАЗМЕРНОГО НУЛЯ")
    for named, null in [("заряд (20)", "нуль 20 из DESC"), ("EState (21)", "нуль 21 из DESC"),
                        ("мостик MECH (10)", "нуль 10 из MECH"), ("гем MECH (3)", "нуль 3 из MECH")]:
        if named in p.index and null in p.index:
            d = p.loc[named] - p.loc[null]
            print(f"  {named:18s} " + "  ".join(f"{c} {d[c]:+.4f}" for c in CYPS))
    print(f"\nсохранено: {a.out}", flush=True)


if __name__ == "__main__":
    main()
