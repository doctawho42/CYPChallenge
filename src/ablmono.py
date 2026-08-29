"""Constrain the model where biochemistry already knows the sign.

Why this is not decoration. `monotonic_cst` costs nothing and is usually sold as a mild
regulariser, but the reason it belongs in THIS project is different: a monotone function
extrapolates safely and a free one does not. Two months of this repository say the test set
sits where the training set is thin, and in a sparse region an unconstrained tree can go any
direction it likes while a constrained one cannot. It is the only knob here that joins the
shift work to the chemistry inside a single parameter.

The sign is known for a handful of features, not for 2295, and it is known per enzyme:

  CYP1A2  narrow planar cavity - more aromatic rings help, sp3 character hurts;
  CYP2C9  Arg108 binds anions - anionic character helps;
  CYP2D6  Asp301 salt bridge to protonated amines - protonated fraction helps. This one is
          not borrowed from the literature: item 48 measured it here, bases are +0.55 more
          active on CYP2D6 while they are 0.24 to 0.45 LESS active on the other three;
  CYP3A4  large hydrophobic site - lipophilicity helps.

Three arms, and the third is the control that makes the second interpretable:

  free      - no constraints. The control, which must reproduce 0.7673;
  correct   - the signs above;
  flipped   - every sign reversed.

The flipped arm is not a curiosity. A monotonic constraint pointing against the data does not
bias the fit slightly, it ERASES the feature: fitting y = -x under a +1 constraint gives a
model whose prediction range is exactly zero. So if `correct` and `flipped` score the same, the
constrained features carry nothing and the whole exercise is inert; if `flipped` is clearly
worse, the biochemistry is doing real work and the direction is right.

Tails are reported separately, because that is where the argument lives. Macro is dominated by
the bulk, and the claim is about the thin region, so ST-RAE is also computed on the top and
bottom deciles of each enzyme's own constrained feature.

Reads data/feats.npz, data/rows.csv. Writes results/preds/oof_mono.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)

# Знак при pIC50: +1 значит "больше признака --- потенциальнее". Имя ведущего признака
# фермента стоит первым, по нему считаются хвосты.
SIGNS = {
    "CYP1A2": [("NumAromaticRings", +1), ("FractionCSP3", -1)],
    "CYP2C9": [("ph74_n_anion", +1)],
    "CYP2D6": [("frac_prot_74", +1)],
    "CYP3A4": [("MolLogP", +1)],
}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--out", default=RES + "preds/oof_mono.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    dn = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
    mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    nfp, nd = z["FP"].shape[1], len(dn)
    pos = {**{n: nfp + i for i, n in enumerate(dn)},
           **{n: nfp + nd + i for i, n in enumerate(mn)}}

    cst, lead = {}, {}
    for c in CYPS:
        v = np.zeros(X.shape[1], dtype=int)
        for name, s in SIGNS[c]:
            if name not in pos:
                raise SystemExit(f"признак {name} не найден")
            v[pos[name]] = s
        cst[c], lead[c] = v, pos[SIGNS[c][0][0]]
        print(f"{c}: " + ", ".join(f"{n} {s:+d}" for n, s in SIGNS[c]))
    print()

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for arm in ("свободно", "по знаку", "знак наоборот"):
            t0 = time.time()
            r = {"seed": seed, "рука": arm}
            tails = []
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                Xi, fi = X[m], fold[m]
                kw = {}
                if arm == "по знаку":
                    kw["monotonic_cst"] = cst[c]
                elif arm == "знак наоборот":
                    kw["monotonic_cst"] = -cst[c]
                p = np.zeros_like(y)
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    p[te] = HistGradientBoostingRegressor(**KW, **kw).fit(Xi[trn], y[trn]).predict(Xi[te])
                out[f"{seed}|{arm}|{c}"] = p.tolist()
                r[c] = round(float(strae(y, p, y_true_upper=hi, y_true_lower=lo)), 4)
                # Хвосты по ведущему признаку: там, где обучение разрежено и живёт сдвиг.
                v = Xi[:, lead[c]]
                q = np.quantile(v, [0.1, 0.9])
                t = (v <= q[0]) | (v >= q[1])
                tails.append(float(strae(y[t], p[t], y_true_upper=hi[t], y_true_lower=lo[t])))
            r["MACRO"] = round(float(np.mean([r[c] for c in CYPS])), 4)
            r["ХВОСТЫ"] = round(float(np.mean(tails)), 4)
            table.append(r)
            print(f"  сид {seed} {arm:15s} макро {r['MACRO']:.4f}  хвосты {r['ХВОСТЫ']:.4f}  "
                  f"({time.time()-t0:.0f} с)", flush=True)

    print()
    print(pd.DataFrame(table)[["seed", "рука", *CYPS, "MACRO", "ХВОСТЫ"]].to_string(index=False))
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Сначала третью строку, а не вторую. Если "знак наоборот" стоит там же, где
"свободно", ограниченные признаки моделью не используются вовсе, и разность между первыми
двумя строками --- шум, чем бы она ни оказалась. Только если перевёрнутый знак заметно хуже,
второй строке есть чему верить.

Потом столбец ХВОСТЫ отдельно от МАКРО. Довод в пользу монотонности --- про безопасную
экстраполяцию в разреженную область, а не про среднее. Выигрыш в макро при нуле на хвостах
означал бы, что сработала регуляризация, а не то, ради чего это делалось.""")


if __name__ == "__main__":
    main()
