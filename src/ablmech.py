"""Which third of the mechanistic block does the work, now that we know the block does some.

Why this exists. Two facts sat next to each other without being compared. The block's ST-RAE
gain was recorded as indistinguishable from zero, so it looked decorative; and `frac_prot_74`,
the protonated fraction on which this document's whole CYP2D6 argument rests, turns out never to
be split on at all, which looked like confirmation. Measuring the block after the affine pair
rather than before it reverses the first fact: on CYP2D6 the block is worth -0.022 with
p = 0.0004 and adds 0.045 of rank correlation, more than any other intervention here. So the
chemistry does enter the predictor - just not through the feature the text credits.

The block splits into three groups that match how the document thinks about it:

  groups     16 counts of nitrogen and acid functionality - what the molecule contains;
  state       8 features of protonation at pH 7.4 - pKa, protonated fraction, charges;
  geometry    6 topological distances from basic nitrogen to an aromatic ring, plus the
              explicit CYP2D6 pharmacophore flags.

Prediction, written before the run. `frac_prot_74` sits in `state` and is inert, and the salt
bridge is a geometric constraint rather than a counting one, so `geometry` should carry the
CYP2D6 gain and `state` should carry little. If instead `state` carries it, the inert feature is
merely redundant with its neighbours and the text is right after all. If `groups` carries it,
the effect is ordinary functional-group information and the CYP2D6 story is decoration on top of
something plainer.

Scored the way item 80 says to score: raw, after the affine pair, and as a change in rank
correlation. Reads data/feats.npz and data/rows.csv. Writes results/preds/oof_mech.json.
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
STATE = ["pka_max_basic", "frac_prot_74", "is_base_74", "n_basic_sites",
         "formal_charge", "ph74_net_charge", "ph74_n_cation", "ph74_n_anion"]
GEOM = ["topo_bN_to_arom_min", "topo_bN_to_arom_mean", "n_bN_geom", "topo_mbc_to_arom",
        "pharm_2d6", "pharm_2d6_x_prot"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--cyps", default="CYP2D6,CYP2C9")
    ap.add_argument("--out", default=RES + "preds/oof_mech.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    cyps = a.cyps.split(",")

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    for n in STATE + GEOM:
        if n not in mn:
            raise SystemExit(f"признак {n} не найден в блоке")
    idx = {"состояние": [mn.index(n) for n in STATE],
           "геометрия": [mn.index(n) for n in GEOM]}
    idx["группы"] = [i for i in range(len(mn)) if i not in set(idx["состояние"] + idx["геометрия"])]
    base = np.hstack([z["FP"], z["DESC"]])
    print("блок разбит: " + ", ".join(f"{k} {len(v)}" for k, v in idx.items()))

    SETS = {"без блока": None, "весь блок": list(range(len(mn)))}
    SETS.update({f"только {k}": v for k, v in idx.items()})

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for name, cols in SETS.items():
            t0 = time.time()
            X = base if cols is None else np.hstack([base, z["MECH"][:, cols]])
            r = {"seed": seed, "набор": name}
            for c in cyps:
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
                r[f"{c} сырой"] = round(float(strae(y, p, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            table.append(r)
            print(f"  сид {seed} {name:18s} "
                  + "  ".join(f"{c}: {r[f'{c} пара']:.4f} rho {r[f'{c} rho']:.4f}" for c in cyps)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("набор").mean(numeric_only=True).drop(columns="seed").to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Смотреть надо на столбцы «пара» и «rho», а не на «сырой»: у блока сигнал ранговый,
и сырая мера его занижает --- на CYP2D6 она даёт p = 0.09 там, где после пары выходит 0.0004.

Ответ --- та группа, чей rho ближе всего к rho всего блока. Если это «геометрия», значит
работает расстояние от основного азота до ароматики, то есть солевой мостик в той форме, в
какой его вообще можно закодировать признаком. Если «состояние» --- значит инертность
frac_prot_74 означала лишь избыточность внутри группы, а не отсутствие вклада.""")


if __name__ == "__main__":
    main()
