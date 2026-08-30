"""Are some of our strongest inhibitors not inhibitors at all but aggregates?

The chemistry. Compounds that appear to inhibit everything at 50 micromolar frequently are not
inhibiting: they aggregate or fall out of solution and sequester the enzyme. The signature is
known and it is structural - high lipophilicity, several aromatic rings, no ionisable group to
keep them dissolved - plus a behavioural tell, activity on all four enzymes at once.

Why it matters here rather than in general. Item 86 established that the compounds the organisers
declined to run curves on are, on CYP3A4, precisely the promiscuous multi-enzyme hits. If some
fraction of those are a solubility artefact rather than chemistry, they teach the model that
"greasy and flat means inhibitor" when the truth is "greasy and flat means precipitated", and
they teach it on all four enzymes simultaneously.

That last word is what makes the hypothesis testable rather than a story. A solubility artefact
is a property of the compound, not of the enzyme, so removing its influence should raise rank on
**all four enzymes at once**. If the gain lands on one enzyme it is not solubility; it is
something enzymatic that happens to correlate with being greasy.

An asymmetry the design has to respect. The promiscuity tell comes from the screening file, which
exists for training compounds and not for the blinded test. So it may be used as a WEIGHT, which
only touches training, and may not be used as a FEATURE, which would have to be computed for a
test molecule that has no screen. The structural part can be either.

Four arms:

  base            unchanged;
  flag            the structural aggregator score as one extra column;
  down-weight     the top decile by structural score trained at quarter weight;
  down-weight+    the same by structural score combined with screening promiscuity.

Scored as item 80 requires, per enzyme, since the prediction is about all four at once.

Reads data/feats.npz, the screening file. Writes results/preds/oof_agg.json.
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
LOW = 0.25       # вес подозреваемых
FRAC = 0.10      # какая доля попадает под подозрение


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--out", default=RES + "preds/oof_agg.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    dn = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
    mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    g = lambda blk, nm, names: blk[:, names.index(nm)]
    logp = g(z["DESC"], "MolLogP", dn)
    arom = g(z["DESC"], "NumAromaticRings", dn)
    base_f = g(z["MECH"], "is_base_74", mn)
    anion = g(z["MECH"], "ph74_n_anion", mn)
    neutral = (base_f < 0.5) & (anion < 0.5)

    zs = lambda v: (v - v.mean()) / (v.std() + 1e-9)
    struct = zs(logp) + zs(arom) + 1.5 * neutral.astype(float)

    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    sig = (sc.assign(hit=sc.log2fc_fdr < 0.05)
             .pivot_table(index="Molecule_Name", columns="enzyme", values="hit", aggfunc="first"))
    prom = sig.reindex(rows.Molecule_Name).sum(1).fillna(0).to_numpy(float)
    struct_p = struct + zs(prom)

    thr_s = np.quantile(struct, 1 - FRAC)
    thr_p = np.quantile(struct_p, 1 - FRAC)
    flag_s, flag_p = struct >= thr_s, struct_p >= thr_p
    print(f"подозреваемых по структуре {int(flag_s.sum())}, "
          f"со скринингом {int(flag_p.sum())}, пересечение {int((flag_s & flag_p).sum())}")
    print(f"у подозреваемых по структуре: logP {logp[flag_s].mean():.2f} против "
          f"{logp.mean():.2f}, колец {arom[flag_s].mean():.2f} против {arom.mean():.2f}, "
          f"ферментов со значимым скрином {prom[flag_s].mean():.2f} против {prom.mean():.2f}")
    for e, c in enumerate(CYPS):
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy()
        print(f"  {c}: медиана pIC50 у подозреваемых {np.median(y[flag_s[m]]):.2f} "
              f"против {np.median(y):.2f} по ферменту")
    print()

    XF = np.hstack([X, struct[:, None].astype(np.float32)])
    ARMS = {"база": (X, None), "флаг": (XF, None),
            "вес по структуре": (X, flag_s), "вес со скринингом": (X, flag_p)}

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for nm, (XX, fl) in ARMS.items():
            t0 = time.time()
            r = {"seed": seed, "рука": nm}
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                Xi, fi = XX[m], fold[m]
                w = None if fl is None else np.where(fl[m], LOW, 1.0)
                p = np.zeros_like(y)
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    kw = {} if w is None else {"sample_weight": w[trn]}
                    p[te] = HistGradientBoostingRegressor(**KW).fit(
                        Xi[trn], y[trn], **kw).predict(Xi[te])
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
    print(df.groupby("рука")[["MACRO пара", "MACRO rho"] +
                             [f"{c} rho" for c in CYPS]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Смотреть на ЗНАК прироста ранга по всем четырём ферментам, а не на макро.

Растворимость --- свойство соединения, не фермента. Значит если подозреваемые действительно
артефакт, понижение их веса обязано поднять ранг ВЕЗДЕ. Прирост на одном ферменте при потерях
на других означает, что найдена не растворимость, а что-то ферментное, коррелирующее с
жирностью, --- и тогда трогать веса нельзя, потому что мы выбросим настоящий сигнал.""")


if __name__ == "__main__":
    main()
