"""The acidic half of the mechanistic block, which was never built.

The asymmetry. For bases the block carries strength, not just presence: `pka_max_basic`,
`frac_prot_74`, `is_base_74`, `n_basic_sites`, and the topological distance from a basic nitrogen
to an aromatic ring. For acids it carries three counters - `n_acid`, `n_tetrazole`,
`n_sulfonamide` - and nothing about how acidic they are.

That is the wrong way round for one enzyme. CYP2C9 binds anions through Arg108, so acidity plays
the role there that basicity plays on CYP2D6. Item 85 found that the block does work on CYP2C9
(+0.0138 of rank) while pooling does not, which says the block is already catching something real
there with a blunt instrument.

The same dimorphite table that supplies the basic pKa supplies the acidic one - `feats.py` already
loads it and already knows which sites are bases, so everything not in that set is an acid. Five
features, mirroring the basic five: the lowest acidic pKa, the deprotonated fraction at 7.4, a
flag, a site count, and the distance from the acid to an aromatic ring.

The prediction, per enzyme and therefore falsifiable: a gain on CYP2C9, nothing on CYP2D6 or
CYP1A2. A uniform gain across all four would mean five extra columns helped, not acidity.

Reads data/feats.npz, data/rows.csv. Writes data/acid.npz and results/preds/oof_acid.json.
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
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from feats import PKA, BASIC
from cypsplit import butina_folds
from shrinkchoice import fit_apply

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
NAMES = ["pka_min_acid", "frac_deprot_74", "is_acid_74", "n_acid_sites",
         "topo_acid_to_arom_min", "topo_acid_to_arom_mean"]


def acid_block(smi):
    """Кислотные SMARTS у dimorphite требуют ЯВНОГО водорода: кислота его теряет, и шаблон
    записан как ...-[O]-[#1]. Основные шаблоны его не требуют, поэтому существующий код в
    feats.py работает без AddHs, а этот без него дал бы блок из одних нулей."""
    m0 = Chem.MolFromSmiles(smi)
    if m0 is None:
        return [0.0] * len(NAMES)
    m = Chem.AddHs(m0)
    # Шаблоны dimorphite перекрываются --- Thioic_acid (pKa 0.68) и Alcohol (14.8) оба
    # ложно матчатся на уксусной кислоте вместе с верным Carboxyl (3.46). Сам dimorphite
    # разрешает это порядком применения; здесь --- по специфичности: за атомом остаётся
    # pKa самого специфичного шаблона, где специфичность есть число атомов в запросе.
    site = {}
    for name, patt, pkas in PKA:
        if name in BASIC:
            continue
        spec = patt.GetNumAtoms()
        for h in m.GetSubstructMatches(patt):
            a = h[0]
            if a not in site or spec > site[a][0]:
                site[a] = (spec, min(pkas))
    if not site:
        return [0.0] * len(NAMES)
    idxs = list(site)
    best = min(v[1] for v in site.values())
    # Доля депротонированных для КИСЛОТЫ: зеркало формулы для основания.
    frac = 1.0 / (1.0 + 10.0 ** (best - 7.4))
    dm = Chem.GetDistanceMatrix(m)
    rings = [r for r in m.GetRingInfo().AtomRings()
             if all(m.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    ds = [dm[i, j] for i in set(idxs) for r in rings for j in r]
    return [best, frac, 1.0 if best < 7.4 else 0.0, float(len(set(idxs))),
            float(min(ds)) if ds else 0.0, float(np.mean(ds)) if ds else 0.0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/oof_acid.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    base = np.hstack([z["FP"], z["DESC"], z["MECH"]])

    cache = _pl.Path(D + "acid.npz")
    if cache.exists():
        A = np.load(cache)["A"]
    else:
        print("считаю кислотный блок", flush=True)
        t0 = time.time()
        A = np.array([acid_block(s) for s in rows.SMILES], np.float32)
        np.savez_compressed(cache, A=A)
        print(f"  {A.shape} за {time.time()-t0:.0f} с", flush=True)
    got = A[:, 2] > 0.5
    print(f"анионных при 7.4: {int(got.sum())} из {len(A)} ({got.mean():.1%}); "
          f"медиана самого кислого pKa у них {np.median(A[got, 0]):.2f}")

    SETS = {"без кислот": base, "с кислотами": np.hstack([base, A])}
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
            print(f"  сид {seed} {nm:13s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
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
