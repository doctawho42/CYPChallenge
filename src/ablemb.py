"""One more row of the ablation table, with the representation swapped for a pretrained one.

The question this answers, and the reason it is worth an hour of compute. Every row of
src/ablate.py is Morgan counts, RDKit descriptors and the mechanistic block in some
combination. When the document concluded that the choice of model is worth 0.005 while the
post-processing is worth 0.051, that conclusion was drawn entirely inside that one
representation - and written as though it held generally. If a different representation moves
the number, the reading changes: the search was local, and we mistook it for global.

Nothing is changed except the features. Same learner and the same settings as ablate.py, the
same Butina folds, the same masks, the same metric. The embedding comes from
src/embed.py, which runs the organisers' pretrained chemprop encoder in a separate
environment and leaves an array behind.

Three rows are computed and two of them are controls:

  EMB          - the embedding alone, against FP+DESC alone;
  EMB+MECH     - with the mechanistic block, against FP+DESC+MECH, the headline row;
  FP+DESC+MECH - recomputed here rather than read from oof.json, so that any difference in
                 environment or library version shows up as a discrepancy in a number we
                 already know rather than as a silent bias in the new one.

That third row is the one that makes the comparison trustworthy. If it does not reproduce
0.767 the rest of the table means nothing.

Reads data/feats.npz, data/emb_chemprop_rdkit2d.npz, results/preds/oof.json. Writes
results/preds/oof_emb.json.
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
# Точно те же настройки, что в src/ablate.py. Меняется представление, не ученик.
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--emb", default=D + "emb_chemprop_rdkit2d.npz")
    ap.add_argument("--out", default=RES + "preds/oof_emb.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    E = np.load(a.emb)["train"]
    if len(E) != len(rows):
        raise SystemExit(f"эмбеддинг {len(E)} строк против {len(rows)} в rows.csv")

    SETS = {
        "EMB":          E,
        "EMB+MECH":     np.hstack([E, z["MECH"]]),
        "FP+DESC+MECH": np.hstack([z["FP"], z["DESC"], z["MECH"]]),
    }
    print(f"эмбеддинг {E.shape}, для сравнения FP+DESC+MECH {SETS['FP+DESC+MECH'].shape}\n")

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for name, X in SETS.items():
            t0 = time.time()
            r = {"seed": seed, "set": name}
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
                out[f"{seed}|{name}|{c}"] = p.tolist()
                r[c] = round(float(strae(y, p, y_true_upper=hi, y_true_lower=lo)), 4)
            r["MACRO"] = round(float(np.mean([r[c] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {name:14s} макро {r['MACRO']:.4f}  ({time.time()-t0:.0f} с)",
                  flush=True)

    df = pd.DataFrame(table)[["seed", "set", *CYPS, "MACRO"]]
    print()
    print(df.to_string(index=False))

    base = [r for r in table if r["set"] == "FP+DESC+MECH" and r["seed"] == seeds[0]]
    if base:
        b = base[0]["MACRO"]
        print(f"\nКонтроль воспроизводимости: FP+DESC+MECH здесь {b:.4f}, "
              f"в документе 0.7673, расхождение {abs(b - 0.7673):.4f}.")
        print("Если оно больше пары тысячных, сравнивать строки нельзя и надо разбираться.")

    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать результат. Если эмбеддинг не двигает макро за пределы разброса по сидам,
это ОТРИЦАТЕЛЬНЫЙ результат, и он ценен: значит вывод про 0.005 против 0.051 переносится
за пределы одного представления, а не только внутри него. Если двигает - значит поиск был
локальным, и половину выводов раздела про постобработку надо перечитывать с этой поправкой.""")


if __name__ == "__main__":
    main()
