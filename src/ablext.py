"""One more row of the ablation table, with external labels added to the training folds.

The reason this is worth doing and the encoder was not enough. Swapping the representation
moved macro by 0.010 in the wrong direction, so the bottleneck is not how molecules are
described. It is more likely how many labelled ones there are: R^2 by enzyme is 0.59 / 0.37 /
0.17 / 0.23 for 3A4 / 2C9 / 2D6 / 1A2, and no amount of post-processing repairs an R^2 of
0.17. CYP2D6 has 1493 labels and the worst fit of the four, and it is the enzyme the document
is built around.

Where the data comes from. The organisers published a CheMeleon baseline for exactly these
four enzymes, and shipped its training set with it: 8068 compounds with pIC50 curated from
ChEMBL, Apache-2.0, columns named per enzyme. It roughly triples the labels on CYP2C9 and
CYP2D6 and doubles them on CYP3A4. The competition rules permit public data without
restriction.

Two checks had to pass before any of it could be used, and both are re-run here rather than
asserted. Overlap with the BLINDED TEST is zero of 750 - there is no leakage of any kind.
Overlap with our own training set is 64 of 4905, and those compounds are dropped so that no
molecule appears with two labels.

The part that decides the design. The two label sets are not on the same scale. On the 64
shared compounds the external values run higher by +0.22 / +0.35 / +0.60 / +0.56, and the
marginal gap is larger still - up to +1.36 on CYP3A4 - because ChEMBL is enriched in actives,
which people publish. That is a shift of the same size as the test-set shift this document
spent a month measuring. Appending rows without correcting it would inject exactly the bias
we have been chasing.

So three arms are run, and the middle one is the point:

  raw       - external rows appended as they are, no correction. The wrong thing to do,
              included because its size shows what the correction is worth;
  paired    - each external label shifted by the offset measured on the shared compounds.
              An independent estimate, but a thin one: n = 15 / 6 / 41 / 9;
  marginal  - shifted so that the external mean matches ours per enzyme. Robust to small n
              and wrong if the external set is genuinely more active rather than differently
              calibrated, which it probably partly is.

The truth is between the last two, and reporting both is more honest than choosing.

Held-out folds contain only OUR compounds throughout: external rows go into the training
folds and never into the fold being scored, so the metric is the one the rest of the document
uses. Same learner, same settings, same Butina folds as src/ablate.py.

Reads data/feats.npz, data/rows.csv and the external CSVs. Writes results/preds/oof_ext.json.
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
from rdkit import Chem
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

import feats as F
from cypsplit import butina_folds

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)


def canon(s):
    try:
        m = Chem.MolFromSmiles(s)
        return Chem.MolToSmiles(m) if m else None
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ext-x", required=True)
    ap.add_argument("--ext-y", required=True)
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--out", default=RES + "preds/oof_ext.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])

    ext = pd.concat([pd.read_csv(a.ext_x), pd.read_csv(a.ext_y)], axis=1)
    ext["k"] = [canon(s) for s in ext.OPENADMET_CANONICAL_SMILES]
    ours = set(filter(None, (canon(s) for s in rows.SMILES)))
    te_keys = set(filter(None, (canon(s) for s in
                                pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv").SMILES)))

    leak = ext.k.isin(te_keys).sum()
    print(f"пересечение внешних данных с ТЕСТОМ: {leak}")
    if leak:
        raise SystemExit("внешние данные пересекаются с тестом, использовать нельзя")
    dup = ext.k.isin(ours)
    print(f"пересечение с обучением: {int(dup.sum())}, эти строки отброшены")

    # Смещения по парным соединениям считаются ДО того, как дубликаты выброшены.
    shared = ext[dup].merge(
        tr.assign(k=[canon(s) for s in rows.SMILES]), on="k", how="inner")
    off_pair, off_marg = {}, {}
    print(f"\n{'фермент':8s} {'парных':>7s} {'сдвиг парн.':>12s} {'сдвиг маргин.':>14s}")
    for c in CYPS:
        xc, yc = f"{c}_pIC50_direct_inhibition", f"OPENADMET_LOGAC50_{c.lower()}"
        ok = shared[xc].notna() & shared[yc].notna()
        off_pair[c] = float((shared[yc][ok] - shared[xc][ok]).mean()) if ok.sum() >= 5 else np.nan
        off_marg[c] = float(ext[yc].dropna().mean() - tr[xc].dropna().mean())
        print(f"{c:8s} {int(ok.sum()):7d} {off_pair[c]:+12.2f} {off_marg[c]:+14.2f}")

    ext = ext[~dup].reset_index(drop=True)
    print(f"\nвнешних строк после отбрасывания дубликатов: {len(ext)}")

    print("считаю признаки для внешних соединений", flush=True)
    dn = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
    mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    FP, dsc, M, ok_idx = F.build(list(ext.OPENADMET_CANONICAL_SMILES), dn, mn)
    Xe = np.hstack([FP, dsc.to_numpy(np.float32), M.to_numpy(np.float32)])
    ext = ext.iloc[ok_idx].reset_index(drop=True)
    print(f"  разобралось {len(ok_idx)} из {len(ok_idx) + (len(Xe) - len(ok_idx))}, "
          f"матрица {Xe.shape}", flush=True)
    if Xe.shape[1] != X.shape[1]:
        raise SystemExit(f"ширина не совпала: {Xe.shape[1]} против {X.shape[1]}")

    ARMS = {"без внешних": None, "внешние как есть": {c: 0.0 for c in CYPS},
            "сдвиг по парным": off_pair, "сдвиг по маргинали": off_marg}

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for name, off in ARMS.items():
            t0 = time.time()
            r = {"seed": seed, "рука": name}
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                Xi, fi = X[m], fold[m]
                ye = ext[f"OPENADMET_LOGAC50_{c.lower()}"].to_numpy(float)
                em = ~np.isnan(ye)
                p = np.zeros_like(y)
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    if off is None:
                        Xt, yt = Xi[trn], y[trn]
                    else:
                        Xt = np.vstack([Xi[trn], Xe[em]])
                        yt = np.concatenate([y[trn], ye[em] - off[c]])
                    p[te] = HistGradientBoostingRegressor(**KW).fit(Xt, yt).predict(Xi[te])
                out[f"{seed}|{name}|{c}"] = p.tolist()
                r[c] = round(float(strae(y, p, y_true_upper=hi, y_true_lower=lo)), 4)
            r["MACRO"] = round(float(np.mean([r[c] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {name:20s} макро {r['MACRO']:.4f}  ({time.time()-t0:.0f} с)",
                  flush=True)

    print()
    print(pd.DataFrame(table)[["seed", "рука", *CYPS, "MACRO"]].to_string(index=False))
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Первая строка --- контроль: она обязана воспроизвести 0.7673 из документа, иначе сравнивать
остальные не с чем. Вторая показывает, во что обходится игнорировать разницу шкал. Третья и
четвёртая --- две честные попытки её учесть, и расхождение между ними есть оставшаяся
неопределённость, а не выбор в пользу лучшей.""")


if __name__ == "__main__":
    main()
