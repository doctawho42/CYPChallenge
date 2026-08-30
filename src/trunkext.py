"""The shared trunk again, with external CYP labels in the auxiliary head instead of the screen.

Why this exists. The joint-likelihood trunk was built, measured and shelved: the screening
channel it was given carries little about pIC50, and the noise ladder showed that destroying
that channel costs almost nothing on the free-head arm. The machinery was never the problem -
it was handed a channel with no information in it.

The external CYP labels are a channel with information in it. 8004 compounds with pIC50 on the
same four enzymes, curated from ChEMBL by the organisers, zero overlap with the blinded test.
Appending them as extra training rows to the boosting is worth -0.013 macro (src/ablext.py),
and the largest share of that lands on CYP2D6, the enzyme with the worst R^2.

But appending forced a choice nobody could make well. The two label sets differ in scale by
+0.22 to +0.60 on shared compounds, and all three corrections tried - none, paired offset,
marginal offset - gave different answers, with the crudest actively harmful. An auxiliary head
removes the choice rather than answering it: each source gets its own output, the trunk is
shared, and nothing has to be assumed about how the two scales relate.

The construction needs no change to src/trunk.py, which matters because that file's numbers
are published. The external rows are appended with fold index -1, so `fold == f` is never true
for them and they sit in every training fold and no held-out one. Their pIC50 entries are NaN,
so the primary head gets no gradient from them; our own rows have NaN in the auxiliary block,
so the auxiliary head gets no gradient from us. run_fold is imported and used exactly as it
stands.

Mode is twohead deliberately. The calibrated mode routes the auxiliary target through the
instrument's Hill curve, which describes the screening readout and has nothing to do with a
pIC50 measured elsewhere.

What the arms mean, after the first attempt at a control turned out to be wrong. It looked
obvious that lambda_scr = 0 here must reproduce trunk_twohead.json exactly, since external
rows carry no pIC50 label and so contribute no gradient to the primary head. The check was
written to assert it and the run refuted it: maximum absolute difference 6.49, not zero.

The reason is in run_fold, and it is not a leak in the fold logic. Feature standardisation is
computed over the training rows, and there are now 12909 of them instead of 4905, so the input
scaling changes. So does the number of gradient steps: epochs are fixed, batches are not, so
more rows means more updates. External rows change the model without contributing a single
gradient through either head.

That makes the effect decompose into two parts rather than one, which is more informative than
the control it replaced:

  trunk_twohead at lambda 0   - no external rows at all, the published baseline;
  this file at lambda 0       - external rows present, auxiliary term off. The difference from
                                the line above is what the rows do through normalisation and
                                step count alone, with no information transfer;
  this file at lambda 3       - the auxiliary gradient on top. The difference from the line
                                above is what the channel actually carries.

Only the third minus the second is the auxiliary head's contribution. Comparing the third to
the published baseline would credit the channel with an effect that is partly bookkeeping.

Reads data/feats.npz, the external CSVs, results/preds/trunk_twohead.json for the control
comparison. Writes results/preds/trunk_ext.json.
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

import feats as F
import trunk as T
from cypsplit import butina_folds

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]


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
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--lams", default="0,3.0")
    ap.add_argument("--device", default="mps")
    ap.add_argument("--out", default=RES + "preds/trunk_ext.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    lams = [float(x) for x in a.lams.split(",")]

    X, y, lo, hi, _, smiles = T.load(T.BLOCKS)
    rows = pd.read_csv(D + "rows.csv")

    ext = pd.concat([pd.read_csv(a.ext_x), pd.read_csv(a.ext_y)], axis=1)
    ext["k"] = [canon(s) for s in ext.OPENADMET_CANONICAL_SMILES]
    ours = set(filter(None, (canon(s) for s in rows.SMILES)))
    te_keys = set(filter(None, (canon(s) for s in
                                pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv").SMILES)))
    if ext.k.isin(te_keys).any():
        raise SystemExit("внешние данные пересекаются с тестом")
    ext = ext[~ext.k.isin(ours)].reset_index(drop=True)

    dn = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
    mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    print(f"признаки для {len(ext)} внешних соединений", flush=True)
    FP, dsc, M, ok = F.build(list(ext.OPENADMET_CANONICAL_SMILES), dn, mn)
    ext = ext.iloc[ok].reset_index(drop=True)
    parts = {"FP": np.log1p(FP), "DESC": dsc.to_numpy(np.float32), "MECH": M.to_numpy(np.float32)}
    Xe = np.hstack([parts[b] for b in T.BLOCKS.split("+")]).astype(np.float32)
    Xe = np.nan_to_num(Xe, nan=0.0, posinf=0.0, neginf=0.0)
    if Xe.shape[1] != X.shape[1]:
        raise SystemExit(f"ширина не совпала: {Xe.shape[1]} против {X.shape[1]}")

    ye = np.stack([ext[f"OPENADMET_LOGAC50_{c.lower()}"].to_numpy(np.float32) for c in CYPS], 1)
    n, m = len(X), len(Xe)
    XX = np.vstack([X, Xe])
    # Наши строки: метка pIC50 есть, вспомогательной нет. Внешние: наоборот.
    yy = np.vstack([y, np.full_like(ye, np.nan)])
    ss = np.vstack([np.full_like(y, np.nan), ye])
    print(f"строк: наших {n}, внешних {m}, всего {len(XX)}")
    print(f"вспомогательных меток по ферментам: "
          + ", ".join(f"{c} {int(np.isfinite(ye[:, e]).sum())}" for e, c in enumerate(CYPS)))

    base = None
    bp = _pl.Path(RES + "preds/trunk_twohead.json")
    if bp.exists():
        base = json.load(open(bp))["preds"]

    table, saved = [], {}
    for seed in seeds:
        fold, _ = butina_folds(smiles, seed=seed)
        # -1 значит «всегда в обучении»: fold == f для такой строки не выполняется никогда.
        ff = np.concatenate([fold, np.full(m, -1, dtype=fold.dtype)])
        for lam in lams:
            t0 = time.time()
            pred = np.full_like(yy, np.nan)
            for f in range(5):
                te = ff == f
                if te.sum() == 0:
                    continue
                pred[te] = T.run_fold(XX, yy, ss, ff, f, lam, seed, a.device, "twohead")
            p = pred[:n]
            r = T.evaluate(y, lo, hi, p)
            r["seed"], r["lambda"] = seed, lam
            table.append(r)
            saved[f"{seed}|{lam}"] = np.where(np.isnan(y), np.nan, p).tolist()
            note = ""
            if lam == 0 and base is not None and f"twohead|{seed}|0.0" in base:
                bt = pd.DataFrame(json.load(open(bp))["table"])
                row = bt[(bt.seed == seed) & (bt["lambda"] == 0.0)]
                if len(row):
                    note = (f"   без внешних строк было {float(row.MACRO.iloc[0]):.4f}, "
                            f"то есть сами строки стоят {r['MACRO'] - float(row.MACRO.iloc[0]):+.4f}")
            print(f"  сид {seed} lambda {lam:<4} макро {r['MACRO']:.4f}  "
                  f"({time.time()-t0:.0f} с){note}", flush=True)

    df = pd.DataFrame(table)[["seed", "lambda", *CYPS, "MACRO"]]
    print()
    print(df.to_string(index=False))
    json.dump({"table": table, "preds": saved}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Эффект раскладывается надвое, и путать половины нельзя.

Первая --- что делают сами внешние СТРОКИ при выключенном вспомогательном слагаемом. Они
не дают ни одного градиента, но меняют стандартизацию входа и число шагов, потому что эпох
фиксированное число, а батчей стало больше. Это чистая бухгалтерия, и она печатается рядом
с каждой строкой lambda = 0.

Вторая --- что даёт сам КАНАЛ: разность между lambda = 3 и lambda = 0 в этом файле, при
одних и тех же строках. Только она и есть вклад вспомогательной головы.

Сравнивать вторую надо не с нулём, а с -0.013: во столько обходится то же самое, поданное
дописыванием строк в бустинг. Вспомогательная голова оправдана, только если снимает разницу
шкал не хуже, чем её снимает отсутствие всякой поправки.""")


if __name__ == "__main__":
    main()
