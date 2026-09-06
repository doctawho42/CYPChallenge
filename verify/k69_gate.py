"""How much of the gate do we already recover? The fork item 234 opened and nobody has measured.

Item 234 folded the label into a conjunction:

    is_TDI  <=>  (Delta > log10 2)  AND  (pi_TDI > 4.301)

and measured the second conjunct alone, as an oracle, at MCC +0.5667 on CYP3A4 -- larger than the
whole deployed classifier's 0.315. The gate is a threshold on pi_TDI, a quantity we predict
directly with a regressor at RMSE 0.75. Nobody has ever asked what our PREDICTED gate scores.

That question forks the remaining effort and nothing else here does:

  * recover most of 0.567 and the gate is solved, every remaining hour belongs to Delta;
  * recover half and there is better than 0.1 of MCC sitting in a better potency model on the
    gate-closed population -- far cheaper than a better Delta model, which item 229 measured at
    R^2 near zero.

Item 229's decomposition asked the same kind of question in the WRONG coordinates -- level and
shift -- because the piecewise form of the rule suggested them. In conjunction coordinates the
decomposition is gate and shift, and the mixed oracles below are the ones that correspond to
something the rule actually does.

**One thing here is not a detail.** The gate is `pi_TDI > 4.301`, but thresholding a SHRUNK
prediction at the true threshold is not the right rule: a regressor pulled toward the mean crosses
4.301 at the wrong place. So the cut on pi_hat is also fitted out of fold, and the difference
between the two is reported. If it is large, that is a free correction that costs one number.

Reads data/feats.npz, data/rows.csv, both label files. CYP3A4 and CYP2D6 only. Same folds and
learner as src/ablate.py and verify/k67_tdiceiling.py, so the numbers compose with item 229's.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES

import argparse
import json
import time

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import matthews_corrcoef, roc_auc_score

from cypsplit import butina_folds

TDI_CYPS = ["CYP3A4", "CYP2D6"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
L = np.log10(2.0)
GATE = 4.0 + L


def oof_reg(X, y, fold, tag=""):
    p = np.full(len(y), np.nan)
    obs = np.isfinite(y)
    for f in range(5):
        trn, te = obs & (fold != f), fold == f
        if te.sum() == 0 or trn.sum() < 50:
            continue
        t0 = time.time()
        p[te] = HistGradientBoostingRegressor(**KW).fit(X[trn], y[trn]).predict(X[te])
        print(f"      {tag} фолд {f}: {int(trn.sum())} обуч, {time.time()-t0:.0f} с", flush=True)
    return p


def fitted_cut(pred, truth_bool, fold):
    """Порог на предсказании, подобранный ВНЕ ФОЛДА по эмпирическому MCC."""
    grid = np.quantile(pred[np.isfinite(pred)], np.linspace(0.02, 0.98, 97))
    out = np.zeros(len(pred), bool)
    cuts = []
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0:
            continue
        best = max(grid, key=lambda t: matthews_corrcoef(truth_bool[trn], pred[trn] >= t))
        out[te] = pred[te] >= best
        cuts.append(float(best))
    return out, cuts


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/gate.json")
    a = ap.parse_args()

    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tdi = pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv").set_index("Molecule_Name")
    inh = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv").set_index("Molecule_Name")
    T = tdi.reindex(rows.Molecule_Name)
    I = inh.reindex(rows.Molecule_Name)

    recs = []
    for seed in [int(s) for s in a.seeds.split(",")]:
        fold, ncl = butina_folds(list(rows.SMILES), seed=seed)
        print(f"\n=== сид {seed}: {ncl} кластеров ===", flush=True)
        for c in TDI_CYPS:
            lab = T[f"{c}_is_TDI"]
            dr = I[f"{c}_pIC50_direct_inhibition"].to_numpy(float)
            ta = T[f"{c}_pIC50_TDI_condition"].to_numpy(float)
            m = (lab.notna() & np.isfinite(dr) & np.isfinite(ta)).to_numpy()
            y = lab[m].astype(bool).to_numpy()
            Xi, fi = X[m], fold[m]
            dr_, ta_ = dr[m], ta[m]
            dl = ta_ - dr_
            print(f"  {c}: n {len(y)}, положительных {int(y.sum())} ({y.mean():.3f})", flush=True)

            ta_hat = oof_reg(Xi, ta_, fi, f"{c} плечо TDI")
            dl_hat = oof_reg(Xi, dl, fi, f"{c} Delta")

            gate_t, shift_t = ta_ > GATE, dl > L
            gate_h, shift_h = ta_hat > GATE, dl_hat > L
            gate_f, cuts = fitted_cut(ta_hat, gate_t, fi)

            assert (gate_t & shift_t == y).all(), "правило не воспроизвелось -- стенд сломан"

            arms = {
                "ворота ИСТИННЫЕ (оракул)":        gate_t,
                "ворота предск., порог 4.301":     gate_h,
                "ворота предск., порог подогнан":  gate_f,
                "сдвиг ИСТИННЫЙ (оракул)":         shift_t,
                "сдвиг предсказан":                shift_h,
                "оба предсказаны (4.301)":         gate_h & shift_h,
                "оба предсказаны (подогн.)":       gate_f & shift_h,
                "ворота ист. + сдвиг предск.":     gate_t & shift_h,
                "ворота предск. + сдвиг ист.":     gate_f & shift_t,
                "правило целиком (оракул)":        gate_t & shift_t,
            }
            print(f"    RMSE плеча TDI {np.sqrt(np.mean((ta_hat-ta_)**2)):.3f}, "
                  f"Delta {np.sqrt(np.mean((dl_hat-dl)**2)):.3f}; "
                  f"AUC предсказанного плеча по воротам {roc_auc_score(gate_t, ta_hat):.4f}", flush=True)
            print(f"    подогнанные пороги на предсказанном плече: "
                  f"{' '.join(f'{x:.2f}' for x in cuts)} (истинный 4.301)", flush=True)
            print(f"    ворота: истинных открыто {gate_t.mean():.3f}, "
                  f"предсказанных при 4.301 {gate_h.mean():.3f}, при подогнанном {gate_f.mean():.3f}",
                  flush=True)
            print(f"    восстановление самих ворот: MCC(предск., ист.) при 4.301 "
                  f"{matthews_corrcoef(gate_t, gate_h):.4f}, при подогнанном "
                  f"{matthews_corrcoef(gate_t, gate_f):.4f}", flush=True)
            for k, v in arms.items():
                mcc = matthews_corrcoef(y, v)
                print(f"      {k:32s} MCC {mcc:+.4f}  положительных {v.mean():.3f}", flush=True)
                recs.append({"seed": seed, "cyp": c, "arm": k, "mcc": float(mcc),
                             "rate": float(v.mean()), "n": int(len(y))})
        json.dump(recs, open(a.out, "w"))

    df = pd.DataFrame(recs)
    print("\n\nMCC по плечам, среднее по сидам")
    print(df.pivot_table(index="arm", columns="cyp", values="mcc").round(4).to_string())
    print(f"\nсохранено: {a.out}", flush=True)


if __name__ == "__main__":
    main()
