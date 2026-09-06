"""Three cheap things the conjunction makes obvious, none of which has been tried.

Item 238 measured that our CYP3A4 TDI classifier is, to within the floor, a potency threshold:
the predicted gate alone scores 0.3043 against the deployed classifier's 0.315, and the predicted
shift adds +0.027. It also left 53.7 per cent recovery of a gate whose oracle is +0.5667. This file
attacks that gap three ways, and each is free in the sense that it needs no new chemistry.

**1. The 1238 rows that are poison for the label and clean for the arm.** Item 231 closed the
molecules in the TDI table that carry no direct-inhibition arm: 1238 of them are labelled
`is_TDI = False`, 84.7 per cent clear the gate, and zero are positive against a base rate of 0.326,
so the label is a placeholder and training a classifier on it would be training on a fiction.

But the poison is in the LABEL, not in the arm. Their `pi_TDI` is a real measurement, and the gate
is a threshold on `pi_TDI`. They sit outside `data/rows.csv` -- outside the 4905 molecules the
project has features for -- which is why nothing has used them. All 1238 parse under RDKit and none
appears in the test set by name or by SMILES. Featurising them adds **53 per cent** to the training
set of the one regressor that carries CYP3A4's whole classification score.

They are not a free lunch: their median `pi_TDI` is 5.40 against 4.63 for the modelled set and
their gate is open 84.7 per cent against 60.2, so they are a shifted population. That is measured
here rather than assumed, and they enter TRAINING folds only -- never a test fold, since we do not
trust anything about them except the arm.

**2. The shift's cut has never been fitted, and item 238 proved that class of defect is real.**
Item 229 measured the Delta predictor shrunk to 0.42-0.46 of the true spread. Thresholding a shrunk
prediction at the true threshold is the wrong rule -- exactly what item 238 found for the gate, where
CYP2D6's out-of-fold optimum was 4.69-4.98 against a true 4.301. Item 238 fitted the gate's cut and
left the shift's at log10(2). One number, never measured.

**3. The rule is a conjunction and the classifier does not know it.** Instead of two hard
thresholds, form P(Delta > log10 2) and P(pi_TDI > 4.301) and multiply. Under the arms' error
correlation of 0.84 (item 229) that product is not the joint probability, so it is a score to be
thresholded rather than a calibrated probability -- but it uses the label's structure, which the
structural classifier throws away entirely.

Cuts are fitted OUT OF FOLD throughout; the joint arm fits both on a 2-D grid on training folds only.
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
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.metrics import matthews_corrcoef, roc_auc_score

import feats as F
from cypsplit import butina_folds

KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
L = np.log10(2.0)
GATE = 4.0 + L


def extra_rows(desc_names, mech_names, cyp="CYP3A4"):
    """Молекулы файла TDI, которых нет в rows.csv, с измеренным плечом TDI."""
    tdi = pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv")
    rows = pd.read_csv(D + "rows.csv")
    out = tdi[~tdi.Molecule_Name.isin(rows.Molecule_Name)].reset_index(drop=True)
    out = out[np.isfinite(out[f"{cyp}_pIC50_TDI_condition"].to_numpy(float))].reset_index(drop=True)
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    assert not out.Molecule_Name.isin(te.Molecule_Name).any(), "внешние строки пересеклись с тестом"
    assert not out.SMILES.isin(te.SMILES).any(), "внешние SMILES пересеклись с тестом"
    FP, dsc, M, ok = F.build(list(out.SMILES), desc_names, mech_names)
    if len(ok) != len(out):
        out = out.iloc[ok].reset_index(drop=True)
    X = np.hstack([FP, dsc.to_numpy(float), M.to_numpy(float)])
    return X, out[f"{cyp}_pIC50_TDI_condition"].to_numpy(float)


def fit_cut(pred, truth, fold, lo=None, hi=None):
    g = np.quantile(pred[np.isfinite(pred)], np.linspace(0.02, 0.98, 97))
    if lo is not None:
        g = g[(g >= lo) & (g <= hi)] if ((g >= lo) & (g <= hi)).any() else g
    out = np.zeros(len(pred), bool)
    cuts = []
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0:
            continue
        b = max(g, key=lambda t: matthews_corrcoef(truth[trn], pred[trn] >= t))
        out[te] = pred[te] >= b
        cuts.append(float(b))
    return out, cuts


def fit_cut2(pa, pb, y, fold):
    """Совместная подгонка ДВУХ порогов на обучающих фолдах."""
    ga = np.quantile(pa, np.linspace(0.05, 0.95, 25))
    gb = np.quantile(pb, np.linspace(0.05, 0.95, 25))
    out = np.zeros(len(y), bool)
    cuts = []
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0:
            continue
        best, bv = (ga[0], gb[0]), -2
        for ta in ga:
            for tb in gb:
                v = matthews_corrcoef(y[trn], (pa[trn] >= ta) & (pb[trn] >= tb))
                if v > bv:
                    bv, best = v, (ta, tb)
        out[te] = (pa[te] >= best[0]) & (pb[te] >= best[1])
        cuts.append([float(best[0]), float(best[1])])
    return out, cuts


def oof(X, y, fold, Xex=None, yex=None, clf=False, tag=""):
    p = np.full(len(y), np.nan)
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0:
            continue
        Xt, yt = X[trn], y[trn]
        if Xex is not None:
            Xt, yt = np.vstack([Xt, Xex]), np.concatenate([yt, yex])
        t0 = time.time()
        if clf:
            p[te] = (HistGradientBoostingClassifier(**KW).fit(Xt, yt.astype(int))
                     .predict_proba(X[te])[:, 1])
        else:
            p[te] = HistGradientBoostingRegressor(**KW).fit(Xt, yt).predict(X[te])
        print(f"      {tag} фолд {f}: {len(yt)} обуч, {time.time()-t0:.0f} с", flush=True)
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/gatedata.json")
    a = ap.parse_args()

    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    dn = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
    mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    tdi = pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv").set_index("Molecule_Name").reindex(rows.Molecule_Name)
    inh = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv").set_index("Molecule_Name").reindex(rows.Molecule_Name)

    print("строю признаки для внешних строк (переиндексация по именам)", flush=True)
    Xex, yex = extra_rows(dn, mn)
    if Xex.shape[1] != X.shape[1]:
        raise SystemExit(f"ширина не совпала: {Xex.shape[1]} против {X.shape[1]}")
    print(f"  внешних {Xex.shape}, ворота открыты у {(yex>GATE).mean():.3f}, "
          f"медиана pi_TDI {np.median(yex):.2f}", flush=True)

    recs = []
    for seed in [int(s) for s in a.seeds.split(",")]:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        print(f"\n=== сид {seed} ===", flush=True)
        for c in ["CYP3A4", "CYP2D6"]:
            lab = tdi[f"{c}_is_TDI"]
            ta = tdi[f"{c}_pIC50_TDI_condition"].to_numpy(float)
            dr = inh[f"{c}_pIC50_direct_inhibition"].to_numpy(float)
            m = (lab.notna() & np.isfinite(ta) & np.isfinite(dr)).to_numpy()
            y = lab[m].astype(bool).to_numpy()
            Xi, fi, ta_, dr_ = X[m], fold[m], ta[m], dr[m]
            dl = ta_ - dr_
            gate_t = ta_ > GATE
            use_ex = (c == "CYP3A4")
            print(f"  {c}: n {len(y)}, положительных {y.mean():.3f}"
                  f"{'  (+внешние)' if use_ex else ''}", flush=True)

            ta_hat = oof(Xi, ta_, fi, tag=f"{c} плечо")
            dl_hat = oof(Xi, dl, fi, tag=f"{c} Delta")
            gp = oof(Xi, gate_t, fi, clf=True, tag=f"{c} ворота-клф")
            ta_hat_x = oof(Xi, ta_, fi, Xex, yex, tag=f"{c} плечо+внеш") if use_ex else None
            gp_x = oof(Xi, gate_t, fi, Xex, yex > GATE, clf=True,
                       tag=f"{c} ворота-клф+внеш") if use_ex else None
            dp = oof(Xi, dl > L, fi, clf=True, tag=f"{c} сдвиг-клф")

            sh_f, sh_cuts = fit_cut(dl_hat, y, fi)
            arms = {
                "ворота: регр, срез 4.301":        ta_hat > GATE,
                "ворота: регр, срез подогнан":     fit_cut(ta_hat, gate_t, fi)[0],
                "ворота: классификатор":           fit_cut(gp, gate_t, fi)[0],
                "сдвиг: регр, срез 0.301":         dl_hat > L,
                "сдвиг: регр, срез подогнан":      sh_f,
                "сдвиг: классификатор":            fit_cut(dp, y, fi)[0],
                "оба, срезы книжные":              (ta_hat > GATE) & (dl_hat > L),
                "оба, срезы подогнаны раздельно":  fit_cut(ta_hat, gate_t, fi)[0] & sh_f,
                "оба, срезы подогнаны СОВМЕСТНО":  fit_cut2(ta_hat, dl_hat, y, fi)[0],
                "произведение вероятностей":       fit_cut(gp * dp, y, fi)[0],
            }
            if use_ex:
                arms["ворота: регр+внешние, 4.301"] = ta_hat_x > GATE
                arms["ворота: регр+внешние, подогн"] = fit_cut(ta_hat_x, gate_t, fi)[0]
                arms["ворота: клф+внешние"] = fit_cut(gp_x, gate_t, fi)[0]
                arms["оба: ворота+внешние, совместно"] = fit_cut2(ta_hat_x, dl_hat, y, fi)[0]
                print(f"    RMSE плеча: без внешних {np.sqrt(np.mean((ta_hat-ta_)**2)):.4f}, "
                      f"с внешними {np.sqrt(np.mean((ta_hat_x-ta_)**2)):.4f}; "
                      f"AUC по воротам {roc_auc_score(gate_t,ta_hat):.4f} -> "
                      f"{roc_auc_score(gate_t,ta_hat_x):.4f}", flush=True)
            print(f"    подогнанные срезы сдвига: {' '.join(f'{x:.3f}' for x in sh_cuts)} "
                  f"(книжный 0.301); усадка sd {np.std(dl_hat)/np.std(dl):.3f}", flush=True)
            for k, v in arms.items():
                mcc = matthews_corrcoef(y, v)
                print(f"      {k:34s} MCC {mcc:+.4f}  доля {v.mean():.3f}", flush=True)
                recs.append({"seed": seed, "cyp": c, "arm": k, "mcc": float(mcc),
                             "rate": float(v.mean())})
        json.dump(recs, open(a.out, "w"))

    df = pd.DataFrame(recs)
    print("\n\nMCC, среднее по сидам")
    print(df.pivot_table(index="arm", columns="cyp", values="mcc").round(4).to_string())
    print(f"\nсохранено: {a.out}", flush=True)


if __name__ == "__main__":
    main()
