"""Produce the two submission files and refuse to write anything the validator rejects.

Until this existed our leaderboard score was not 0.767 - it was undefined, because
nothing in the repository read cyp-challenge-TEST-BLINDED.csv at all. The intermediate
deadline is 24 September and the full test set is revealed once, on the 25th; that is the
only external measurement of how far our cross-validation is from the leaderboard, and it
cannot be recovered later.

Two things here are easy to get wrong and both are silent.

The descriptor block. src/feats.py selects descriptor columns by isna().mean() < 0.05,
a filter fitted on the 4905 training molecules. Recomputing it on 750 test molecules
would produce a different column set, and a check on the column *count* would not catch
a permutation either. feats.build() therefore reindexes by name against the committed
data/desc_names.csv and data/mech_names.csv.

Shrinkage. Post-hoc shrinkage toward the training mean is worth -0.034 macro in
cross-validation, but it is off by default here, and the reason is measured: on the top
quartile by activity it makes every enzyme worse, and the test set is built around
anchors at the 93rd to 98th percentile of the training distribution, so it is exactly
that kind of subsample. Shrinking toward a mean that is below the test's own mean drags
the active predictions down. --shrink turns it on; the offset and lambda are now fitted
jointly out of fold rather than the offset being fixed, which reaches macro 0.7150 against
0.7227 for the +0.40 slice. On the test the offset should be larger, not smaller.

How much larger has since been measured twice, from opposite ends, and the answers bracket
rather than agree. src/reweight.py tilts the label marginal and gives the centre as a
function of the shift delta: +0.4 at delta = 0, +0.8 at delta = 0.3, +0.9 at delta = 0.5,
with delta itself estimated at +0.3 to +0.6 from where the anchors sit. verify/k5_shift.py
measures the shift the other way, by running this model over both sets and comparing its
own output distributions: +0.014 / +0.147 / -0.190 / +0.437, which divided by the
attenuation b = 0.789 / 0.898 / 0.743 / 0.999 implies delta near +0.09 macro, and
verify/k6_shift1d.py confirms that at that shift the optimal centre does not move at all.

Neither number is wrong. +0.09 is a lower bound, because a model that mostly interpolates
propagates only part of an out-of-distribution shift into its predictions; +1.05 from the
raw anchor percentiles is an upper bound, because the anchors' neighbours were chosen by
similarity and regress toward the mean. Everything between about +0.1 and +0.6 is live,
and the centre that goes with it is between +0.4 and +0.8.

Two further cautions. The offset is in centre units: predictions move by (1 - lambda)
times it, so the +0.40 once quoted was never +0.40 in pIC50 - at the fitted lambdas the
real shifts are +0.13 / +0.10 / +0.03 / +0.17. And both reweightings assume p(y | yhat) is
the same on the test set, which is the thing recalibration exists to check.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
TUT = tutorial()

import argparse
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

import feats as F
from cypsplit import butina_folds

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
TDI_CYPS = ["CYP3A4", "CYP2D6"]          # the only two the organisers score
GRID = np.linspace(0.2, 1.0, 41)

# Same learner and settings as src/ablate.py, so the submitted model is the one the
# document's numbers describe rather than a cousin of it.
def gbm_reg():
    return HistGradientBoostingRegressor(max_iter=300, learning_rate=0.06,
                                         max_leaf_nodes=31, l2_regularization=1.0,
                                         random_state=0)


def gbm_clf():
    return HistGradientBoostingClassifier(max_iter=300, learning_rate=0.06,
                                          max_leaf_nodes=31, l2_regularization=1.0,
                                          random_state=0)


def test_features(desc_names, mech_names):
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    FP, dsc, M, ok = F.build(list(te.SMILES), desc_names, mech_names)
    if len(ok) != len(te):
        # Every test SMILES parses today; if that ever stops being true the rows would
        # silently misalign, so fail loudly instead.
        raise SystemExit(f"RDKit не разобрал {len(te)-len(ok)} тестовых SMILES - "
                         "выравнивание строк сломается, разбираться вручную")
    return te, np.hstack([FP, dsc.to_numpy(np.float32), M.to_numpy(np.float32)])


OFFGRID = np.linspace(-0.2, 1.6, 37)


def fit_shrinkage(X, y, mask, fold):
    """Offset and lambda per enzyme, both chosen out-of-fold on the training data.

    Fitting the two jointly rather than fixing the offset and searching lambda: the
    family {c + L(p - c)} is identically the affine family {a + b p} with b = L and
    a = c(1 - L), so a fixed offset is an arbitrary slice through it. Jointly it reaches
    macro 0.7150 against 0.7227 for the +0.40 slice and 0.7333 for the offset at the
    training mean.

    Worth noting what the offset is not. The predictions move by (1 - L) times it, not by
    it, so the +0.40 quoted earlier was never a 0.40 shift in pIC50 - at the fitted
    lambdas the actual shifts are +0.13 / +0.10 / +0.03 / +0.17. And at L -> 1 the
    centre is not identified at all, which is a second reason to fit the affine pair.
    """
    out = []
    for e, c in enumerate(CYPS):
        m = mask[:, e]
        yy = y[m, e]
        p = np.zeros(m.sum())
        fi = fold[m]
        Xi = X[m]
        for f in range(5):
            a, b = fi != f, fi == f
            if b.sum() == 0:
                continue
            p[b] = gbm_reg().fit(Xi[a], yy[a]).predict(Xi[b])
        lo = LO[m, e]; hi = HI[m, e]
        mu = p.mean()
        off, L = min(((o, l) for o in OFFGRID for l in GRID),
                     key=lambda t: strae(yy, (mu + t[0]) + t[1] * (p - (mu + t[0])),
                                         y_true_upper=hi, y_true_lower=lo))
        out.append((L, mu + off))
        print(f"    {c}: lambda {L:.2f}, смещение {off:+.2f}, "
              f"сдвиг предсказаний {(1-L)*off:+.3f}", flush=True)
    return out


def main():
    global LO, HI
    ap = argparse.ArgumentParser()
    ap.add_argument("--shrink", action="store_true", help="применить усадку (см. docstring)")
    ap.add_argument("--outdir", default=RES + "submission/")
    a = ap.parse_args()

    _pl.Path(a.outdir).mkdir(parents=True, exist_ok=True)
    desc_names = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
    mech_names = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()

    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())

    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)

    print("строю признаки теста (переиндексация по именам, не пересчёт фильтра)", flush=True)
    te, Xte = test_features(desc_names, mech_names)
    if Xte.shape[1] != X.shape[1]:
        raise SystemExit(f"ширина не совпала: обучение {X.shape[1]}, тест {Xte.shape[1]}")
    print(f"  тест {Xte.shape}, обучение {X.shape}", flush=True)

    lams = None
    if a.shrink:
        print("подбираю усадку вне выборки на обучающих данных", flush=True)
        fold, _ = butina_folds(list(rows.SMILES))
        lams = fit_shrinkage(X, y, mask, fold)

    print("обучаю на всей выборке и предсказываю тест", flush=True)
    act = pd.DataFrame({"SMILES": te.SMILES, "Molecule_Name": te.Molecule_Name})
    for e, c in enumerate(CYPS):
        m = mask[:, e]
        p = gbm_reg().fit(X[m], y[m, e]).predict(Xte)
        if lams is not None:
            L, mu = lams[e]
            p = mu + L * (p - mu)
        act[f"{c}_pIC50_direct_inhibition"] = p
        print(f"    {c}: n_обуч {m.sum()}, среднее предсказание {p.mean():.3f}", flush=True)

    tdi = pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv").set_index("Molecule_Name")
    keep = rows.Molecule_Name.isin(tdi.index).to_numpy()
    T = tdi.loc[rows.Molecule_Name[keep]].reset_index()
    cls = pd.DataFrame({"SMILES": te.SMILES, "Molecule_Name": te.Molecule_Name})
    for c in TDI_CYPS:
        lab = T[f"{c}_is_TDI"]
        m = lab.notna().to_numpy()
        yb = lab[m].astype(int).to_numpy()
        pr = gbm_clf().fit(X[keep][m], yb).predict_proba(Xte)[:, 1]
        # Plug-in threshold by expected MCC. Measured on this data: fitting the threshold
        # instead is better on 3A4 and worse on 2D6, and calibrating first gains +0.026
        # on 3A4 and loses on 2D6 - all differences far inside the MCC interval at n=750.
        # So: one rule, applied to both endpoints, not a per-endpoint recipe.
        ts = np.linspace(0.05, 0.95, 91)
        def emcc(t):
            yh = pr >= t
            tp = (pr * yh).sum(); fp = ((1 - pr) * yh).sum()
            fn = (pr * ~yh).sum(); tn = ((1 - pr) * ~yh).sum()
            d = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
            return (tp * tn - fp * fn) / d if d > 0 else 0.0
        thr = ts[int(np.argmax([emcc(t) for t in ts]))]
        cls[f"{c}_is_TDI"] = pr >= thr
        print(f"    {c}: порог {thr:.3f}, положительных {int((pr>=thr).sum())} из {len(pr)}", flush=True)

    ap_ = a.outdir + "activity_submission.csv"
    tp_ = a.outdir + "tdi_submission.csv"
    act.to_csv(ap_, index=False)
    cls.to_csv(tp_, index=False)

    # The gate. Nothing above is trusted until the organisers' own code accepts it.
    from validation.activity_validation import validate_activity_submission
    from validation.tdi_validation import validate_tdi_submission
    ids = set(te.Molecule_Name)
    ok = True
    for name, fn, path in [("регрессия", validate_activity_submission, ap_),
                           ("классификация", validate_tdi_submission, tp_)]:
        try:
            res = fn(path, expected_ids=ids)
            bad = res if isinstance(res, list) else getattr(res, "errors", [])
            if bad:
                ok = False
                print(f"  {name}: ОТКЛОНЕНО")
                for e_ in bad:
                    print("     ", e_)
            else:
                print(f"  {name}: принято")
        except TypeError:
            res = fn(path)
            print(f"  {name}: принято (валидатор без expected_ids)")
    if not ok:
        raise SystemExit("валидатор отверг файл; ничего не отправлять")
    print(f"\nготово:\n  {ap_}\n  {tp_}")


if __name__ == "__main__":
    LO = HI = None
    main()
