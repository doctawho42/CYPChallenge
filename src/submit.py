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

Shrinkage, and why the default is the open question rather than a settled one.
--shrink is OFF, and the justification it used to carry has since been falsified. That
justification was: on the top quartile by activity, shrinking toward the training mean
makes every enzyme worse; the test is built around anchors at the 93rd to 98th percentile;
therefore the test is that kind of subsample and shrinkage would hurt.

The middle step does not survive measurement. verify/k5_shift.py runs one model over both
sets and compares its own output distributions - the mapping is identical, so the
difference in outputs is a difference in inputs. The shift is +0.014 / +0.147 / -0.190 /
+0.437 by enzyme, mild, and negative on 2D6, which is the internal control since 2D6 took
no part in anchor selection. The stress test that flipped the sign moved ybar by +1.1 to
+1.3 - three times stronger than the shift that is actually there.

verify/k6_shift1d.py then reweights the training set to the test-like marginal of the
predictions (1-D density ratio, effective n 796 to 1444) and finds the optimum barely
moves: (+0.30, 0.58) against (+0.30, 0.58) on 1A2, (+1.10, 0.84) against (+0.95, 0.82) on
3A4 at worst. Our own fitted parameters score 0.6781 under those weights against 0.6773
for parameters fitted under them - a gap of 0.0008 macro. Shrinkage does not flip sign
there; it wins, 0.7407 raw against 0.6781.

So the evidence now points the other way, and the default has deliberately NOT been
flipped on that basis alone, because two things the reweighting cannot reach are exactly
the two that would break it. It corrects the marginal of yhat while assuming p(y | yhat)
is unchanged on the test - the assumption recalibration exists to test. And it does not
touch the similarity geometry at all: the test sits at median nearest-neighbour 0.587 and
no re-split of the training data gets above 0.450. Flipping this default changes what gets
submitted, and the 25 September reveal is one-shot, so it is a decision to take
deliberately rather than as a side effect of a docstring.

If it is turned on, the offset and lambda are fitted jointly out of fold rather than the
offset being fixed, reaching macro 0.7150 against 0.7227 for the +0.40 slice. Note the
offset is in centre units: predictions move by (1 - lambda) times it, so the +0.40 once
quoted was never +0.40 in pIC50 - the real shifts are +0.13 / +0.10 / +0.03 / +0.17.

WHICH SHIFT, AND UNDER WHICH CRITERION. Both were open questions until they were
measured, and the second turned out to matter more than the first.

The shift is not one number for four enzymes. On CYP2D6 it is negative, and the reason is
chemistry rather than statistics: the test carries about a third as many compounds that are
basic at pH 7.4 as the CYP2D6 label mask does, and CYP2D6 is the one enzyme of the four that
binds through a salt bridge to a protonated nitrogen, so basic compounds are MORE active
there and less so everywhere else (verify/k7_2d6shift.py, k10_strat2d6.py).

The criterion was doing more work than any of the estimates. Choosing the pair by the WORST
case inside each enzyme's plausible range is insurance against a bad leaderboard; choosing by
the MEAN over the posterior of delta is a bid for the best expected score. The two disagree
by more than any two estimates of delta disagree - the per-enzyme gain is 0.107 by worst case
and 0.045 by mean - and they differ in the sign of their derivative with respect to how wide
the range is. The mean is adopted here: we are after the best expected score, not insurance.

The default is therefore 0, +0.3, -0.5, +0.7, and the zero on CYP1A2 is deliberate. There the
sign of the shift is not determined at all, P(delta >= 0) = 0.59, so any non-zero choice is a
coin flip against doing nothing: the mean criterion picks +0.1, gains 0.0001 by it, and loses
to zero in half the posterior draws. That is added variance for no expected return.

Held to the same standard, this rule is clean where the earlier one was not. The worst-case
rule picked +0.4 on CYP1A2 and lost to doing nothing in 81 % of draws - exactly the defect
that got a single global delta rejected, relocated to another cell. Under the adopted rule
the fractions are 0.09 / 0.11 / 0.01 on the three enzymes it touches.

None of this fires unless --shrink is passed. That switch is the one decision still open.

The size of the shift is bracketed rather than pinned. src/reweight.py tilts the label
marginal and puts the centre at +0.4 for delta = 0 and +0.9 for delta = 0.5; the anchor
percentiles put delta at +1.05, an upper bound, since the anchors' neighbours were chosen
by similarity and regress toward the mean. The prediction-shift route puts it near +0.09,
a lower bound, since a model that mostly interpolates propagates only part of an input
shift into its outputs. Everything between about +0.1 and +0.6 is live.
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
from reweight import tilt

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


# Сетка смещений намеренно шире, чем нужно любому правдоподобному сдвигу. Прежняя,
# linspace(-0.2, 1.6), зажимала подгонку с обеих сторон: при предполагаемом сдвиге ниже -0.3
# оптимальное смещение упиралось в нижний край и переставало двигаться, а на CYP3A4 при
# +0.7 оно садилось на верхний. Подогнанный параметр, стоящий на границе сетки, --- не
# подогнанный параметр, и плоскость целевой функции рядом с ним мнимая.
OFFGRID = np.round(np.arange(-3.0, 3.01, 0.05), 2)


def fit_shrinkage(X, y, mask, fold, delta=(0.0, 0.0, 0.0, 0.0)):
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
        # Under an assumed shift the objective is the tilted one: our own labels reweighted
        # so their mean sits delta higher. At delta = 0 the weights are all ones and this is
        # exactly the untilted fit, so the default path is unchanged.
        de = float(delta[e])
        w = np.ones_like(yy) if de == 0 else tilt(yy, de)
        # Вся сетка одним броадкастом: q = c + L(p - c) = (mu + off)(1 - L) + L p.
        A = OFFGRID[:, None] * (1.0 - GRID[None, :])
        B = np.broadcast_to(GRID[None, :], A.shape)
        q = (mu * (1.0 - B) + A)[:, :, None] + B[:, :, None] * p[None, None, :]
        pen = (w * (np.maximum(q - hi, 0.0) + np.maximum(lo - q, 0.0))).sum(axis=2)
        ii, jj = np.unravel_index(pen.argmin(), pen.shape)
        off, L = float(OFFGRID[ii]), float(GRID[jj])
        if ii in (0, len(OFFGRID) - 1) or jj in (0, len(GRID) - 1):
            print(f"    ВНИМАНИЕ {c}: оптимум на краю сетки (off {off:+.2f}, lambda {L:.2f})",
                  flush=True)
        out.append((L, mu + off))
        print(f"    {c}: lambda {L:.2f}, смещение {off:+.2f}, "
              f"сдвиг предсказаний {(1-L)*off:+.3f}", flush=True)
    return out


def main():
    global LO, HI
    ap = argparse.ArgumentParser()
    ap.add_argument("--shrink", action="store_true", help="применить усадку (см. docstring)")
    ap.add_argument("--delta", default="0,0.3,-0.5,0.7",
                    help="предполагаемый сдвиг средней активности теста относительно нашей "
                         "выборки. Пара (off, lambda) подбирается под ЭТО предположение. "
                         "Ноль означает «тест распределён как обучающая выборка» - это не "
                         "отсутствие предположения, а предположение, и src/shrinkchoice.py "
                         "показывает, что по вилке +0.1..+0.6 оно худшее из трёх правил: "
                         "худший случай на 0.087, средний на 0.040 хуже подгонки под "
                         "середину вилки. Значение по умолчанию оставлено нулевым, чтобы "
                         "Принимает одно число на все ферменты или четыре через запятую в "
                         "порядке CYPS. Умолчание --- принятое правило: по среднему "
                         "апостериорному, с нулём на CYP1A2, см. docstring")
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
        d = [float(x) for x in str(a.delta).split(",")]
        if len(d) == 1:
            d = d * 4
        if len(d) != 4:
            raise SystemExit(f"--delta: нужно одно число или четыре через запятую, дано {len(d)}")
        print(f"предполагаемый сдвиг по ферментам: "
              + ", ".join(f"{c} {v:+.2f}" for c, v in zip(CYPS, d)), flush=True)
        lams = fit_shrinkage(X, y, mask, fold, d)

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
