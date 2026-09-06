"""The TDI classifier has never seen our own pIC50 predictions, and the label is a function of them.

Classification is **a third of the leaderboard** -- two endpoints of six, by the organisers' own
`CLASSIFICATION_ENDPOINTS` -- and almost nothing in this file is about it. The regression track
carries a five-member ensemble, the dead zone, pooling, per-enzyme member selection and the affine
pair. The TDI track is one `HistGradientBoostingClassifier` on the feature matrix plus a plug-in
threshold. That asymmetry is the opportunity, and this file takes the most specific part of it.

**The observation this rests on.** Item 197 established that `is_TDI` is a *deterministic function*
of the two potency arms: applying the organisers' rule to the measured labels reproduces the flag on
**100.0 per cent** of compounds with zero errors either way on both enzymes. The label is not merely
correlated with what we model -- it is computed from it.

**And the classifier does not see it.** `src/submit.py` fits on the structural matrix alone. Nothing
in this file has ever given the classifier the regression track's own output.

**Why this is not item 197 again, which is the first objection.** That item applied the *rule* to
predicted quantities, integrating it over their joint distribution, and got calibration three times
better with MCC and AUC worse. Imposing the rule's functional form pushes prediction error through a
hard threshold, where it is amplified. Handing the same quantities over as **features** is a
different object: the classifier learns how much to trust them, including where their error is
large, instead of being told to believe them exactly.

Arms:

    структура                 контроль --- то, что подаётся сейчас
    + пи_прям                 плюс предсказание прямого плеча, вне фолда
    + пи_прям + пи_tdi + дельта   плюс оба плеча и их разность, вне фолда
    + правило                 плюс детерминированный выход правила на предсказаниях
    перемешан                 контроль: те же признаки, связь со строкой разорвана

**Nesting, because the trap here is obvious and fatal.** The regression predictions that become
features must be out of fold with respect to the classifier's own fold. They are computed once per
fold from the other four, so no compound's feature was built by a model that saw it.

Metric is MCC, which is what the leaderboard scores, with AUC alongside for diagnosis --- MCC
depends on a threshold and AUC does not, so a gain in one without the other says where it came from.
The threshold is the same plug-in expected-MCC rule `src/submit.py` uses, so the comparison is
against the submission rather than against an idealised classifier.

Writes results/preds/oof_tdif.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.metrics import matthews_corrcoef, roc_auc_score

from cypsplit import butina_folds

TDI_CYPS = ["CYP2D6", "CYP3A4"]
KWR = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
           l2_regularization=1.0, random_state=0)
KWC = dict(KWR)


def rule(direct, tdiarm):
    """Правило организаторов, воспроизводящее is_TDI. Проверяется на метках ниже."""
    return np.where(direct > 4.0, (tdiarm - direct) > np.log10(2.0), tdiarm > 4.301)


def plugin_threshold(pr):
    """Тот же подстановочный порог по ожидаемому MCC, что стоит в src/submit.py."""
    ts = np.linspace(0.05, 0.95, 91)
    best, bt = -2.0, 0.5
    for t in ts:
        yh = pr >= t
        tp = (pr * yh).sum(); fp = ((1 - pr) * yh).sum()
        fn = (pr * ~yh).sum(); tn = ((1 - pr) * ~yh).sum()
        d = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        v = (tp * tn - fp * fn) / d if d > 0 else 0.0
        if v > best:
            best, bt = v, t
    return bt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--arms", default="структура|+пи_прям|+оба плеча|+правило|перемешан")
    ap.add_argument("--out", default=RES + "preds/oof_tdif.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    arms = a.arms.split("|")

    rows = pd.read_csv(D + "rows.csv")
    inh = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
             .set_index("Molecule_Name").reindex(rows.Molecule_Name).reset_index())
    tdi = (pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv")
             .set_index("Molecule_Name").reindex(rows.Molecule_Name).reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)

    print("ПРЕДУСЛОВИЕ: воспроизводит ли правило метку на ИЗМЕРЕННЫХ плечах")
    for c in TDI_CYPS:
        d = inh[f"{c}_pIC50_direct_inhibition"].to_numpy(float)
        t = tdi[f"{c}_pIC50_TDI_condition"].to_numpy(float)
        f = tdi[f"{c}_is_TDI"]
        m = f.notna().to_numpy() & np.isfinite(d) & np.isfinite(t)
        agree = (rule(d[m], t[m]) == f[m].astype(bool).to_numpy()).mean()
        print(f"  {c}: {m.sum():5d} строк, совпадение {agree:.4f}, "
              f"положительных {int(f[m].sum())}")
    print()

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        rng = np.random.default_rng(5000 + seed)
        # Предсказания плеч вне фолда: считаются ОДИН раз и переиспользуются всеми руками,
        # чтобы руки различались только составом признаков.
        PRED = {}
        for c in TDI_CYPS:
            for nm, col, src in (("прям", f"{c}_pIC50_direct_inhibition", inh),
                                 ("tdi", f"{c}_pIC50_TDI_condition", tdi)):
                t0 = time.time()
                yv = src[col].to_numpy(float)
                obs = np.isfinite(yv)
                p = np.full(len(rows), np.nan)
                for f in range(5):
                    trn = obs & (fold != f)
                    te = fold == f
                    if te.sum() == 0 or trn.sum() < 50:
                        continue
                    p[te] = HistGradientBoostingRegressor(**KWR).fit(
                        X[trn], yv[trn]).predict(X[te])
                PRED[(c, nm)] = p
                print(f"    {c} {nm:5s} вне фолда посчитано ({time.time()-t0:.0f} с)", flush=True)

        for arm in arms:
            t0 = time.time()
            r = {"seed": seed, "рука": arm}
            for c in TDI_CYPS:
                lab = tdi[f"{c}_is_TDI"]
                m = lab.notna().to_numpy()
                yb = lab[m].astype(int).to_numpy()
                pd_, pt_ = PRED[(c, "прям")][m], PRED[(c, "tdi")][m]
                extra = []
                if arm in ("+пи_прям",):
                    extra = [pd_]
                elif arm in ("+оба плеча", "перемешан"):
                    extra = [pd_, pt_, pt_ - pd_]
                elif arm == "+правило":
                    extra = [pd_, pt_, pt_ - pd_, rule(pd_, pt_).astype(float)]
                if extra and arm == "перемешан":
                    idx = rng.permutation(len(yb))
                    extra = [v[idx] for v in extra]
                Xi = X[m] if not extra else np.hstack(
                    [X[m], np.column_stack(extra).astype(np.float32)])
                fi = fold[m]
                pr = np.zeros(len(yb))
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    pr[te] = HistGradientBoostingClassifier(**KWC).fit(
                        Xi[trn], yb[trn]).predict_proba(Xi[te])[:, 1]
                thr = plugin_threshold(pr)
                r[f"{c} MCC"] = round(float(matthews_corrcoef(yb, pr >= thr)), 4)
                r[f"{c} AUC"] = round(float(roc_auc_score(yb, pr)), 4)
                out[f"{seed}|{arm}|{c}"] = pr.tolist()
            r["MACRO MCC"] = round(float(np.mean([r[f"{c} MCC"] for c in TDI_CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {arm:14s} MCC макро {r['MACRO MCC']:.4f}  "
                  + "  ".join(f"{c[3:]} MCC {r[f'{c} MCC']:.3f} AUC {r[f'{c} AUC']:.3f}"
                              for c in TDI_CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("рука", sort=False)[["MACRO MCC"]
          + [f"{c} {k}" for c in TDI_CYPS for k in ("MCC", "AUC")]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Классификация --- треть лидерборда, и подаётся на ней один классификатор с порогом,
тогда как на регрессии стоит пять членов, мёртвая зона, пулирование и аффинная пара.

MCC и AUC вместе, потому что они отвечают на разные вопросы: AUC не зависит от порога, MCC
зависит. Выигрыш в MCC без выигрыша в AUC означает, что признаки помогли порогу, а не
упорядочению, и это слабее.

«перемешан» решает: те же добавленные столбцы, разорвана связь со строкой. Если рука с обоими
плечами идёт наравне с перемешанной, работает не предсказание потенции, а лишняя ёмкость.""")


if __name__ == "__main__":
    main()
