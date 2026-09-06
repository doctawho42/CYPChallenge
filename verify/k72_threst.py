"""The threshold is an estimator, and argmax is a bad one. Item 235 left +0.047 on the table.

Item 235 deployed Platt calibration and took +0.0133 of macro MCC. The threshold ORACLE -- a cut
chosen knowing the scored fold's own labels -- reached +0.0472 with the sign in 8 cells of 8. So
three quarters of the available threshold gain is still unclaimed, and the arm that was supposed to
claim it, a cut fitted by empirical MCC on training labels, went 4 of 8 and was not adopted.

**Why it failed is not a shortage of data.** In item 235 that cut was fitted on the WHOLE training
part of each outer fold -- about 1877 rows, four fifths of the set, not a fifth. It still came out
at 0.05, 0.70, 0.05, 0.41, 0.14 across five folds of one seed on CYP2D6. At AUC 0.588 the empirical
MCC curve is flat near its top, so its `argmax` is a high-variance statistic no matter how many rows
go into it. **The fix is a better estimator of the same quantity, not more rows.**

Seven estimators of one scalar, all fitted on training labels only, all applied out of fold:

    plug-in         максимум ОЖИДАЕМОГО MCC при собственных вероятностях (меток не требует)
    argmax          максимум эмпирического MCC -- то, что провалилось в пункте 235
    сглаженный      argmax эмпирической кривой после скользящего среднего по сетке
    центроид        середина области, лежащей в одной бутстрап-ошибке от максимума
    бэггинг         среднее argmax по 200 бутстрап-перевыборкам
    квантиль        срез, дающий на обучении ту же долю положительных, что и базовая
    оракул          argmax на метках самого оцениваемого фолда -- потолок, не для подачи

Each is run on raw probabilities and on Platt-calibrated ones, since item 235 put Platt into the
submission and the question is what to put on top of it.

**One honesty note about the design.** The out-of-fold probabilities here come from a single level
of cross-validation, so a row's probability was produced by a model that saw the fold being scored.
That leak inflates every arm equally, and this file compares ESTIMATORS OF A THRESHOLD against each
other rather than reporting a deployable MCC -- item 235's fully nested numbers remain the ones to
quote. What transfers from here is the ranking and the size of the gaps, not the absolute level.
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
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import matthews_corrcoef

from cypsplit import butina_folds
from submit import gbm_clf, plugin_threshold, TDI_CYPS

GRID = np.linspace(0.02, 0.95, 187)


def mcc_curve(p, y, grid=GRID):
    return np.array([matthews_corrcoef(y, p >= t) for t in grid])


def est_argmax(p, y):
    return float(GRID[int(np.argmax(mcc_curve(p, y)))])


def est_smooth(p, y, w=15):
    c = mcc_curve(p, y)
    k = np.ones(w) / w
    s = np.convolve(np.pad(c, w // 2, mode="edge"), k, mode="valid")[:len(c)]
    return float(GRID[int(np.argmax(s))])


def est_centroid(p, y, n=200, seed=0):
    c = mcc_curve(p, y)
    rng = np.random.default_rng(seed)
    peaks = []
    for _ in range(n):
        k = rng.integers(0, len(y), len(y))
        peaks.append(mcc_curve(p[k], y[k]).max())
    se = float(np.std(peaks))
    ok = c >= c.max() - se
    return float(GRID[ok].mean())


def est_bagged(p, y, n=200, seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(n):
        k = rng.integers(0, len(y), len(y))
        out.append(GRID[int(np.argmax(mcc_curve(p[k], y[k])))])
    return float(np.mean(out))


def est_quantile(p, y):
    return float(np.quantile(p, 1.0 - y.mean()))


ESTS = {"plug-in": lambda p, y: plugin_threshold(p, GRID),
        "argmax": est_argmax, "сглаженный": est_smooth,
        "центроид": est_centroid, "бэггинг": est_bagged, "квантиль": est_quantile}


def logit(q):
    q = np.clip(q, 1e-6, 1 - 1e-6)
    return np.log(q / (1 - q)).reshape(-1, 1)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/threst.json")
    a = ap.parse_args()

    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tdi = pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv").set_index("Molecule_Name")
    keep = rows.Molecule_Name.isin(tdi.index).to_numpy()
    T = tdi.loc[rows.Molecule_Name[keep]].reset_index()
    Xt, smi = X[keep], list(rows.SMILES[keep])

    recs = []
    for seed in [int(s) for s in a.seeds.split(",")]:
        fold, _ = butina_folds(smi, seed=seed)
        print(f"\n=== сид {seed} ===", flush=True)
        for c in TDI_CYPS:
            lab = T[f"{c}_is_TDI"]
            m = lab.notna().to_numpy()
            y = lab[m].astype(int).to_numpy()
            Xi, fi = Xt[m], fold[m]
            p = np.zeros(len(y))
            t0 = time.time()
            for f in range(5):
                trn, te = fi != f, fi == f
                if te.sum() == 0:
                    continue
                p[te] = gbm_clf().fit(Xi[trn], y[trn]).predict_proba(Xi[te])[:, 1]
            print(f"  {c}: n {len(y)}, положительных {y.mean():.3f}, "
                  f"{time.time()-t0:.0f} с", flush=True)

            for cal in ("сырые", "Платт"):
                dec = {k: np.zeros(len(y), bool) for k in list(ESTS) + ["оракул"]}
                cuts = {k: [] for k in ESTS}
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    if cal == "сырые":
                        ptr, pte = p[trn], p[te]
                    else:
                        lr = LogisticRegression(C=1e6).fit(logit(p[trn]), y[trn])
                        ptr = lr.predict_proba(logit(p[trn]))[:, 1]
                        pte = lr.predict_proba(logit(p[te]))[:, 1]
                    for k, fn in ESTS.items():
                        t = fn(ptr, y[trn])
                        dec[k][te] = pte >= t
                        cuts[k].append(t)
                    dec["оракул"][te] = pte >= est_argmax(pte, y[te])
                base = matthews_corrcoef(y, dec["plug-in"])
                for k in dec:
                    mcc = matthews_corrcoef(y, dec[k])
                    sp = (f"  срезы {' '.join(f'{x:.2f}' for x in cuts[k])}"
                          if k in cuts else "")
                    print(f"    {cal:6s} {k:12s} MCC {mcc:+.4f}  против plug-in "
                          f"{mcc-base:+.4f}{sp}", flush=True)
                    recs.append({"seed": seed, "cyp": c, "cal": cal, "est": k,
                                 "mcc": float(mcc), "vs_plug": float(mcc - base)})
        json.dump(recs, open(a.out, "w"))

    df = pd.DataFrame(recs)
    print("\n\nMCC, среднее по сидам")
    print(df.pivot_table(index=["cal", "est"], columns="cyp", values="mcc").round(4).to_string())
    print("\nпротив plug-in, среднее и знак по 8 клеткам (сид x фермент)")
    for (cal, est), g in df.groupby(["cal", "est"]):
        if est == "plug-in":
            continue
        print(f"  {cal:6s} {est:12s} {g.vs_plug.mean():+.4f}  знак "
              f"{int((g.vs_plug > 0).sum())}/{len(g)}")
    print(f"\nсохранено: {a.out}", flush=True)


if __name__ == "__main__":
    main()
