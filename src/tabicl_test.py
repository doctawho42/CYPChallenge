"""Fit the in-context learner on ALL training rows of each enzyme and predict the 750. Runs on the
GPU box; scores nothing.

The out-of-fold path is `src/tabicl_run.py`, and the two must agree in everything except what
changes: there the scaler and the PCA are fitted on each fold's four fitting folds, here on every
labelled row of the enzyme, which is what `submit.py:_probe_test_pred` does for the shipped probe
and what every other member's test path does. The estimator, the cap, the alpha grid and the random
state are identical, so a difference between the two cannot come from the recipe.

A RidgeCV on the identical post-PCA matrix is fitted in the same loop, exactly as out of fold, for
two reasons: it keeps item 372's preprocessing control alive on the test path, and it gives
`submit.py` a cheap arm to fall back to if the learner's test predictions ever look wrong.

THE CONTROL THAT MUST SUCCEED, and it is a question rather than an action, so it is not silenced.
A model fitted on 1285--2335 rows and asked for 750 predictions can fail two ways that an exit code
does not catch: it can collapse toward a constant, or it can leave the label range. So the per-enzyme
mean and standard deviation of the 750 predictions are printed beside the training labels' own, and
the run REFUSES if the spread falls below a fifth of the labels' or the mean moves by more than one
label standard deviation. Item 366's spread went to the board unchecked and cost 46 places.

    python src/tabicl_test.py --inp data/tabicl_test_in.npz --out out_test

Writes <out>/tabicl_test.json: {"preds": {rep: [[750 floats] x 4]}, "ridge": same, "pca": int}.

THE RECIPE, written out because a log whose script cannot be re-run is item 345's k123 defect and a
dependency that lives only in one shell's history is the same defect one step out. Measured on the
GPU box, 2026-10-23:

    uv pip install --python <env>/bin/python tabicl          # 2.2.0, no account, no licence

    tabicl 2.2.0   torch 2.13.0+cu130   scikit-learn 1.9.1   numpy 2.5.4   scipy 1.18.1   py3.12

The checkpoint downloads itself on first use from Hugging Face `jingang/TabICL`, snapshot
4dcd344ece2c00be9e831fdd35bed57b5ad83e19, file `tabicl-regressor-v2-20260212.ckpt`, **114324594
bytes**. No token and no accepted licence are involved, which is the whole reason this family is
TabICL and not TabPFN: TabPFN 9.2 gates its weights behind a registered account and an accepted
licence, and accepting one is the user's act, not this repository's.

**Note the scikit-learn: 1.9.1, not the pinned 1.3.2.** That is why this file writes a file and
scores nothing. `CLAUDE.md` records that 1.9.0 moves every prediction in this repository by ~0.135
pIC50 and doubles the measured effect of the mechanistic block, so every number that reaches the
journal is computed under the pin, by `verify/k147_tabicl2.py` and the scripts around it.
"""

import argparse
import json
import pathlib
import time

import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inp", default="/tmp/d42/data/tabicl_test_in.npz")
    ap.add_argument("--out", default="/tmp/d42/out_test")
    ap.add_argument("--reps", default="E,S")
    ap.add_argument("--pca", type=int, default=100)
    ap.add_argument("--device", default="cuda")
    a = ap.parse_args()
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import RidgeCV
    from tabicl import TabICLRegressor
    import torch

    z = np.load(a.inp, allow_pickle=True)
    Y = z["Y"].astype(np.float64)
    M = np.isfinite(Y)
    CY = [str(c) for c in z["cyps"]]
    n_te = z["E_te"].shape[0]
    out, ridge_out, stats = {}, {}, {}
    for rep in [x for x in a.reps.split(",") if x.strip()]:
        Rtr = z["%s_tr" % rep].astype(np.float64)
        Rte = z["%s_te" % rep].astype(np.float64)
        if Rtr.shape[1] != Rte.shape[1]:
            raise SystemExit("%s: ширина %d обучения против %d теста"
                             % (rep, Rtr.shape[1], Rte.shape[1]))
        P, Q, st = [], [], []
        t0 = time.time()
        for e in range(len(CY)):
            m = M[:, e]
            A, ym = Rtr[m], Y[m, e]
            sc = StandardScaler().fit(A)
            A, B = sc.transform(A), sc.transform(Rte)
            k = min(a.pca, A.shape[1], A.shape[0] - 1)
            pc = PCA(n_components=k, random_state=0).fit(A)
            A, B = pc.transform(A), pc.transform(B)
            mdl = TabICLRegressor(device=a.device, random_state=0)
            mdl.fit(A, ym)
            p = np.asarray(mdl.predict(B), dtype=np.float64)
            q = RidgeCV(alphas=[100., 1000., 10000.]).fit(A, ym).predict(B)
            del mdl
            torch.cuda.empty_cache()
            if not np.isfinite(p).all() or not np.isfinite(q).all():
                raise SystemExit("%s/%s: не-конечные предсказания" % (rep, CY[e]))
            P.append(p)
            Q.append(q)
            st.append(dict(n_fit=int(m.sum()), y_mean=float(ym.mean()), y_sd=float(ym.std()),
                           p_mean=float(p.mean()), p_sd=float(p.std()),
                           q_mean=float(q.mean()), q_sd=float(q.std()),
                           p_min=float(p.min()), p_max=float(p.max())))
            print("    %s %s: обучено на %d, 750 предсказаний за %.0f с"
                  % (rep, CY[e], m.sum(), time.time() - t0), flush=True)
        out[rep] = [v.tolist() for v in P]
        ridge_out[rep] = [np.asarray(v).tolist() for v in Q]
        stats[rep] = st

    print("\n  КОНТРОЛЬ разброса: 750 предсказаний против обучающих меток")
    print("  %-4s %-8s %7s %9s %9s %9s %9s %9s" % ("рек", "фермент", "n", "метки ср",
                                                   "метки sd", "пред ср", "пред sd", "доля sd"))
    bad = []
    for rep, st in stats.items():
        for e, s in enumerate(st):
            frac = s["p_sd"] / s["y_sd"]
            print("  %-4s %-8s %7d %9.3f %9.3f %9.3f %9.3f %9.2f"
                  % (rep, CY[e], s["n_fit"], s["y_mean"], s["y_sd"], s["p_mean"], s["p_sd"], frac))
            if frac < 0.2:
                bad.append("%s/%s разброс %.2f от метки" % (rep, CY[e], frac))
            if abs(s["p_mean"] - s["y_mean"]) > s["y_sd"]:
                bad.append("%s/%s среднее ушло на %.2f при sd метки %.2f"
                           % (rep, CY[e], s["p_mean"] - s["y_mean"], s["y_sd"]))
    if bad:
        print("\nОТКАЗ: " + "; ".join(bad))
        return 2
    print("  пройден: ни один арм не сжался и не ушёл по уровню")

    pathlib.Path(a.out).mkdir(parents=True, exist_ok=True)
    fn = "%s/tabicl_test.json" % a.out
    json.dump({"preds": out, "ridge": ridge_out, "pca": a.pca, "n_test": int(n_te),
               "stats": stats}, open(fn, "w"))
    print("\n  записано %s" % fn)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
