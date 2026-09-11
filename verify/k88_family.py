"""Estimator-family screen. Item 288 (pre-registered). Sieve one of three, not an adoption.

Decorrelation in this project comes from the estimator family, not features (item 280: least
correlated pair GP-ствол at 0.913, both on DESC+MECH). This builds OOF on seed 0, per enzyme, for
families NOT in the ensemble, and measures each one's standalone rank and the correlation of its
out-of-fold error with the five-member mean's error. Passes the gate iff mean error-correlation
< 0.93 AND standalone macro rank within 0.02 of the weakest current member (item 288). Reads the
cached members_seed0.json for the ensemble mean and the current members' own error-correlations as
the calibration. Writes results/logs/k88_family.json. Minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json, time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestRegressor, ExtraTreesRegressor
from sklearn.svm import SVR
from sklearn.kernel_ridge import KernelRidge
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

from cypsplit import butina_folds
from submit import CYPS

MEM = ["поферментно", "пул", "GP", "гребневая", "ствол"]


def desc_mech(X):
    return X[:, -247:]                        # DESC+MECH block (kernel methods use it, like GP/ridge)


def oof(model_fn, X, y, fold, transform=None):
    p = np.zeros(len(y))
    for f in range(5):
        a, b = fold != f, fold == f
        if b.sum() == 0:
            continue
        Xa, Xb = (X[a], X[b]) if transform is None else (transform(X[a]), transform(X[b]))
        # clip infinities in DESC for kernel/linear stability
        Xa = np.nan_to_num(Xa, posinf=0.0, neginf=0.0)
        Xb = np.nan_to_num(Xb, posinf=0.0, neginf=0.0)
        p[b] = model_fn().fit(Xa, y[a]).predict(Xb)
    return p


CANDS = {
    "случ.лес":      (lambda: RandomForestRegressor(n_estimators=400, n_jobs=-1, random_state=0), None),
    "экстра-деревья": (lambda: ExtraTreesRegressor(n_estimators=400, n_jobs=-1, random_state=0), None),
    "SVR-RBF":       (lambda: make_pipeline(StandardScaler(), SVR(C=5.0, gamma="scale")), desc_mech),
    "ядр.гребн-RBF": (lambda: make_pipeline(StandardScaler(), KernelRidge(alpha=1.0, kernel="rbf",
                                                                         gamma=None)), desc_mech),
}


def main():
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(Y)
    fold, _ = butina_folds(list(rows.SMILES))
    mem = json.load(open(RES + "preds/members_seed0.json"))

    out = {"по_ферментам": {}, "current_members": {}}
    # reference: five-member mean and each current member's error-corr with it
    print("Калибровка: текущие члены, корреляция ошибки со средним пятёрки, и их ранг\n", flush=True)
    ens = {}
    for e, c in enumerate(CYPS):
        m = mask[:, e]; yy = Y[m, e]
        M = np.stack([np.asarray(mem[k][e], float) for k in MEM])
        ensmean = M.mean(0); ens[c] = (yy, ensmean, X[m], fold[m])
        errE = ensmean - yy
        row = {}
        for i, k in enumerate(MEM):
            row[k] = {"rho_err": float(np.corrcoef(M[i] - yy, errE)[0, 1]),
                      "rank": float(spearmanr(yy, M[i]).statistic)}
        out["current_members"][c] = row
    for k in MEM:
        rr = np.mean([out["current_members"][c][k]["rho_err"] for c in CYPS])
        rk = np.mean([out["current_members"][c][k]["rank"] for c in CYPS])
        print(f"  {k:12s}  ρ(ошибки, ансамбль) {rr:+.3f}   ранг {rk:.4f}", flush=True)
    weakest = min(np.mean([out["current_members"][c][k]["rank"] for c in CYPS]) for k in MEM)
    ensrank = np.mean([spearmanr(*ens[c][:2][::-1]).statistic for c in CYPS])
    print(f"\n  слабейший текущий член по макро-рангу: {weakest:.4f}; ансамбль: {ensrank:.4f}\n", flush=True)

    print(f"  {'семейство':14s} {'ρ ошибки ср.':>12s} {'ранг макро':>11s} {'ниже слаб.':>10s} {'вердикт':>10s}",
          flush=True)
    for name, (fn, tf) in CANDS.items():
        t0 = time.time(); rhos = []; ranks = []
        for e, c in enumerate(CYPS):
            yy, ensmean, Xm, fm = ens[c]
            p = oof(fn, Xm, yy, fm, tf)
            rhos.append(float(np.corrcoef(p - yy, ensmean - yy)[0, 1]))
            ranks.append(float(spearmanr(yy, p).statistic))
        rho = float(np.mean(rhos)); rk = float(np.mean(ranks))
        gate = rho < 0.93 and rk >= weakest - 0.02
        out["по_ферментам"][name] = {"rho_err": rhos, "rank": ranks, "rho_mean": rho,
                                     "rank_macro": rk, "gate": bool(gate)}
        print(f"  {name:14s} {rho:12.3f} {rk:11.4f} {weakest-rk:+10.4f} "
              f"{'ПРОХОДИТ' if gate else 'закрыт':>10s}   ({time.time()-t0:.0f} c)", flush=True)
        json.dump(out, open(RES + "logs/k88_family.json", "w"), ensure_ascii=False, indent=1)

    passed = [n for n in CANDS if out["по_ферментам"][n]["gate"]]
    print(f"\n  прошли ворота (ρ<0.93 при ранге в 0.02 от слабейшего): {passed or 'ни одного'}", flush=True)
    print(f"  условие пункта 288: {'есть кандидат на плечо' if passed else 'ось закрыта скрином, без ансамблевого плеча'}",
          flush=True)


if __name__ == "__main__":
    main()
