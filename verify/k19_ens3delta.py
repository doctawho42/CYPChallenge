"""The shift on the three-member ensemble, and whether a different model family widens the spread.

Two questions in one run. The submitted model changed again - the Gaussian process joined the
ensemble in item 92 - so its own delta has to be estimated, since item 87 established that the
point estimate must come from the model whose predictions are being corrected.

And the spread across models is now testable against a real objection. Items 87 and 88 measured
it on three models that were all gradient boosting on the same features, so the growth from two
models to three could have been three variants measuring themselves. The GP is a different
family with a different failure mode. If the spread widens again when it is added, the quantity
is still not converged; if it does not, the boosting family was the source and the estimate is
closer to complete than it looked.

Four candidate models enter the mixture: the two components of the boosting ensemble, the GP,
and the three-member ensemble that is actually submitted. The two-member ensemble is left out
deliberately - it is no longer a model anyone would submit, and including it would count the
boosting pair twice.

Writes results/preds/oof_ens3.json, delta_draws_ens3.json, delta_draws_mix4.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json

import numpy as np
import pandas as pd

import feats as F
from covshift import estimate
from cypsplit import cluster_ids
from gp import prepare as gp_prepare, gp_predict
from submit import gbm_reg, pooled_design

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    P = json.load(open(RES + "preds/oof_pool_all.json"))["preds"]
    G = json.load(open(RES + "preds/oof_gp.json"))["preds"]

    ens = {f"{s}|ансамбль3|{c}": ((np.asarray(P[f"{s}|независимо|{c}"], float)
                                  + np.asarray(P[f"{s}|пул|{c}"], float)
                                  + np.asarray(G[f"{s}|GP|{c}"], float)) / 3.0).tolist()
           for s in range(4) for c in CYPS}
    json.dump({"preds": ens}, open(RES + "preds/oof_ens3.json", "w"))

    dn = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
    mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    print("признаки теста", flush=True)
    FP, dsc, M, ok = F.build(list(te.SMILES), dn, mn)
    Xte = np.hstack([FP, dsc.to_numpy(np.float32), M.to_numpy(np.float32)])
    te = te.iloc[ok].reset_index(drop=True)

    mask = np.stack([tr[f"{c}_pIC50_direct_inhibition"].notna().to_numpy() for c in CYPS], 1)
    yall = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    print("обучаю три базы на всей выборке", flush=True)
    Xs = [pooled_design(X[mask[:, e]], e) for e in range(4)]
    ys = [yall[mask[:, e], e] for e in range(4)]
    shared = gbm_reg().fit(np.vstack(Xs), np.concatenate(ys))
    PT = {}
    for e, c in enumerate(CYPS):
        m = mask[:, e]
        a = gbm_reg().fit(X[m], yall[m, e]).predict(Xte)
        b = shared.predict(pooled_design(Xte, e))
        tf = gp_prepare(X[m])
        g = gp_predict(tf(X[m]), yall[m, e], tf(Xte))
        PT[("GP", c)] = g
        PT[("ансамбль3", c)] = (a + b + g) / 3.0

    cid = np.asarray(cluster_ids(list(rows.SMILES))[0])
    tcid = np.asarray(cluster_ids(list(te.SMILES), threshold=0.50)[0])
    di = json.load(open(RES + "preds/delta_draws.json"))
    dp = json.load(open(RES + "preds/delta_draws_pool.json"))

    DR, DG, mix = {}, {}, {}
    print(f"\n{'фермент':8s} {'раздельно':>10s} {'пул':>8s} {'GP':>8s} {'ансамбль3':>10s} | "
          f"{'размах 4':>9s} {'было на 3':>10s} | {'смесь':>8s}")
    WAS = {"CYP1A2": 0.354, "CYP2C9": 0.518, "CYP2D6": 0.762, "CYP3A4": 0.230}
    for e, c in enumerate(CYPS):
        m = mask[:, e]
        y = yall[m, e]
        sg = estimate(y, np.asarray(G[f"0|GP|{c}"]), PT[("GP", c)],
                      groups=cid[m], blind_groups=tcid, draws=1500)
        se = estimate(y, np.asarray(ens[f"0|ансамбль3|{c}"]), PT[("ансамбль3", c)],
                      groups=cid[m], blind_groups=tcid, draws=1500)
        DG[c] = [float(x) for x in sg.draws]
        DR[c] = [float(x) for x in se.draws]
        mix[c] = list(di[c]) + list(dp[c]) + DG[c] + DR[c]
        v = [float(np.median(di[c])), float(np.median(dp[c])), sg.delta, se.delta]
        print(f"{c:8s} {v[0]:+10.3f} {v[1]:+8.3f} {v[2]:+8.3f} {v[3]:+10.3f} | "
              f"{max(v)-min(v):9.3f} {WAS[c]:10.3f} | {np.median(mix[c]):+8.3f}")

    json.dump(DG, open(RES + "preds/delta_draws_gp.json", "w"))
    json.dump(DR, open(RES + "preds/delta_draws_ens3.json", "w"))
    json.dump({c: [float(x) for x in mix[c]] for c in CYPS},
              open(RES + "preds/delta_draws_mix4.json", "w"))
    print("\nсохранено: delta_draws_gp.json, delta_draws_ens3.json, delta_draws_mix4.json")
    print("""
Как читать. Столбец «было на 3» --- размах по трём бустинговым моделям из пункта 88. Если
четвёртая, из другого семейства, его заметно раздвинула, значит прежние три мерили сами себя и
систематика по-прежнему не сошлась. Если не раздвинула --- источник разброса внутри бустинга, и
оценка ближе к полной, чем казалось.""")


if __name__ == "__main__":
    main()
