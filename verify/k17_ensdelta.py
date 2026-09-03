"""The shift estimated on the model that is actually submitted, and the posterior widened by three.

Where this stands. Item 87 found delta moves by 0.11 to 0.51 when the model changes and adopted
a rule under an equal-weight mixture of two posteriors. Since then the ensemble - the average of
the per-enzyme and the pooled model - turned out to beat both by 0.023 after post-processing, so
it is what gets submitted, and neither of the two posteriors is about it.

The scheme this settles on, and the reasoning for it. The point estimate should come from the
model whose predictions are being corrected, because that is whose kernel the inversion uses.
The width should come from the spread across models, because the ensemble is one model among
many and nothing guarantees its kernel invariance any more than the others'. So: estimate delta
on the ensemble, and mix its draws with those of the two components rather than replacing them.

Averaging the models would NOT be a way to resolve the disagreement, and it is worth saying why,
since the ensemble looks like exactly that. Seed-to-seed noise in the delta estimate - the folds
change, the test predictions do not - has sd 0.005 to 0.048, while the model-to-model difference
is 0.11 to 0.51, nine to eighty-eight times larger. The disagreement is systematic. Collapsing it
into a third model produces a narrow interval around an object whose own bias is unmeasured.

Writes results/preds/oof_ensemble.json, delta_draws_ens.json and delta_draws_mix3.json.
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
from submit import gbm_reg, pooled_design

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    P = json.load(open(RES + "preds/oof_pool_all.json"))["preds"]

    # Предсказания ансамбля вне фолда --- среднее двух рук на одних фолдах.
    ens = {f"{s}|ансамбль|{c}": ((np.asarray(P[f"{s}|независимо|{c}"], float)
                                 + np.asarray(P[f"{s}|пул|{c}"], float)) / 2.0).tolist()
           for s in range(4) for c in CYPS}
    json.dump({"preds": ens}, open(RES + "preds/oof_ensemble.json", "w"))
    print(f"предсказания ансамбля вне фолда сохранены: {len(ens)} ключей")

    dn = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
    mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    print("признаки теста", flush=True)
    FP, dsc, M, ok = F.build(list(te.SMILES), dn, mn)
    Xte = np.hstack([FP, dsc.to_numpy(np.float32), M.to_numpy(np.float32)])
    te = te.iloc[ok].reset_index(drop=True)

    mask = np.stack([tr[f"{c}_pIC50_direct_inhibition"].notna().to_numpy() for c in CYPS], 1)
    yall = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    print("обучаю обе базы на всей выборке", flush=True)
    Xs = [pooled_design(X[mask[:, e]], e) for e in range(4)]
    ys = [yall[mask[:, e], e] for e in range(4)]
    shared = gbm_reg().fit(np.vstack(Xs), np.concatenate(ys))
    PT = {}
    for e, c in enumerate(CYPS):
        m = mask[:, e]
        a = gbm_reg().fit(X[m], yall[m, e]).predict(Xte)
        b = shared.predict(pooled_design(Xte, e))
        PT[c] = (a + b) / 2.0

    cid = np.asarray(cluster_ids(list(rows.SMILES))[0])
    tcid = np.asarray(cluster_ids(list(te.SMILES), threshold=0.50)[0])
    di = json.load(open(RES + "preds/delta_draws.json"))
    dp = json.load(open(RES + "preds/delta_draws_pool.json"))

    DR, mix = {}, {}
    print(f"\n{'фермент':8s} {'раздельно':>10s} {'пул':>8s} {'ансамбль':>9s} | "
          f"{'размах':>7s} | {'смесь из 3':>11s} {'95%':>18s}")
    for e, c in enumerate(CYPS):
        m = mask[:, e]
        s = estimate(yall[m, e], np.asarray(ens[f"0|ансамбль|{c}"]), PT[c],
                     groups=cid[m], blind_groups=tcid, draws=1500)
        DR[c] = [float(x) for x in s.draws]
        mix[c] = list(di[c]) + list(dp[c]) + DR[c]
        v = [float(np.median(di[c])), float(np.median(dp[c])), s.delta]
        q = np.percentile(mix[c], [2.5, 97.5])
        print(f"{c:8s} {v[0]:+10.3f} {v[1]:+8.3f} {v[2]:+9.3f} | {max(v)-min(v):7.3f} | "
              f"{np.median(mix[c]):+11.3f} [{q[0]:+7.3f},{q[1]:+7.3f}]")

    json.dump(DR, open(RES + "preds/delta_draws_ens.json", "w"))
    json.dump({c: [float(x) for x in mix[c]] for c in CYPS},
              open(RES + "preds/delta_draws_mix3.json", "w"))
    print(f"\nсохранено: delta_draws_ens.json, delta_draws_mix3.json")
    print("""
Как читать. Оценка ансамбля не обязана лежать между двумя другими: обращение ядра нелинейно,
и усреднение предсказаний не есть усреднение решений уравнения. Если легла между --- расхождение
ведёт себя как гладкая функция модели; если вышла за --- ансамбль отличается от обеих баз
сильнее, чем они друг от друга, и тогда трёх моделей мало даже для нижней границы.""")


if __name__ == "__main__":
    main()
