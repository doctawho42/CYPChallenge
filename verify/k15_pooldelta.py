"""Re-estimate the test shift on the pooled base, which item 84 left as the open question.

Why it was not settled by re-running the rule. src/shrinkchoice.py takes the posterior over
delta as a fixed input from results/preds/delta_draws.json, and those draws were produced by
inverting the kernel E[yhat|y] built on the PER-ENZYME model. Feeding pooled predictions to the
chooser therefore changes which rule is optimal - it turned out not to, only CYP1A2 moves from
+0.1 to 0.0 and that enzyme is skipped anyway - but leaves the estimate of delta itself untouched.

The prediction this tests, recorded when the chooser was run. Inversion amplifies an error in the
observed blind mean by the reciprocal of the kernel slope, which was 3.4 / 2.4 / 4.5 / 1.7. A
model with better rank has a steeper kernel, so pooling should shrink those factors, most on
CYP2D6 where it gained the most rank (+0.042) and where the interval was the one that straddled
zero. If the factors do not move, the rank gain does not reach the inversion and the delta work
is unaffected by the model - which would itself be worth knowing, since it is the assumption
under which the published intervals stay valid.

Both sides need predictions from the SAME model, so the blind predictions are recomputed here
with the pooled learner rather than read from data/test_pred.npz, which came from the per-enzyme
one. Everything else is src/covshift.py unchanged.

Reads data/feats.npz, results/preds/oof_pool_all.json, the blinded test. Prints only.
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
PUB = {"CYP1A2": 0.045, "CYP2C9": 0.362, "CYP2D6": -0.917, "CYP3A4": 0.740}


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    oofp = json.load(open(RES + "preds/oof_pool_all.json"))["preds"]
    oofi = json.load(open(RES + "preds/oof.json"))

    dn = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
    mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    print("признаки теста", flush=True)
    FP, dsc, M, ok = F.build(list(te.SMILES), dn, mn)
    Xte = np.hstack([FP, dsc.to_numpy(np.float32), M.to_numpy(np.float32)])
    te = te.iloc[ok].reset_index(drop=True)
    if Xte.shape[1] != X.shape[1]:
        raise SystemExit(f"ширина {Xte.shape[1]} против {X.shape[1]}")

    mask = np.stack([tr[f"{c}_pIC50_direct_inhibition"].notna().to_numpy() for c in CYPS], 1)
    yall = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    print("обучаю пул на всей выборке и предсказываю тест", flush=True)
    Xs = [pooled_design(X[mask[:, e]], e) for e in range(4)]
    ys = [yall[mask[:, e], e] for e in range(4)]
    shared = gbm_reg().fit(np.vstack(Xs), np.concatenate(ys))
    PT = {c: shared.predict(pooled_design(Xte, e)) for e, c in enumerate(CYPS)}

    cid = np.asarray(cluster_ids(list(rows.SMILES))[0])
    tcid = np.asarray(cluster_ids(list(te.SMILES), threshold=0.50)[0])
    PTi = np.load(D + "test_pred.npz")

    DR = {}
    print(f"\n{'фермент':8s} | {'раздельно':>28s} | {'пул':>28s}")
    print(f"{'':8s} | {'delta':>7s} {'95%':>18s} {'усил':>4s} | "
          f"{'delta':>7s} {'95%':>18s} {'усил':>4s}")
    for e, c in enumerate(CYPS):
        m = mask[:, e]
        y = yall[m, e]
        g = cid[m]
        si = estimate(y, np.asarray(oofi[f"FP+DESC+MECH|{c}"]), PTi[c],
                      groups=g, blind_groups=tcid, draws=600)
        sp = estimate(y, np.asarray(oofp[f"0|пул|{c}"]), PT[c],
                      groups=g, blind_groups=tcid, draws=1500)
        DR[c] = sp.draws
        print(f"{c:8s} | {si.delta:+7.3f} [{si.lo:+7.3f},{si.hi:+7.3f}] {si.amplification():4.1f} | "
              f"{sp.delta:+7.3f} [{sp.lo:+7.3f},{sp.hi:+7.3f}] {sp.amplification():4.1f}")

    json.dump({c: [float(x) for x in DR[c]] for c in CYPS},
              open(RES + "preds/delta_draws_pool.json", "w"))
    print(f"\nрозыгрыши на пулированном ядре сохранены: {RES}preds/delta_draws_pool.json")
    print("""
Как читать. Левая половина обязана воспроизвести опубликованное (+0.045 / +0.362 / -0.917 /
+0.740), иначе правую сравнивать не с чем.

Столбец «усил» --- во сколько раз обращение ядра множит ошибку в наблюдаемом среднем по тесту.
Предсказание: у пула он должен упасть, сильнее всего на CYP2D6. Если не упал --- прирост ранга
до обращения не доходит, и опубликованные интервалы по delta остаются в силе при любой модели,
что тоже стоит знать.""")


if __name__ == "__main__":
    main()
