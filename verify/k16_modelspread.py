"""How far apart three models put the same test shift, and what that does to the uncertainty.

Item 87 found that delta moves by 0.11 to 0.51 when the model changes, which is as large as the
whole bootstrap interval. Two models is a thin basis for calling that systematic, so this file
adds a third and separates the three sources of uncertainty that were previously conflated:

  kernel estimation   how much delta moves when only the folds change. Test predictions do not
                      depend on the split seed, so this isolates noise in fitting E[yhat|y];
  sampling            the cluster bootstrap over both compound sets, which is what the published
                      intervals report;
  model choice        the spread across model families, which nothing has ever counted.

The third model is L1 - the same learner with absolute error instead of squared. It is the most
structurally distant of the three from either: the per-enzyme L2 model and the pooled model share
a loss and differ in what they are trained on, while L1 differs in what it estimates, the
conditional median rather than the mean. Item 75 measured it at -0.044 raw and -0.002 after
post-processing, so it is a model we would never submit - which is the point, since the question
is how much delta depends on a choice that barely changes the predictions' quality.

Every model needs blind predictions of its own, so all three are refitted on the full training
set here rather than read from data/test_pred.npz, which belongs to the first of them.

Reads data/feats.npz, the saved out-of-fold predictions, the blinded test. Prints only.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor

import feats as F
from covshift import estimate
from cypsplit import cluster_ids
from submit import pooled_design

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    dn = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
    mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    print("признаки теста", flush=True)
    FP, dsc, M, ok = F.build(list(te.SMILES), dn, mn)
    Xte = np.hstack([FP, dsc.to_numpy(np.float32), M.to_numpy(np.float32)])
    te = te.iloc[ok].reset_index(drop=True)

    OOF = {
        "поферментно L2": (json.load(open(RES + "preds/oof.json")), "FP+DESC+MECH|{c}"),
        "пул":            (json.load(open(RES + "preds/oof_pool_all.json"))["preds"], "0|пул|{c}"),
        "поферментно L1": (json.load(open(RES + "preds/oof_l1_4seed.json"))["preds"],
                           "0|L1 по метке|{c}"),
    }
    mask = np.stack([tr[f"{c}_pIC50_direct_inhibition"].notna().to_numpy() for c in CYPS], 1)
    yall = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)

    print("обучаю три модели на всей выборке", flush=True)
    PT = {}
    Xs = [pooled_design(X[mask[:, e]], e) for e in range(4)]
    ys = [yall[mask[:, e], e] for e in range(4)]
    shared = HistGradientBoostingRegressor(**KW).fit(np.vstack(Xs), np.concatenate(ys))
    for e, c in enumerate(CYPS):
        m = mask[:, e]
        PT[("поферментно L2", c)] = HistGradientBoostingRegressor(**KW).fit(
            X[m], yall[m, e]).predict(Xte)
        PT[("пул", c)] = shared.predict(pooled_design(Xte, e))
        PT[("поферментно L1", c)] = HistGradientBoostingRegressor(
            **KW, loss="absolute_error").fit(X[m], yall[m, e]).predict(Xte)

    cid = np.asarray(cluster_ids(list(rows.SMILES))[0])
    tcid = np.asarray(cluster_ids(list(te.SMILES), threshold=0.50)[0])

    print(f"\n{'фермент':8s} " + " ".join(f"{k:>16s}" for k in OOF)
          + f" {'размах':>8s} {'выборочн.':>10s} {'отнош.':>7s}")
    for e, c in enumerate(CYPS):
        m = mask[:, e]
        y = yall[m, e]
        vals, half = [], None
        for k, (src, key) in OOF.items():
            s = estimate(y, np.asarray(src[key.format(c=c)]), PT[(k, c)],
                         groups=cid[m], blind_groups=tcid, draws=400)
            vals.append(s.delta)
            if half is None:
                half = (s.hi - s.lo) / 2
        vals = np.array(vals)
        rng = vals.max() - vals.min()
        print(f"{c:8s} " + " ".join(f"{v:+16.3f}" for v in vals)
              + f" {rng:8.3f} {half:10.3f} {rng / half:7.2f}")

    print("""
Как читать. Столбец «размах» --- систематика от выбора модели, «выборочн.» --- полуширина
кластерного бутстрапа, то есть то, что документ до сих пор публиковал как всю неопределённость.

Отношение около единицы значит, что систематика сопоставима с выборочной и интервалы занижены
примерно в корень из двух. Заметно больше единицы --- что выбор модели доминирует, и тогда
публиковать бутстрапный интервал как неопределённость delta нельзя вовсе.

Оговорка, которая не снимается добавлением моделей: три --- это всё ещё выборка из семейства
градиентного бустинга на одних признаках. Настоящий размах по всем разумным моделям не меньше
измеренного, а насколько больше --- неизвестно.""")


if __name__ == "__main__":
    main()
