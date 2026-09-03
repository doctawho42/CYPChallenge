"""A Free-Wilson member: additive in fragment indicators, and stratified against the GP's profile.

Why this repository is the case Free-Wilson was invented for, and why it is nonetheless not the
classical method.

**The test set is 132 congeneric series** (item 133's construction reads across "the members of each
of the 132 test series"). Free-Wilson -- activity equals a scaffold contribution plus a sum of
substituent contributions, fitted by regression on fragment indicators -- was devised for exactly
that data shape.

**The basis already exists and this model found it.** Item 150 measured that fifty fingerprint bits
carry 80.7 per cent of the whole block's effect and item 156 decoded them: pyrimidine C2, N-aryl
azole, pyridine nitrogen, azine carbon, tertiary aliphatic amine. That is a Free-Wilson basis
obtained by boosting rather than by a chemist.

**And it evades what killed the series layer.** Substituent contributions are estimated across the
whole training set -- a fragment appears in hundreds of molecules, not only inside one series -- so
"there are no series in training" is not an objection here. It is, however, why this is not the
classical method: with 4520 Murcko scaffolds over 4905 molecules there is no scaffold term to fit,
so the intercept absorbs it and the model is additive in fragments alone.

**The honest objection, and its measurable form.** The boosting already learns piecewise-constant
fragment contributions; Free-Wilson adds additivity by construction. Additivity holds within a
congeneric series and breaks between scaffolds, so cluster cross-validation puts it in its worst
regime and the test in its best. That is exactly the shape of a convenient excuse, so it is measured
rather than asserted:

    вклад аддитивного члена по стратам сходства ближайшего соседа
    растёт с близостью соседа  ->  аддитивность работает там, где сидит тест
    не растёт                  ->  отговорка, закрыто

The prediction is **opposite** to the GP's measured profile: item 97 found the GP's contribution is
-0.0161 in the far stratum and -0.0084 in the middle, that is it helps most where neighbours are
distant. Two members with opposite stratum profiles would be a real complementarity rather than one
more correlated member -- and items 176 and 182 measured that correlated members are exactly what
this ensemble no longer benefits from.

Arms. The fourth exists because the third would otherwise be unreadable: if all 2048 bits do as well
as the top fifty, the selection is doing nothing and the claim about a compact basis is empty.

    гребневая, дескрипторы     нынешний член ансамбля, 247 колонок DESC+MECH
    Фри-Вилсон, top-50         аддитивная модель на пятидесяти битах
    Фри-Вилсон, top-200
    Фри-Вилсон, все 2048

Bits are ranked **on training folds only**, by the same split-count procedure as
`verify/k40_topk.py`, so the basis is never selected on the answer.

Writes results/preds/oof_fw.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.linear_model import RidgeCV
from sklearn.tree import DecisionTreeRegressor
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
NTREE, LR, DEPTH, MF = 150, 0.06, 5, 0.3          # для ранжирования битов, как в k40
ALPHAS = np.logspace(-2, 4, 15)
STRATA = [(0.0, 0.35), (0.35, 0.45), (0.45, 0.55), (0.55, 1.01)]


def bit_counts(X, y, nfp, seed):
    """Сколько раз бит выбран сплитом. Копия процедуры verify/k40_topk.py."""
    rng = np.random.default_rng(seed)
    s = np.full(len(y), float(y.mean()))
    cnt = np.zeros(nfp, int)
    for _ in range(NTREE):
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MF,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(X, y - s)
        s = s + LR * t.predict(X)
        f = t.tree_.feature
        f = f[(f >= 0) & (f < nfp)]
        np.add.at(cnt, f, 1)
    return cnt


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/oof_fw.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    FP = z["FP"].astype(np.float64)
    DM = np.hstack([z["DESC"], z["MECH"]]).astype(np.float64)
    nfp = FP.shape[1]

    # Сходство с ближайшим соседом ЧЕРЕЗ фолд --- страты пункта 97.
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]

    def scaled(A):
        A = np.nan_to_num(A, posinf=0.0, neginf=0.0)
        return np.clip((A - A.mean(0)) / (A.std(0) + 1e-9), -5.0, 5.0)

    arms = ["гребневая, дескрипторы", "Фри-Вилсон, top-50",
            "Фри-Вилсон, top-200", "Фри-Вилсон, все 2048"]
    out, table, strat = {}, [], []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        nn = np.zeros(len(rows))
        for f in range(5):
            te, trn = np.where(fold == f)[0], np.where(fold != f)[0]
            ref = [fps[j] for j in trn]
            for i in te:
                nn[i] = max(DataStructs.BulkTanimotoSimilarity(fps[i], ref))
        for arm in arms:
            t0 = time.time()
            r = {"seed": seed, "рука": arm}
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                fi = fold[m]
                p = np.zeros(len(y))
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    if arm.startswith("гребнев"):
                        A, B = scaled(DM[m][trn]), scaled(DM[m][te])
                    else:
                        if "все" in arm:
                            keep = np.arange(nfp)
                        else:
                            k = int(arm.split("top-")[1])
                            cnt = bit_counts(np.hstack([FP[m][trn], DM[m][trn]]),
                                             y[trn], nfp, seed * 10 + f)
                            keep = np.argsort(-cnt)[:k]
                        A, B = FP[m][trn][:, keep], FP[m][te][:, keep]
                    p[te] = RidgeCV(alphas=ALPHAS).fit(A, y[trn]).predict(B)
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{arm}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
                for j, (a0, b0) in enumerate(STRATA):
                    sel = (nn[m] >= a0) & (nn[m] < b0)
                    if sel.sum() >= 40:
                        strat.append(dict(seed=seed, рука=arm, фермент=c, страта=j,
                                          n=int(sel.sum()),
                                          rho=float(spearmanr(y[sel], p[sel]).statistic)))
            for t in ("пара", "rho"):
                r[f"MACRO {t}"] = round(float(np.mean([r[f"{c} {t}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {arm:24s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}"
                  f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print("\n" + df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
          + [f"{c} rho" for c in CYPS]].mean().round(4).to_string())

    ds = pd.DataFrame(strat)
    piv = ds.pivot_table(index="страта", columns="рука", values="rho")
    cnt = ds[ds.рука == arms[0]].groupby("страта").n.sum()
    piv.index = [f"{STRATA[i][0]:.2f}-{STRATA[i][1]:.2f}  n={cnt.get(i, 0)}" for i in piv.index]
    print("\n=== ранг ПО СТРАТАМ сходства ближайшего соседа ===")
    print(piv.round(4).to_string())
    print("\n=== разность Фри-Вилсон минус гребневая, по стратам ===")
    for a_ in arms[1:]:
        if a_ in piv.columns:
            d = piv[a_] - piv[arms[0]]
            print(f"  {a_:22s} " + "  ".join(f"{v:+.4f}" for v in d.to_numpy()))
    json.dump({"table": table, "preds": out, "strat": strat}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Уровень аддитивной модели сам по себе не важен --- она заведомо слабее бустинга,
как и гауссов процесс, который слабее и всё равно полезен в ансамбле (пункт 92).

Решает ПРОФИЛЬ ПО СТРАТАМ. Предсказано: разность растёт слева направо, то есть аддитивность
работает тем лучше, чем ближе сосед --- а тест сидит на медианном сходстве 0.587 против
0.435 у отложенных строк. Не растёт --- отговорка закрыта.

И профиль обязан быть ОБРАТНЫМ к профилю GP из пункта 97, где вклад больше всего в дальней
страте. Два члена с противоположным профилем --- настоящая дополнительность, а не ещё один
коррелированный член, чего пункты 176 и 182 намерили как бесполезное.""")


if __name__ == "__main__":
    main()
