"""Negative correlation learning: train the members to be wrong differently, not just to be right.

The defect, and it is measured twice with the same shape. Item 176: the trunk improved by +0.031 on
its own and contributed **nothing** to the five-member ensemble. Item 182: the per-enzyme screening
arm improved by +0.029 on its own and contributed +0.0058 to the plain ensemble and +0.0014 to the
best one -- with the diagnostic item 176 lacked, the correlation between the new member's errors and
the ensemble's, at **0.963, 0.935, 0.969 and 0.945**. The members are wrong on the same compounds.
An average pays for disagreement and there is almost none left to pay for.

Why this is not item 100 again, which is the first objection. Item 100 proposed choosing members by
residual decorrelation **measured before inclusion**, and found the criterion picks the wrong
candidate. That is selection, and selection cannot create diversity that is not already in the
candidate pool. Negative correlation learning **trains for it**: each member's own loss carries a
term that rewards departing from the ensemble mean, so the diversity is manufactured during fitting
rather than shopped for afterwards. Item 100's failure is a reason to try this, not a reason not to.

The method (Liu and Yao). Member m minimises

    (1/2)(f_m - y)^2 + lambda * p_m,   p_m = (f_m - F) * sum_{j != m}(f_j - F)

with F the ensemble mean. Since the deviations sum to zero, `p_m = -(f_m - F)^2`, so the negative
gradient each tree is fitted to becomes

    (y - f_m) + 2 * lambda * (f_m - F)

-- the ordinary residual plus a push away from the consensus. At lambda = 0 the members are
independent and the arm must reproduce plain bagging exactly, which is the control.

**Members are trained jointly, round by round**, because F has to be the *current* ensemble mean; a
member trained to completion against a frozen ensemble is a different and much weaker thing.

Pre-registered, and both outcomes are worth having:

    корреляция ошибок падает И ансамбль растёт  -> разнообразие было связывающим ограничением,
                                                   и пункты 176 и 182 объяснены
    корреляция падает, ансамбль НЕ растёт       -> насыщение не в разнообразии; тогда мы узнали,
                                                   где оно НЕ находится, что дороже прироста

A third outcome closes the method rather than the question: if the correlation does not fall as
lambda rises, the penalty is not reaching the fit and the implementation is wrong, not the idea.

Writes results/preds/oof_ncl.json.
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
from sklearn.tree import DecisionTreeRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
NTREE, LR, DEPTH, MAXFEAT = 200, 0.06, 5, 0.3      # пины из src/ablpairloss.py
NMEM = 4
SUBROW = 0.7      # доля строк на члена: без неё члены стартуют почти одинаковыми
FLOOR = {"CYP1A2": 0.0061, "CYP2C9": 0.0071, "CYP2D6": 0.0049, "CYP3A4": 0.0033}


def ncl(Xtr, ytr, Xte, seed, lam):
    """NMEM членов, обучаемых СОВМЕСТНО. Возвращает (ансамбль, матрица членов на тесте).

    Каждому члену --- своя подвыборка СТРОК, фиксированная на всё обучение. Без неё дымовой
    прогон дал корреляцию ошибок 0.9952 при lam = 0 и 0.9950 при lam = 0.25: штраф добавляет
    к остатку 2*lam*(f_m - F), и если члены почти совпадают, эта добавка равна нулю. NCL ---
    положительная обратная связь, которой нужна начальная асимметрия, чтобы было что
    усиливать. Одного random_state при подвыборке колонок для этого мало.
    """
    rngs = [np.random.default_rng(seed * 100 + m) for m in range(NMEM)]
    subs = [rngs[m].choice(len(ytr), max(20, int(SUBROW * len(ytr))), replace=False)
            for m in range(NMEM)]
    s = np.tile(float(ytr.mean()), (NMEM, len(ytr)))
    p = np.tile(float(ytr.mean()), (NMEM, len(Xte)))
    for _ in range(NTREE):
        F = s.mean(0)                       # текущее среднее ансамбля
        for m in range(NMEM):
            # обычный остаток плюс отталкивание от консенсуса
            r = (ytr - s[m]) + 2.0 * lam * (s[m] - F)
            t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MAXFEAT,
                                      random_state=int(rngs[m].integers(1 << 30)))
            t.fit(Xtr[subs[m]], r[subs[m]])
            s[m] += LR * t.predict(Xtr)
            p[m] += LR * t.predict(Xte)
    return p.mean(0), p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1")
    ap.add_argument("--lams", default="0,0.25,0.5,0.75")
    ap.add_argument("--out", default=RES + "preds/oof_ncl.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    lams = [float(x) for x in a.lams.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)
    print(f"членов {NMEM}, деревьев {NTREE}, штрафы {lams}")
    print("цель: корреляция ошибок членов; пункты 176 и 182 намерили 0.94-0.97 "
          "между новым членом и ансамблем\n")

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for lam in lams:
            t0 = time.time()
            r = {"seed": seed, "lam": lam}
            corrs, mem_rho = [], []
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                Xi, fi = X[m], fold[m]
                p = np.zeros(len(y)); PM = np.zeros((NMEM, len(y)))
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    pe, pm = ncl(Xi[trn], y[trn], Xi[te], seed * 10 + f, lam)
                    p[te] = pe; PM[:, te] = pm
                # средняя ПОПАРНАЯ корреляция ошибок членов --- диагностика штрафа
                E = PM - y
                cs = [np.corrcoef(E[i], E[j])[0, 1]
                      for i in range(NMEM) for j in range(i + 1, NMEM)]
                corrs.append(float(np.mean(cs)))
                mem_rho.append(float(np.mean([spearmanr(y, PM[i]).statistic
                                              for i in range(NMEM)])))
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{lam}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for t in ("пара", "rho"):
                r[f"MACRO {t}"] = round(float(np.mean([r[f"{c} {t}"] for c in CYPS])), 4)
            r["корр ошибок"] = round(float(np.mean(corrs)), 4)
            r["ранг члена"] = round(float(np.mean(mem_rho)), 4)
            table.append(r)
            print(f"  сид {seed} lam {lam:<5} пара {r['MACRO пара']:.4f} "
                  f"ансамбль rho {r['MACRO rho']:.4f}  член rho {r['ранг члена']:.4f}  "
                  f"корр ошибок {r['корр ошибок']:.4f}  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    g = df.groupby("lam")[["MACRO пара", "MACRO rho", "ранг члена", "корр ошибок"]].mean()
    print("\n" + g.round(4).to_string())
    if 0.0 in g.index:
        print(f"\n{'lam':>6s} {'Δ ансамбль':>12s} {'Δ член':>10s} {'Δ корреляция':>14s}")
        for lam in g.index:
            print(f"{lam:6} {g.loc[lam,'MACRO rho']-g.loc[0.0,'MACRO rho']:+12.4f} "
                  f"{g.loc[lam,'ранг члена']-g.loc[0.0,'ранг члена']:+10.4f} "
                  f"{g.loc[lam,'корр ошибок']-g.loc[0.0,'корр ошибок']:+14.4f}")
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Сначала столбец «корр ошибок»: если он не падает с ростом lam, штраф не доходит
до подгонки, и закрыта реализация, а не идея. Если падает --- смотреть на «Δ ансамбль».

Рост ансамбля при падающей корреляции означает, что разнообразие и было связывающим
ограничением, и пункты 176 и 182 получают объяснение. Отсутствие роста означает, что
насыщение не в разнообразии --- и это дороже прироста, потому что закрывает целое направление
поиска членов.

Столбец «ранг члена» --- цена: NCL всегда ухудшает отдельного члена, вопрос лишь в том,
окупается ли это в среднем.""")


if __name__ == "__main__":
    main()
