"""Splitting the fit by binding mode, where item 184 licensed it and where it did not.

The licence, and it is per enzyme rather than general. `verify/k56_modes.py` ran the gate before
anything was built: a **small** level difference between modes (else a greedy criterion would take
the split anyway) together with a **large** transfer gap at matched training size (the divergence
greedy cannot see one step ahead).

    фермент   разница уровня   разрыв, коорд   разрыв, прочие   лицензия
    CYP2D6            -0.031          +0.063           +0.175   да
    CYP3A4            +0.089          +0.067           +0.040   да
    CYP2C9            +0.096          -0.016           +0.110   половина
    CYP1A2            +0.262          +0.006           -0.120   нет

The two halves came out inversely related across the four, which is the signature: where greedy can
see the split it does not need help. **CYP1A2 is therefore a negative control that lives in the
data** -- if the mode split helps there, the harness is wrong rather than the idea.

The mode is `[nX2]` against everything else: a pyridine-type aromatic nitrogen has its lone pair in
the ring plane and available to the haem iron, a pyrrole-type `[nX3;H1]` has it in the pi system and
does not coordinate. 3505 of 4905 molecules carry the coordinating type.

**Four arms, and the second is what makes the third readable.** The outside argument had two layers
and they predict different things:

    один                    одна модель на всё --- эталон
    +индикатор              мода одной колонкой. Если этого ДОСТАТОЧНО, то дело было в
                            РЕПРЕЗЕНТАЦИОННОЙ цене --- дизъюнкция четырёх колонок стоила
                            дереву до четырёх сплитов, и хватило её предвычислить
    по модам                отдельные модели внутри мод. Если помогает ТОЛЬКО это, дело было
                            в БЛИЗОРУКОСТИ --- индикатор не берётся в корень, потому что его
                            немедленный выигрыш мал, и никакая колонка этого не чинит
    по модам, перемешанным  моды переставлены случайно, размеры те же. Отделяет "деление
                            помогает" от "деление ИМЕННО по этой моде помогает"

The honest cost is visible in the third arm and is not hidden: splitting means the minority model
sees 278 to 707 rows instead of the full table. The gate measured transfer at matched sizes; this
measures whether the divergence is worth the rows it costs.

Read per enzyme against item 165's floors, four seeds. Writes results/preds/oof_mode.json.
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
from rdkit import Chem, RDLogger
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
NTREE, LR, DEPTH, MAXFEAT = 200, 0.06, 5, 0.3      # пины из src/ablpairloss.py
FLOOR = {"CYP1A2": 0.0061, "CYP2C9": 0.0071, "CYP2D6": 0.0049, "CYP3A4": 0.0033}
LICENCE = {"CYP1A2": "нет", "CYP2C9": "половина", "CYP2D6": "ДА", "CYP3A4": "ДА"}
COORD = Chem.MolFromSmarts("[nX2]")
MIN_ROWS = 80          # меньше --- отдельную модель не строим, откатываемся к общей


def boost(Xtr, ytr, Xte, seed):
    rng = np.random.default_rng(seed)
    s = np.full(len(ytr), float(ytr.mean())); p = np.full(len(Xte), s[0])
    for _ in range(NTREE):
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MAXFEAT,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(Xtr, ytr - s)
        s = s + LR * t.predict(Xtr); p = p + LR * t.predict(Xte)
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--arms", default="один|+индикатор|по модам|по модам, перемешанным")
    ap.add_argument("--out", default=RES + "preds/oof_mode.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    arms = a.arms.split("|")

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)
    mode = np.array([bool(m and m.HasSubstructMatch(COORD))
                     for m in (Chem.MolFromSmiles(s) for s in rows.SMILES)])
    print(f"координирующих [nX2]: {mode.sum()} из {len(rows)} ({100*mode.mean():.1f} %)")
    print(f"{'фермент':8s} {'коорд':>7s} {'прочие':>7s}   лицензия по пункту 184")
    for c in CYPS:
        m = tr[f"{c}_pIC50_direct_inhibition"].notna().to_numpy()
        print(f"{c:8s} {int((m&mode).sum()):7d} {int((m&~mode).sum()):7d}   {LICENCE[c]}")
    print()

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for arm in arms:
            t0 = time.time()
            r = {"seed": seed, "рука": arm}
            # Перемешанная мода: тот же размер, разорвана связь с химией. Фиксируется на сид,
            # а не на фолд, иначе она перестала бы быть свойством молекулы.
            if "перемешанным" in arm:
                mm = mode[np.random.default_rng(seed * 31 + 5).permutation(len(mode))]
            else:
                mm = mode
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr[col].to_numpy(float)
                lo = tr[col + "_conf_low"].to_numpy(float)
                hi = tr[col + "_conf_high"].to_numpy(float)
                P = np.zeros(len(rows))
                for f in range(5):
                    trn_mol, te_mol = fold != f, fold == f
                    if arm == "один":
                        A, B = m & trn_mol, m & te_mol
                        if B.any():
                            P[B] = boost(X[A], y[A], X[B], seed * 10 + f)
                    elif arm == "+индикатор":
                        Xi = np.hstack([X, mm.astype(np.float32)[:, None]])
                        A, B = m & trn_mol, m & te_mol
                        if B.any():
                            P[B] = boost(Xi[A], y[A], Xi[B], seed * 10 + f)
                    else:
                        for g in (True, False):
                            A, B = m & trn_mol & (mm == g), m & te_mol & (mm == g)
                            if not B.any():
                                continue
                            if A.sum() < MIN_ROWS:      # мало строк --- общая модель
                                A = m & trn_mol
                            P[B] = boost(X[A], y[A], X[B], seed * 10 + f)
                yy, fi = y[m], fold[m]
                p = P[m]
                q = fit_apply(p, lo[m], hi[m], fi, np.ones(int(m.sum())) / int(m.sum()))
                out[f"{seed}|{arm}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(yy, q, y_true_upper=hi[m],
                                                   y_true_lower=lo[m])), 4)
                r[f"{c} rho"] = round(float(spearmanr(yy, p).statistic), 4)
            for t in ("пара", "rho"):
                r[f"MACRO {t}"] = round(float(np.mean([r[f"{c} {t}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {arm:24s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    g = df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                       + [f"{c} rho" for c in CYPS]].mean()
    print("\n" + g.round(4).to_string())
    if "один" in g.index:
        print(f"\n{'рука':24s} " + "".join(f"{c[3:]+' Δ':>13s}" for c in CYPS))
        for arm in g.index:
            if arm == "один":
                continue
            line = f"{arm:24s} "
            for c in CYPS:
                d = g.loc[arm, f"{c} rho"] - g.loc["один", f"{c} rho"]
                line += f"{d:+11.4f}{'*' if abs(d) > FLOOR[c] else ' '}{LICENCE[c][0]}"
            print(line)
        print("  * --- больше собственного пола (165); последняя буква --- лицензия 184: "
              "Д=да, п=половина, н=нет")
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Три вопроса, и порядок обязателен.

Первый: совпадает ли выигрыш с лицензией пункта 184. Он предсказал CYP2D6 и CYP3A4, половину
CYP2C9 и ОТСУТСТВИЕ на CYP1A2. Выигрыш на CYP1A2 означает сломанный стенд, а не находку ---
там отрицательный контроль живёт в самих данных.

Второй: «+индикатор» против «по модам». Хватает индикатора --- дело было в репрезентационной
цене, дизъюнкцию достаточно предвычислить. Помогает только деление --- дело в близорукости,
и никакая колонка этого не чинит, потому что жадность не возьмёт в корень сплит с малым
немедленным выигрышем.

Третий: перемешанная мода. Она обязана НЕ помочь, иначе помогает само деление выборки, а не
эта химия.""")


if __name__ == "__main__":
    main()
