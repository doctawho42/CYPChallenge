"""The external rows again, with a source indicator instead of a scale correction.

The reason a constant shift was the wrong instrument. The gap between the two label sets was
treated as a bias to be subtracted, and three subtractions were tried: none, the offset
measured on shared compounds, the offset that equalises the marginals. The crudest of them
was actively harmful. But the gap is not a bias, it is a SELECTION effect - ChEMBL contains
what people chose to publish, which is what worked - and selection does not act uniformly
across chemical space. Subtracting one number per enzyme is the same mistake this document
already diagnosed for the test set in verify/k9: a shift is not the whole shape.

What this arm does instead. One extra column, one for external rows and zero for ours. Train
on the union, predict with the column set to zero. The learner then decides for itself, region
by region, how much the external labels say about ours, instead of being told a single number
that has to hold everywhere.

The pre-registration matters here, because this arm is a strict generalisation of both
extremes. Split on the indicator at the root and the model reproduces "no external rows";
never split on it and the model reproduces "external as is". So the result is expected to
land BETWEEN those two, and the only informative outcome is if it lands below both. That
would mean the borrowing genuinely has to differ across regions, which no single offset can
express. Landing between them means the indicator bought nothing that a constant did not.

The second arm follows the same premise further than an indicator does. If the gap is
selection on the label, the textbook correction is not on the label at all but on the WEIGHT:
the external sample is ours reweighted by whatever made a result publishable, and the repair
is to divide that back out. Exponential tilting is the minimum-relative-entropy way to move a
mean, it is already written in src/reweight.py for the delta work, and the learner takes
sample weights. So external rows are weighted by exp(theta*y) with theta chosen to bring their
weighted label mean onto ours, and no label is altered.

That is the difference the marginal-shift arm got wrong and paid +0.037 for. Shifting corrupts
every label, including the ones that were right. Tilting leaves them all alone and only
down-weights the actives that are over-represented because actives get published. If the gap
really is selection, tilting should beat both the shift and the indicator; if it is genuine
calibration drift between assays, it should do nothing, because reweighting cannot move a
label that is simply on a different scale.

Everything else is src/ablext.py unchanged: same learner, same Butina folds, same masks, same
metric, external rows in the training folds only.

Reads data/feats.npz, data/rows.csv and the external CSVs. Writes results/preds/oof_src.json.
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
from rdkit import Chem
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

import feats as F
from cypsplit import butina_folds
from reweight import tilt

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
# Сдвиг по парным соединениям (src/ablext.py, n = 15 / 6 / 41 / 9). Тонкая оценка, но
# независимая, и она не наследует обогащение активными, которым раздут маргинальный разрыв.
PAIRED = {"CYP1A2": 0.22, "CYP2C9": 0.35, "CYP2D6": 0.60, "CYP3A4": 0.56}
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)


def canon(s):
    try:
        m = Chem.MolFromSmiles(s)
        return Chem.MolToSmiles(m) if m else None
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ext-x", required=True)
    ap.add_argument("--ext-y", required=True)
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--arms", default="индикатор,наклон,наклон-парный")
    ap.add_argument("--cache", default="")
    ap.add_argument("--out", default=RES + "preds/oof_src.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])

    ext = pd.concat([pd.read_csv(a.ext_x), pd.read_csv(a.ext_y)], axis=1)
    ext["k"] = [canon(s) for s in ext.OPENADMET_CANONICAL_SMILES]
    ours = set(filter(None, (canon(s) for s in rows.SMILES)))
    te_keys = set(filter(None, (canon(s) for s in
                                pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv").SMILES)))
    if ext.k.isin(te_keys).any():
        raise SystemExit("внешние данные пересекаются с тестом")
    ext = ext[~ext.k.isin(ours)].reset_index(drop=True)

    if a.cache and _pl.Path(a.cache).exists():
        Xe = np.load(a.cache)["X"]
        if len(Xe) != len(ext):
            raise SystemExit(f"кэш на {len(Xe)} строк, а внешних {len(ext)}")
    else:
        dn = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
        mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
        print(f"считаю признаки для {len(ext)} внешних соединений", flush=True)
        FP, dsc, M, ok = F.build(list(ext.OPENADMET_CANONICAL_SMILES), dn, mn)
        Xe = np.hstack([FP, dsc.to_numpy(np.float32), M.to_numpy(np.float32)])
        ext = ext.iloc[ok].reset_index(drop=True)
        if a.cache:
            np.savez_compressed(a.cache, X=Xe)
    if Xe.shape[1] != X.shape[1]:
        raise SystemExit(f"ширина не совпала: {Xe.shape[1]} против {X.shape[1]}")

    # Последний столбец --- источник. Обучаемся на обоих значениях, предсказываем при нуле.
    Xi_all = np.hstack([X, np.zeros((len(X), 1), np.float32)])
    Xe_all = np.hstack([Xe, np.ones((len(Xe), 1), np.float32)])
    print(f"строк: наших {len(X)}, внешних {len(Xe)}; ширина {Xi_all.shape[1]}\n")

    arms = a.arms.split(",")
    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for arm in arms:
            t0 = time.time()
            r = {"seed": seed, "рука": arm}
            note = []
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                fi = fold[m]
                ye = ext[f"OPENADMET_LOGAC50_{c.lower()}"].to_numpy(float)
                em = ~np.isnan(ye)

                if arm == "индикатор":
                    Xm, Xa, w = Xi_all[m], Xe_all[em], None
                else:
                    # Наклоняем внешнюю метку, ничего в ней не меняя. Цель --- либо наше
                    # среднее (маргинальная), либо сдвиг по парным соединениям.
                    Xm, Xa = X[m], Xe[em]
                    d = (-PAIRED[c] if arm.startswith("наклон-парный")
                         else float(y.mean() - ye[em].mean()))
                    we = tilt(ye[em], d)
                    if arm.endswith("перемешанный"):
                        # Контроль: тот же эффективный размер выборки, связь с меткой
                        # разрушена. Если он повторяет результат наклона, выигрывала
                        # не форма весов, а просто уменьшение количества данных.
                        we = np.random.default_rng(seed).permutation(we)
                    w = np.concatenate([np.ones(int(m.sum())), we])
                    note.append(f"{c} наклон {d:+.2f}, эфф.размер "
                                f"{we.sum() ** 2 / (we ** 2).sum():.0f} из {int(em.sum())}")

                p = np.zeros_like(y)
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    Xt = np.vstack([Xm[trn], Xa])
                    yt = np.concatenate([y[trn], ye[em]])
                    kw = {}
                    if w is not None:
                        kw["sample_weight"] = np.concatenate([w[:len(y)][trn], w[len(y):]])
                    # При индикаторе предсказываем с источник = наш: столбец в Xm уже нулевой.
                    p[te] = HistGradientBoostingRegressor(**KW).fit(Xt, yt, **kw).predict(Xm[te])
                out[f"{seed}|{arm}|{c}"] = p.tolist()
                r[c] = round(float(strae(y, p, y_true_upper=hi, y_true_lower=lo)), 4)
            r["MACRO"] = round(float(np.mean([r[c] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {arm:12s} макро {r['MACRO']:.4f}  ({time.time()-t0:.0f} с)",
                  flush=True)
            for ln in note:
                print(f"      {ln}")

    print()
    print(pd.DataFrame(table)[["seed", "рука", *CYPS, "MACRO"]].to_string(index=False))
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Сравнивать надо с ДВУМЯ числами сразу, потому что индикатор обобщает обе крайности. На сиде 0
без внешних строк выходит 0.7673, с внешними как есть --- 0.7545. Попадание между ними
означает, что индикатор не купил ничего сверх константы. Ниже обоих --- что заимствование в
разных областях пространства должно быть разным, и никакой один сдвиг этого не выражает.

Рука "наклон-перемешанный" --- обязательный контроль, а не вариант. Наклон уменьшает
эффективный размер внешней выборки (на CYP3A4 с 4773 до 106), и выигрыш мог бы объясняться
просто тем, что внешних данных стало меньше. Перестановка тех же весов по строкам сохраняет
эффективный размер точно и разрушает связь с меткой. Повторит результат наклона --- работала
не форма, а количество, и весь довод про отбор отменяется.

Наклон отвечает на другой вопрос: разрыв шкал --- это отбор или калибровка. Побьёт "как
есть" --- отбор, и его надо было снимать весом, а не вычитанием. Не сдвинет ничего ---
калибровка, и вес тут бессилен по построению: перевзвешивание не двигает метку, лежащую на
другой шкале. Эффективный размер выборки печатается рядом; если он проседает вдвое, часть
проигрыша --- просто потеря данных, а не свойство поправки.""")


if __name__ == "__main__":
    main()
