"""Метрика взвешивает строки неравномерно, а половина ансамбля об этом не знает. Пункт 257.

ST-RAE делит на знаменатель СВОЕГО фермента, а макро усредняет четыре дроби, поэтому вклад
одной строки в отчётное число есть 1/(4*den_e). На обучающих фолдах строка CYP2C9 стоит в
3.03 раза дороже строки CYP3A4. Поферментные члены к этому безразличны --- монотонное
перемасштабирование одной функции потерь не двигает оптимум, --- но ПУЛИРОВАННЫЙ член
складывает все четыре набора меток в одну таблицу и учится без sample_weight.

Плечи:

    равные веса          --- нынешнее поведение, ровно submit._oof_one(pool=True)
    веса 1/den           --- вмешательство, поферментная константа, нормировка на среднее 1
    ПЕРЕСТАНОВКИ 1/den   --- нуль: тот же МУЛЬТИМНОЖЕСТВО весов, другое назначение ферментам

Нуль выбран перестановочный, а не случайный: он согласован по дисперсии ТОЧНО, а не в
среднем, и отвечает ровно на тот вопрос, который нужен, --- работает ли КОНКРЕТНОЕ назначение
весов или сам факт их неравенства. Взяты три перестановки БЕЗ НЕПОДВИЖНЫХ ТОЧЕК, то есть
такие, где ни один фермент не получает свой собственный вес.

Есть и четвёртое плечо, техническое: те же равные веса, но ПЕРЕДАННЫЕ явным sample_weight.
Оно меряется на одном сиде и отвечает на вопрос, воспроизводит ли взвешенная ветка кода
невзвешенную при единичных весах. Если нет --- часть любого «выигрыша» есть смена ветки, а не
веса. Такие подмены этот файл ловил дважды.

Знаменатель считается ПО ОБУЧАЮЩИМ строкам каждого фолда, а не по всем: den зависит от меток,
и подсчёт по всему занёс бы отложенные метки в обучающие веса.

Условия приёмки записаны в пункте 257 ДО прогона и здесь не повторяются.

Читает data/feats.npz, data/rows.csv и обучающую таблицу. Пишет results/logs/k80_denweight.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse, json, time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply
from submit import CYPS, gbm_reg, pooled_design, _oof_one

# Перестановки без неподвижных точек: каждый фермент получает ЧУЖОЙ вес. Зафиксированы
# литералом, а не выбором ГСЧ, чтобы их нельзя было переподобрать по результату.
PERMS = [(3, 2, 1, 0), (1, 0, 3, 2), (2, 3, 0, 1)]


def den(yv, lov, hiv):
    """Знаменатель ST-RAE: константный предсказатель в mean(y) через тот же мягкий порог."""
    m = yv.mean()
    return float(np.maximum(0.0, np.maximum(lov - m, m - hiv)).sum())


def pooled_oof(X, y, mask, fold, wvec=None):
    """Пулированный член, при wvec=None --- копия submit._oof_one(pool=True) строка в строку."""
    P = [np.zeros(int(mask[:, e].sum())) for e in range(len(CYPS))]
    for f in range(5):
        Xs, ys, ws = [], [], []
        for e in range(len(CYPS)):
            sel = mask[:, e] & (fold != f)
            if sel.any():
                Xs.append(pooled_design(X[sel], e))
                ys.append(y[sel, e])
                ws.append(np.full(int(sel.sum()), 1.0 if wvec is None else float(wvec[e])))
        Xtr, ytr = np.vstack(Xs), np.concatenate(ys)
        if wvec is None:
            model = gbm_reg().fit(Xtr, ytr)
        else:
            w = np.concatenate(ws)
            model = gbm_reg().fit(Xtr, ytr, sample_weight=w / w.mean())
        for e in range(len(CYPS)):
            m = mask[:, e]
            b = fold[m] == f
            if b.sum() == 0:
                continue
            P[e][b] = model.predict(pooled_design(X[m][b], e))
    return P


def denom_weights(y, LO, HI, mask, fold, f):
    """Вектор 1/den по ферментам, посчитанный на ОБУЧАЮЩИХ строках фолда f, среднее 1."""
    d = np.array([den(y[mask[:, e] & (fold != f), e],
                      LO[mask[:, e] & (fold != f), e],
                      HI[mask[:, e] & (fold != f), e]) for e in range(len(CYPS))])
    v = 1.0 / d
    return v / v.mean(), d


def pooled_oof_perfold(X, y, LO, HI, mask, fold, perm=None):
    """Как pooled_oof, но веса пересчитываются внутри каждого фолда. perm=None --- порядок как есть."""
    P = [np.zeros(int(mask[:, e].sum())) for e in range(len(CYPS))]
    for f in range(5):
        v, _ = denom_weights(y, LO, HI, mask, fold, f)
        if perm is not None:
            v = v[list(perm)]
        Xs, ys, ws = [], [], []
        for e in range(len(CYPS)):
            sel = mask[:, e] & (fold != f)
            if sel.any():
                Xs.append(pooled_design(X[sel], e))
                ys.append(y[sel, e])
                ws.append(np.full(int(sel.sum()), float(v[e])))
        w = np.concatenate(ws)
        model = gbm_reg().fit(np.vstack(Xs), np.concatenate(ys), sample_weight=w / w.mean())
        for e in range(len(CYPS)):
            m = mask[:, e]
            b = fold[m] == f
            if b.sum() == 0:
                continue
            P[e][b] = model.predict(pooled_design(X[m][b], e))
    return P


def score(P, y, LO, HI, mask, fold):
    """(макро пара, макро ранг, по ферментам). Пара --- ST-RAE ПОСЛЕ аффинной пары по фолдам."""
    pr, rk = [], []
    for e in range(len(CYPS)):
        m = mask[:, e]
        yy, lo, hi, p = y[m, e], LO[m, e], HI[m, e], P[e]
        u = np.ones(len(yy)) / len(yy)
        q = fit_apply(p, lo, hi, fold[m], u)
        pr.append(float(strae(yy, q, y_true_upper=hi, y_true_lower=lo)))
        rk.append(float(spearmanr(yy, p).statistic))
    return float(np.mean(pr)), float(np.mean(rk)), pr, rk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--perms", type=int, default=len(PERMS), help="сколько перестановок нуля")
    ap.add_argument("--branch-check", action="store_true",
                    help="плечо «равные веса явным sample_weight» --- проверка ветки кода")
    a = ap.parse_args()
    seeds = [int(s) for s in a.seeds.split(",") if s.strip()]

    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)

    fold0, _ = butina_folds(list(rows.SMILES))
    v0 = np.array([den(y[mask[:, e], e], LO[mask[:, e], e], HI[mask[:, e], e])
                   for e in range(len(CYPS))])
    print("Знаменатели ST-RAE на всей обучающей выборке и вес строки в макро\n")
    print(f"  {'фермент':10s} {'n':>6s} {'знаменатель':>13s} {'вес строки':>12s}")
    u0 = (1.0 / v0) / (1.0 / v0).mean()
    for e, c in enumerate(CYPS):
        print(f"  {c:10s} {int(mask[:, e].sum()):6d} {v0[e]:13.1f} {u0[e]:12.3f}")
    print(f"  {'':10s} {'':6s} {'размах':>13s} {u0.max() / u0.min():11.2f}x\n")

    out = {"denominators": v0.tolist(), "row_weight": u0.tolist(), "perms": PERMS[:a.perms],
           "seeds": {}}
    for s in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        res = {}
        t0 = time.time()
        P = _oof_one(X, y, mask, fold, True)
        res["равные"] = score(P, y, LO, HI, mask, fold)
        print(f"сид {s}: равные веса   пара {res['равные'][0]:.4f}  ранг {res['равные'][1]:.4f}"
              f"   ({time.time() - t0:.0f} c)")

        if a.branch_check and s == seeds[0]:
            Pb = pooled_oof(X, y, mask, fold, wvec=np.ones(len(CYPS)))
            same = all(np.array_equal(P[e], Pb[e]) for e in range(len(CYPS)))
            res["ветка"] = (score(Pb, y, LO, HI, mask, fold), bool(same))
            print(f"        ветка: единичный sample_weight воспроизводит невзвешенную "
                  f"{'ДА' if same else 'НЕТ'}  пара {res['ветка'][0][0]:.4f} "
                  f"ранг {res['ветка'][0][1]:.4f}")

        Pw = pooled_oof_perfold(X, y, LO, HI, mask, fold)
        res["1/den"] = score(Pw, y, LO, HI, mask, fold)
        print(f"        веса 1/den    пара {res['1/den'][0]:.4f}  ранг {res['1/den'][1]:.4f}"
              f"   dранг {res['1/den'][1] - res['равные'][1]:+.4f}")

        res["перестановки"] = []
        for pm in PERMS[:a.perms]:
            Pp = pooled_oof_perfold(X, y, LO, HI, mask, fold, perm=pm)
            sc = score(Pp, y, LO, HI, mask, fold)
            res["перестановки"].append({"perm": list(pm), "score": sc})
            print(f"        нуль {pm}  пара {sc[0]:.4f}  ранг {sc[1]:.4f}"
                  f"   dранг {sc[1] - res['равные'][1]:+.4f}")
        out["seeds"][str(s)] = res

    if len(seeds) > 0:
        g = [out["seeds"][str(s)]["1/den"][1] - out["seeds"][str(s)]["равные"][1] for s in seeds]
        n = [np.mean([p["score"][1] for p in out["seeds"][str(s)]["перестановки"]])
             - out["seeds"][str(s)]["равные"][1] for s in seeds]
        print(f"\nсводка по {len(seeds)} сидам, макро-ранг относительно равных весов")
        print(f"  1/den          {np.mean(g):+.4f}   знак {sum(x > 0 for x in g)}/{len(g)}")
        print(f"  перестановки   {np.mean(n):+.4f}   знак {sum(x > 0 for x in n)}/{len(n)}")
        print(f"  макро-пол ранга 0.0036 (пункт 165)")
        out["summary"] = {"gain": g, "null": n}

    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    with open(RES + "logs/k80_denweight.json", "w") as fh:
        json.dump(out, fh, ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
