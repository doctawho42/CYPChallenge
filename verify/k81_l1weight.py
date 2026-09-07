"""L1 вместе с весами на пулированном члене. Предрегистрация --- пункт 259.

Пункт 258 нашёл рассогласование: метрика штрафует ошибку кусочно ЛИНЕЙНО, а gbm_reg() учится
под квадратом, где вес c эквивалентен масштабированию ошибки строки на sqrt(c). Значит веса
1/den всё это время взвешивали КВАДРАТИЧНУЮ ошибку и метрику не воспроизводили. Честная пара ---
loss="absolute_error" ВМЕСТЕ с sample_weight=1/den, и на пулированном члене она не мерялась
никогда (поферментные L1 мерялись многократно: пункты 74, 75, 77, 78, 80, 146, 148, 226).

Квадрат 2x2 внутри ОДНОГО скрипта, чтобы сравнение было внутренним:

    потеря  in {squared_error (пины gbm_reg), absolute_error (пины DZ_KW)}
      x  вес  in {явные единицы, 1/den_e нормированные на среднее 1}

ВСЕ ЧЕТЫРЕ ЯЧЕЙКИ ПЕРЕДАЮТ ЯВНЫЙ sample_weight, единичный в том числе. Это не педантизм.
В pinned scikit-learn 1.3.2 AbsoluteError.fit_intercept_only ветвится на `sample_weight is None`:
без весов зовётся np.median, интерполирующая два центральных остатка, с весами ---
_weighted_percentile, возвращающая нижний. Единичный вес под L1 НЕ воспроизводит невзвешенный
вызов, и если бы равновесное плечо звало невзвешенный путь (как в k80), измеренный эффект весов
содержал бы переключение оценщика листа. Режим --switch меряет этот разрыв отдельно.

Десять сидов, а не четыре: при sd парной разности 0.0053 четыре сида дают мощность 0.19 против
эффекта 0.0036. Три перестановочных плеча из k80 выброшены --- в пункте 258 именно это условие
единственным прошло полный провал, --- и счёт куплен сидами при том же времени.

Основная статистика --- ПОСИДОВОЕ ПАРНОЕ ВЗАИМОДЕЙСТВИЕ I = (L1,w - L1,1) - (L2,w - L2,1), и его
разброс считается по n значениям I, а не переносится с двух разностей по отдельности.

МАКРО ПО ТРЁМ ФЕРМЕНТАМ. SOLO = {"CYP3A4": ("GP",)} означает, что пулированный член в подаваемое
плечо CYP3A4 не входит вовсе (пункт 218), поэтому его колонка CYP3A4 для решения об этом члене ---
мёртвый груз. Макро по четырём печатается рядом.

Условия приёмки записаны в пункте 259 до прогона и здесь не повторяются.

Читает data/feats.npz, data/rows.csv и обучающую таблицу. Пишет results/logs/k81_l1weight.json
после КАЖДОГО сида. Полный прогон ~2 часа, строго одним процессом.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse, json, time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, t as tdist
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply
from submit import CYPS, gbm_reg, pooled_design, DZ_KW, SOLO, _oof_one

SHIP = [e for e, c in enumerate(CYPS) if SOLO.get(c) is None or "пул" in SOLO[c]]


def model(loss):
    """Одна из ДВУХ пинованных конфигураций репозитория, без изобретения третьей."""
    return gbm_reg() if loss == "squared_error" else HistGradientBoostingRegressor(**DZ_KW)


def den(yv, lov, hiv):
    """Знаменатель ST-RAE: константный предсказатель в mean(y) через тот же мягкий порог."""
    m = yv.mean()
    return float(np.maximum(0.0, np.maximum(lov - m, m - hiv)).sum())


def wvec(y, LO, HI, mask, fold, f):
    """1/den по ферментам на ОБУЧАЮЩИХ строках фолда f, нормировка на среднее 1."""
    d = np.array([den(y[mask[:, e] & (fold != f), e], LO[mask[:, e] & (fold != f), e],
                      HI[mask[:, e] & (fold != f), e]) for e in range(len(CYPS))])
    v = 1.0 / d
    return v / v.mean()


def pooled(X, y, LO, HI, mask, fold, loss, weighted, raw=False):
    """Пулированный член вне фолда. raw=True --- не передавать sample_weight вовсе."""
    P = [np.zeros(int(mask[:, e].sum())) for e in range(len(CYPS))]
    for f in range(5):
        v = wvec(y, LO, HI, mask, fold, f) if weighted else np.ones(len(CYPS))
        Xs, ys, ws = [], [], []
        for e in range(len(CYPS)):
            sel = mask[:, e] & (fold != f)
            if sel.any():
                Xs.append(pooled_design(X[sel], e))
                ys.append(y[sel, e])
                ws.append(np.full(int(sel.sum()), float(v[e])))
        Xtr, ytr = np.vstack(Xs), np.concatenate(ys)
        m_ = model(loss)
        if raw:
            m_.fit(Xtr, ytr)
        else:
            w = np.concatenate(ws)
            m_.fit(Xtr, ytr, sample_weight=w / w.mean())
        for e in range(len(CYPS)):
            mk = mask[:, e]
            b = fold[mk] == f
            if b.sum() == 0:
                continue
            P[e][b] = m_.predict(pooled_design(X[mk][b], e))
    return P


def score(P, y, LO, HI, mask, fold):
    """(макро3 ранг, макро4 ранг, макро4 пара, ранг по ферментам, пара по ферментам)."""
    pr, rk = [], []
    for e in range(len(CYPS)):
        m = mask[:, e]
        yy, lo, hi, p = y[m, e], LO[m, e], HI[m, e], P[e]
        u = np.ones(len(yy)) / len(yy)
        pr.append(float(strae(yy, fit_apply(p, lo, hi, fold[m], u), y_true_upper=hi, y_true_lower=lo)))
        rk.append(float(spearmanr(yy, p).statistic))
    return (float(np.mean([rk[e] for e in SHIP])), float(np.mean(rk)),
            float(np.mean(pr)), rk, pr)


def lower95(x):
    """Односторонняя нижняя 95 %-граница среднего по t с n-1 степенями."""
    x = np.asarray(x, float)
    n = len(x)
    return float(x.mean() - tdist.ppf(0.95, n - 1) * x.std(ddof=1) / np.sqrt(n))


def load():
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    return X, rows, y, LO, HI, ~np.isnan(y)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3,4,5,6,7,8,9")
    ap.add_argument("--switch", action="store_true",
                    help="только диагностика переключения оценщика листа под L1, ~6 минут")
    a = ap.parse_args()
    X, rows, y, LO, HI, mask = load()
    P = print

    if a.switch:
        fold, _ = butina_folds(list(rows.SMILES))
        P("Разрыв между невзвешенным вызовом и ЕДИНИЧНЫМ sample_weight, сид 0\n", flush=True)
        P(f"  {'потеря':16s} {'макро3':>8s} {'макро4':>8s} {'max|d|':>8s} {'ро':>8s} {'ср. сдвиг':>10s}",
          flush=True)
        for loss in ("squared_error", "absolute_error"):
            A = pooled(X, y, LO, HI, mask, fold, loss, False, raw=True)
            B = pooled(X, y, LO, HI, mask, fold, loss, False, raw=False)
            sa, sb = score(A, y, LO, HI, mask, fold), score(B, y, LO, HI, mask, fold)
            d = np.concatenate([B[e] - A[e] for e in range(len(CYPS))])
            ro = np.mean([spearmanr(A[e], B[e]).statistic for e in range(len(CYPS))])
            P(f"  {loss:16s} {sb[0] - sa[0]:+8.4f} {sb[1] - sa[1]:+8.4f} "
              f"{np.abs(d).max():8.4f} {ro:8.4f} {d.mean():+10.4f}", flush=True)
        P("\n(разность = единичный вес минус невзвешенный; под квадратом обязана быть нулём)",
          flush=True)
        return

    seeds = [int(s) for s in a.seeds.split(",") if s.strip()]
    CELLS = [("L2", "squared_error", False), ("L2w", "squared_error", True),
             ("L1", "absolute_error", False), ("L1w", "absolute_error", True)]
    out = {"ship": [CYPS[e] for e in SHIP], "seeds": {}}
    P(f"Пулированный член ходит в подачу на: {', '.join(out['ship'])}   "
      f"(SOLO выводит CYP3A4 на GP один)\n", flush=True)
    P(f"  {'сид':>3s} {'ячейка':>5s} {'макро3':>8s} {'макро4':>8s} {'пара4':>8s} {'c':>6s}", flush=True)
    for s in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        res = {}
        for name, loss, wt in CELLS:
            t0 = time.time()
            sc = score(pooled(X, y, LO, HI, mask, fold, loss, wt), y, LO, HI, mask, fold)
            res[name] = sc
            P(f"  {s:3d} {name:>5s} {sc[0]:8.4f} {sc[1]:8.4f} {sc[2]:8.4f} {time.time() - t0:6.0f}",
              flush=True)
        if s == 0 and abs(res["L2"][1] - 0.5792) > 5e-4:
            P(f"  ВНИМАНИЕ: якорь не сошёлся, L2 x единицы = {res['L2'][1]:.4f} против 0.5792 "
              f"(пункт 146). Разбиение или пины сдвинулись.", flush=True)
        dl1 = res["L1w"][0] - res["L1"][0]
        dl2 = res["L2w"][0] - res["L2"][0]
        P(f"      dL1 {dl1:+.4f}   dL2 {dl2:+.4f}   I {dl1 - dl2:+.4f}", flush=True)
        out["seeds"][str(s)] = res
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(RES + "logs/k81_l1weight.json", "w"), ensure_ascii=False, indent=1)

    ks = [str(s) for s in seeds]
    I = np.array([out["seeds"][k]["L1w"][0] - out["seeds"][k]["L1"][0]
                  - (out["seeds"][k]["L2w"][0] - out["seeds"][k]["L2"][0]) for k in ks])
    D1 = np.array([out["seeds"][k]["L1w"][0] - out["seeds"][k]["L1"][0] for k in ks])
    D2 = np.array([out["seeds"][k]["L2w"][0] - out["seeds"][k]["L2"][0] for k in ks])
    DL = np.array([out["seeds"][k]["L1"][0] - out["seeds"][k]["L2"][0] for k in ks])
    P(f"\nмакро по трём ферментам, {len(ks)} сидов\n", flush=True)
    P(f"  {'величина':34s} {'среднее':>9s} {'sd':>8s} {'нижн.95%':>9s} {'знак':>6s} {'пол':>7s}",
      flush=True)
    for nm, v, fl in (("веса под L1   (L1w-L1)", D1, 0.0052),
                      ("веса под L2   (L2w-L2)", D2, 0.0052),
                      ("L1 против L2 при единицах", DL, 0.0052),
                      ("ВЗАИМОДЕЙСТВИЕ I", I, 0.0074)):
        P(f"  {nm:34s} {v.mean():+9.4f} {v.std(ddof=1):8.4f} {lower95(v):+9.4f} "
          f"{int((v > 0).sum()):3d}/{len(v):<2d} {fl:7.4f}", flush=True)
    P(f"\n  условия пункта 259: (1) нижн.95% I > 0 --- {'ДА' if lower95(I) > 0 else 'НЕТ'};"
      f"  (2) среднее I > 0.0074 --- {'ДА' if I.mean() > 0.0074 else 'НЕТ'};"
      f"  (3) среднее dL1 > 0.0052 --- {'ДА' if D1.mean() > 0.0052 else 'НЕТ'}", flush=True)
    P(f"\n  по ферментам, среднее движение под весами:", flush=True)
    P(f"  {'фермент':9s} {'подаётся':>9s} {'dL2':>9s} {'dL1':>9s} {'знак dL1':>9s}", flush=True)
    for e, c in enumerate(CYPS):
        a1 = np.array([out["seeds"][k]["L1w"][3][e] - out["seeds"][k]["L1"][3][e] for k in ks])
        a2 = np.array([out["seeds"][k]["L2w"][3][e] - out["seeds"][k]["L2"][3][e] for k in ks])
        P(f"  {c:9s} {('да' if e in SHIP else 'нет'):>9s} {a2.mean():+9.4f} {a1.mean():+9.4f} "
          f"{int((a1 > 0).sum()):5d}/{len(a1):<2d}", flush=True)
    out["summary"] = {"I": I.tolist(), "dL1": D1.tolist(), "dL2": D2.tolist(), "lossOnly": DL.tolist()}
    json.dump(out, open(RES + "logs/k81_l1weight.json", "w"), ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
