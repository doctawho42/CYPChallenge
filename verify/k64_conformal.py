"""The uncertainty the model itself produces: does it predict error, and does it cover?

Section 9 of the document lists three parts of an uncertainty layer. Two are built and negative
-- the heteroscedastic decision layer failed its own pre-registered threshold (item 93), and the
split normal under it failed separately (item 200). The third, conformal calibration, is the one
that has never been run, and its premise as written -- "on a split reproducing the design of the
test" -- was refuted separately: no such split exists in our data (items 123, 129, 206).

But plain split conformal does not need that premise. It needs exchangeability between the
calibration rows and the rows it is applied to, which the Butina folds supply as well as they
supply anything else in this file. So the part can be run as stated minus the impossible clause,
and this file runs it.

**Item 220 makes this worth doing now rather than earlier.** Until it, the aleatoric floor was
unmeasured: the band shipped with the labels is a function of the label (item 114), so it could
not say how much of our error is irreducible. Item 220 measured the assay's own noise from 912
control wells at 0.163 in log2fc, which propagates to roughly 0.03 to 0.05 of pIC50 through a
twelve-point fit. Against a model RMSE near 0.7 that settles the decomposition by arithmetic:
**almost all of our uncertainty is epistemic**, and that agrees with item 194, which found the
ensemble limited by data rather than by model diversity. Two different measurements, one
conclusion.

Three questions, in order of how much they decide.

1. **Does the Gaussian process's own predictive variance predict error?** It is the only member
   producing one natively. If |error| does not rise with the predicted spread, no normalisation
   built on it can help and conditional coverage is out of reach.
2. **Does plain split conformal cover?** Marginal coverage at the nominal level is the weakest
   possible check and it should pass by construction; it is run because a failure would mean the
   folds are not exchangeable, which would matter far beyond this file.
3. **Does coverage hold CONDITIONALLY, by activity zone?** This is the question section 9 raises
   and the one that matters: inactive compounds carry wide bands by construction, so a model that
   misses systematically on the actives can still show correct coverage on average. If conditional
   coverage fails, a marginal guarantee is worthless for the use anyone would put it to.

Normalised conformal -- dividing residuals by the GP's predicted spread before taking the
quantile -- is run alongside, since that is the standard route to conditional coverage and
question 1 decides whether it can work at all.

Reads the member cache written by verify/k58_dzsubmit.py --cache and results/preds/oof_gp.json.
Minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from cypsplit import butina_folds

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
ALPHA = 0.10                      # номинальное покрытие 90 %
FLOOR = 0.04                      # алеаторный пол в pIC50, пункт 220


def conformal(res, score, fold, alpha=ALPHA):
    """Расщеплённый конформный интервал: квантиль берётся на ОСТАЛЬНЫХ фолдах.

    `score` --- нормировка невязки (единицы, если обычный вариант). Возвращает полуширину
    интервала для каждой строки, посчитанную без единого взгляда на её собственный фолд.
    """
    half = np.zeros(len(res))
    for f in np.unique(fold):
        te, trn = fold == f, fold != f
        if te.sum() == 0 or trn.sum() < 20:
            continue
        s = np.abs(res[trn]) / np.maximum(score[trn], 1e-9)
        n = trn.sum()
        # конечновыборочная поправка расщеплённого конформала
        q = np.quantile(s, min(1.0, np.ceil((n + 1) * (1 - alpha)) / n))
        half[te] = q * score[te]
    return half


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    C = json.load(open(a.cache))
    mem = {k: [np.asarray(v, float) for v in P] for k, P in C["dzp"]}
    GV = json.load(open(RES + "preds/oof_gp.json"))["preds"]
    fold_all, _ = butina_folds(list(rows.SMILES), seed=a.seed)

    print(f"сид {a.seed}, номинальное покрытие {1-ALPHA:.0%}, "
          f"алеаторный пол {FLOOR} pIC50 (пункт 220)\n", flush=True)

    T = []
    for e, c in enumerate(CYPS):
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy()
        fold = fold_all[m]
        M = np.stack([P[e] for P in mem.values()])       # (члены, n)
        p = M.mean(0)
        res = y - p
        # Разброс ансамбля по соединению. Пункт 93 пробовал его как ОБУСЛОВЛИВАЮЩУЮ величину
        # для поправки предсказаний и получил +0.0002; здесь употребление другое --- размер
        # интервала, а не сдвиг точки, --- поэтому предусловие меряется заново ниже.
        spr = M.std(0, ddof=1)
        sd = np.sqrt(np.maximum(np.asarray(GV[f"{a.seed}|GPvar|{c}"], float), 1e-12))

        r = {"фермент": c, "n": len(y), "RMSE": float(np.sqrt((res ** 2).mean())),
             "sd GP медиана": float(np.median(sd)),
             "rho(|ошибка|, sd)": float(spearmanr(np.abs(res), sd).statistic),
             "разброс медиана": float(np.median(spr)),
             "rho(|ошибка|, разброс)": float(spearmanr(np.abs(res), spr).statistic)}

        # Нормировка по ПРЕДСКАЗАНИЮ. Пункт 114: ширина полосы есть детерминированная функция
        # МЕТКИ при R^2 0.93--0.98. Метки на тесте нет, но есть предсказание, и оно с меткой
        # связано. Значит подогнанная на обучающих фолдах зависимость «ожидаемая |невязка| от
        # предсказания» --- законный нормировщик: она не смотрит ни на одну отложенную строку.
        byp = np.ones(len(y))
        for f in np.unique(fold):
            te, trn = fold == f, fold != f
            if te.sum() == 0 or trn.sum() < 50:
                continue
            k = np.polyfit(p[trn], np.abs(res[trn]), 2)
            byp[te] = np.maximum(np.polyval(k, p[te]), 0.05)

        # Комбинация: разброс ансамбля, поднятый до масштаба ожидаемой невязки. Обе части
        # подгоняются только на обучающих фолдах.
        comb = np.ones(len(y))
        for f in np.unique(fold):
            te, trn = fold == f, fold != f
            if te.sum() == 0 or trn.sum() < 50:
                continue
            X2 = np.column_stack([p[trn], p[trn] ** 2, spr[trn], np.ones(trn.sum())])
            k, *_ = np.linalg.lstsq(X2, np.abs(res[trn]), rcond=None)
            Z = np.column_stack([p[te], p[te] ** 2, spr[te], np.ones(int(te.sum()))])
            comb[te] = np.maximum(Z @ k, 0.05)

        for nm, sc in (("обычный", np.ones(len(y))), ("нормированный", sd),
                       ("по предсказанию", byp), ("по разбросу", np.maximum(spr, 0.05)),
                       ("предсказание+разброс", comb)):
            half = conformal(res, sc, fold)
            cov = float((np.abs(res) <= half).mean())
            r[f"{nm}: покрытие"] = cov
            r[f"{nm}: мед. полуширина"] = float(np.median(half))
            # условное покрытие по третям активности
            q = np.quantile(y, [1/3, 2/3])
            zone = np.digitize(y, q)
            r[f"{nm}: покр. по зонам"] = " / ".join(
                f"{float((np.abs(res)[zone == z] <= half[zone == z]).mean()):.3f}" for z in (0, 1, 2))
        T.append(r)

    df = pd.DataFrame(T).set_index("фермент")
    print("1. ПРЕДСКАЗЫВАЕТ ЛИ ДИСПЕРСИЯ GP ОШИБКУ")
    print(df[["n", "RMSE", "sd GP медиана", "rho(|ошибка|, sd)",
              "разброс медиана", "rho(|ошибка|, разброс)"]].round(4).to_string())
    print("\n2. ПОКРЫТИЕ, маргинальное")
    print(df[[f"{n}: покрытие" for n in ("обычный", "нормированный", "по предсказанию", "по разбросу", "предсказание+разброс")]
             + [f"{n}: мед. полуширина" for n in ("обычный", "по предсказанию")]].round(4).to_string())
    print("\n3. ПОКРЫТИЕ ПО ТРЕТЯМ АКТИВНОСТИ (слабые / средние / сильные)")
    print(df[[f"{n}: покр. по зонам" for n in
              ("обычный", "нормированный", "по предсказанию")]].to_string())
    VAR = ("обычный", "нормированный", "по предсказанию", "по разбросу",
           "предсказание+разброс")

    def zones(col):
        return np.array([[float(x) for x in v.split(" / ")] for v in df[col]])

    print("\n   СВОДКА: разброс покрытия между зонами (меньше --- лучше)")
    print(f"   {'нормировщик':22s} {'разброс':>8s} {'закрыто':>8s} {'полуширина':>11s}"
          f"   слаб / сред / сильн")
    base = None
    for n in VAR:
        Z = zones(f"{n}: покр. по зонам")
        sp = float(np.mean(Z.max(1) - Z.min(1)))
        if base is None:
            base = sp
        hw = float(df[f"{n}: мед. полуширина"].mean())
        z = Z.mean(0)
        print(f"   {n:22s} {sp:8.3f} {100*(base-sp)/base:7.1f}% {hw:11.3f}   "
              + " / ".join(f"{x:.3f}" for x in z))

    rm = float(df.RMSE.mean()); hw = float(df["обычный: мед. полуширина"].mean())
    print(f"""
4. РАЗЛОЖЕНИЕ, теперь считаемое
   средняя RMSE модели                 {rm:.3f} pIC50
   алеаторный пол (пункт 220)          {FLOOR:.3f} pIC50
   доля дисперсии, объяснимая шумом    {(FLOOR/rm)**2:.1%}
   значит эпистемической                {1-(FLOOR/rm)**2:.1%}

Как читать. Пункт 3 решает. Маргинальное покрытие обязано выйти около номинала по построению
расщеплённого конформала --- если не вышло, невзаимозаменяемы фолды, и это важнее всего
остального в файле. А вот покрытие ПО ЗОНАМ ничем не гарантировано: у слабых соединений полосы
широки по построению, и систематический промах на сильных маргинальная цифра замаскирует.

Нормировка на дисперсию GP имеет смысл ровно настолько, насколько пункт 1 даёт положительную
связь. Если rho около нуля, нормированный вариант --- это деление на шум, и он должен покрытие
по зонам ухудшить, а не улучшить.

Пункт 4 --- арифметика, а не измерение: он делит измеренный в пункте 220 пол на нашу ошибку.
Вывод «неопределённость почти вся эпистемическая» согласуется с пунктом 194, где ансамбль
оказался ограничен данными, а не разнообразием моделей. Два разных измерения, один вывод.""")


if __name__ == "__main__":
    main()
