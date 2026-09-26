"""The dead zone was measured on four ensemble members; the submission has five.

The gap this closes. `src/abldzens.py` reads `oof_pool_all`, `oof_gp` and `oof_weak`, so its
members are the per-enzyme boosting, the pooled boosting, the Gaussian process and the ridge --
**four**. Item 121 wired a fifth, the joint-likelihood trunk, behind `--mode ансамбль5`, and item
120 measured what it is worth:

    состав                пара      ранг
    четыре              0.6819    0.6009
    пять (lam = 3)      0.6758    0.6063

So the best number in this file, item 162's 0.6196, belongs to a configuration that **is not the
one that gets submitted**, and the dead zone has never been measured on the submitted one. That is
the kind of mismatch this file exists to catch, and it was caught by reading rather than by
running.

Refitting the trunk against the reprojected band would mean re-running `src/trunk.py` on torch,
which is an hour and a different environment. It is not needed to answer the question asked. Both
`src/submit.py` (`_combine`) and `src/abldzens.py` (the `arms` selectors, averaged in `main`)
combine members by an **unweighted mean**, so

    пять = (4 * четыре + ствол) / 5

exactly, and the four-member averages are already saved. This script therefore composes rather than
computes, and runs in minutes:

    пять, обычный                 must reproduce item 120's 0.6063, or the composition is wrong
    пять, мёртвая зона в четырёх  the untested configuration: the four refitted members plus the
                                  trunk as it stands
    четыре, мёртвая зона          item 162's 0.6196, for reference

What it cannot answer, stated so it is not read as answered: the trunk itself is not reprojected,
so this is a lower bound on what a fully dead-zone five-member ensemble would score. Every other
member gained from the pass -- +0.0435, +0.0462, +0.0188, +0.0174 -- so the trunk plausibly would
too, and the number here should be read as "at least".

Reads results/preds/oof_dzens.json, oof_dzens123.json and trunk_twohead.json. Writes nothing.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
LAM = "3.0"        # принятая настройка ствола, пункт 120
TRUNK_DEFAULT = "trunk_twohead.json"


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    import argparse as _ap
    _p = _ap.ArgumentParser()
    _p.add_argument("--no-clip", dest="clip", action="store_false",
                    help="НЕ обрезать ствол диапазоном меток. Воспроизводит поведение до "
                         "пункта 215, на котором посчитаны пункты 120, 121 и 164.")
    _p.add_argument("--trunk", default=TRUNK_DEFAULT,
                    help="какой файл ствола подставлять пятым членом")
    _p.add_argument("--lam", default=LAM)
    _a = _p.parse_args()
    lam = _a.lam
    T = json.load(open(RES + "preds/" + _a.trunk))["preds"]
    # ВНИМАНИЕ: имя `tag` в этом файле уже занято циклом по именам метрик
    # (`for tag in ("пара", "rho")`), и переиспользование молча роняло все сиды
    # кроме первого --- таблица при этом печаталась как ни в чём не бывало.
    mode_tag = _a.trunk.replace("trunk_", "").replace(".json", "")
    print(f"ствол: {_a.trunk}, режим {mode_tag}, lambda {lam}")

    P = {}
    for fn in ("oof_dzens.json", "oof_dzens123.json"):
        path = RES + "preds/" + fn
        if os.path.exists(path):
            P.update(json.load(open(path))["preds"])
    seeds = sorted({int(k.split("|")[0]) for k in P})
    print(f"сиды с сохранённым ансамблем: {seeds}\n")

    rec = []
    for seed in seeds:
        tk = f"{mode_tag}|{seed}|{lam}"
        if tk not in T:
            print(f"(нет ствола для сида {seed}, пропускаю)")
            continue
        trunk_all = np.asarray(T[tk], float)          # (4905, 4)
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for arm4, name5 in (("базовый ансамбль", "пять, обычный"),
                            ("мёртвая зона везде", "пять, МЗ в четырёх")):
            for label, use_trunk in ((name5, True), (f"четыре, {arm4}", False)):
                r = {"seed": seed, "состав": label}
                ok = True
                for e, c in enumerate(CYPS):
                    key = f"{seed}|{arm4}|{c}"
                    if key not in P:
                        ok = False
                        break
                    col = f"{c}_pIC50_direct_inhibition"
                    m = tr[col].notna().to_numpy()
                    y = tr.loc[m, col].to_numpy()
                    lo = tr.loc[m, col + "_conf_low"].to_numpy()
                    hi = tr.loc[m, col + "_conf_high"].to_numpy()
                    a4 = np.asarray(P[key], float)
                    if len(a4) != int(m.sum()):
                        ok = False
                        break
                    # Ансамбль --- невзвешенное среднее членов (submit.py, _combine),
                    # поэтому пятичленный собирается из четырёхчленного точно.
                    #
                    # Обрезка ствола диапазоном меток фермента плюс-минус две единицы ---
                    # ровно то, что делает submit._trunk_clip. Без неё этот файл описывал
                    # конфигурацию, СОСЕДНЮЮ с подаваемой, а не подаваемую: пункт 205
                    # намерил, что проход выброс не убирает, а переселяет (на CYP2D6
                    # минимум -360.26 стал максимумом +76.82), так что в одном члене из
                    # пяти усреднялось число порядка семидесяти.
                    tk_e = trunk_all[m, e]
                    if _a.clip:
                        tk_e = np.clip(tk_e, np.nanmin(y) - 2.0, np.nanmax(y) + 2.0)
                    p = (4.0 * a4 + tk_e) / 5.0 if use_trunk else a4
                    q = fit_apply(p, lo, hi, fold[m], np.ones(len(y)) / len(y))
                    r[f"{c} пара"] = float(strae(y, q, y_true_upper=hi, y_true_lower=lo))
                    r[f"{c} rho"] = float(spearmanr(y, p).statistic)
                if not ok:
                    continue
                for tag in ("пара", "rho"):
                    r[f"MACRO {tag}"] = float(np.mean([r[f"{c} {tag}"] for c in CYPS]))
                rec.append(r)

    df = pd.DataFrame(rec)
    g = df.groupby("состав", sort=False)[["MACRO пара", "MACRO rho"]
                                         + [f"{c} rho" for c in CYPS]].mean()
    print(g.round(4).to_string())
    print(f"\nсидов усреднено: {df.groupby('состав').size().to_dict()}")
    print("""
Как читать. «пять, обычный» обязан лечь около 0.6758 / 0.6063 из пункта 120 --- это
проверка самой сборки, а не результат. Если он не ложится, среднее не невзвешенное или
предсказания не те, и остальное читать нельзя.

Решает «пять, МЗ в четырёх» против «пять, обычный»: столько мёртвая зона стоит В ПОДАВАЕМОЙ
конфигурации, а не в измерявшейся. И против «четыре, мёртвая зона везде» --- добавляет ли
ствол что-нибудь поверх мёртвой зоны или она уже забрала то, чем он был полезен.

Ствол здесь НЕ перепроецирован, поэтому число --- нижняя оценка. Остальные четыре члена
выиграли от прохода от +0.0174 до +0.0462, так что и он, вероятно, выиграл бы.""")


if __name__ == "__main__":
    main()
