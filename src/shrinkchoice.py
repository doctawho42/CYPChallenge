"""At which delta should the affine pair be fitted, given that delta is only bracketed?

Why this is now the main question. The comparison of the neural trunk against the boosting
came out at 0.005 over the plausible range, while the post-processing those predictions go
through is worth 0.051 - ten times more (src/trunkdose.py, block 6). So the lever that matters
before submission is not which model to ship but what to do with its output.

And that lever has a free parameter we cannot observe. The affine pair (off, lambda) is fitted
on our own label marginal, that is at delta = 0, while the test's marginal sits somewhere
around delta = +0.1 to +0.6 (src/reweight.py, and verify/k5_shift.py from the other end). The
pair that is optimal at delta = 0 is not optimal at delta = 0.4, and we have to choose one
before seeing the test.

The object that answers it is a matrix, not a number: fit under an ASSUMED delta, score under
a TRUE delta. Its diagonal is the oracle that knows the shift. Its top row is what
src/submit.py does today. Reading a column tells you what a given truth costs each choice;
reading a row tells you how a given choice holds up across truths.

From that matrix three candidate rules are compared, all of them choosing one pair per enzyme:

  fit at zero      - today's behaviour, optimal if the test looks like the training set;
  minimax          - the pair whose WORST case over the plausible bracket is least bad;
  mean over the    - the pair minimising the average over the bracket, which is the Bayes
  bracket            rule under a flat prior on delta.

Everything is fitted out of fold: within each fold the pair is chosen on the other four under
the assumed tilt and applied to the held-out one, so no compound contributes to choosing the
parameters that are then applied to it. Four split seeds, and a rule is only called better
than another where the sign holds on all four.

The assumption none of this escapes, stated once: tilting reweights our own labels and takes
p(y | yhat) to be the same on the test. That is exactly what recalibration exists to check,
so this bounds the damage from guessing delta wrong - it does not validate the tilt itself.

Reads results/preds/oof.json and oof_seeds.json. Writes nothing.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import json

import numpy as np
import pandas as pd

from cypsplit import butina_folds
from reweight import tilt

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
OFFGRID = np.round(np.arange(-0.2, 2.01, 0.05), 2)
LAMGRID = np.round(np.arange(0.20, 1.01, 0.02), 2)
DELTAS = np.round(np.arange(-1.2, 1.01, 0.1), 1)   # до -1.2: столько выходит на 2D6
BRACKET = (0.1, 0.6)          # правдоподобный диапазон сдвига, src/reweight.py и verify/k5
SEEDS = [0, 1, 2, 3]


def wnum(p, lo, hi, w):
    """Weighted ST-RAE numerator; the denominator does not depend on the prediction."""
    return float((w * (np.maximum(p - hi, 0.0) + np.maximum(lo - p, 0.0))).sum())


def wden(y, lo, hi, w):
    mu = float((w * y).sum() / w.sum())
    return float((w * (np.maximum(mu - hi, 0.0) + np.maximum(lo - mu, 0.0))).sum())


def fit_apply(p, lo, hi, fold, w):
    """Choose (off, lambda) per fold on the training folds under weights w, apply held out.

    Whole grid at once: the objective is the weighted numerator, and it is a broadcast over
    (off, lambda, compound). Same trick and same result as src/trunkdose.py.
    """
    o = np.full_like(p, np.nan)
    A = OFFGRID[:, None] * (1.0 - LAMGRID[None, :])
    B = np.broadcast_to(LAMGRID[None, :], A.shape)
    for f in np.unique(fold):
        te, trn = fold == f, fold != f
        if te.sum() == 0 or trn.sum() < 10:
            continue
        mu = p[trn].mean()
        q = (mu * (1.0 - B) + A)[:, :, None] + B[:, :, None] * p[trn][None, None, :]
        pen = (w[trn] * (np.maximum(q - hi[trn], 0.0)
                         + np.maximum(lo[trn] - q, 0.0))).sum(axis=2)
        i, j = np.unravel_index(pen.argmin(), pen.shape)
        c, L = mu + OFFGRID[i], LAMGRID[j]
        o[te] = c + L * (p[te] - c)
    return o


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    oof = json.load(open(RES + "preds/oof.json"))
    oofs = json.load(open(RES + "preds/oof_seeds.json"))

    # M[seed][assumed][true] = макро ST-RAE
    M, MC = {}, {}
    for seed in SEEDS:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        per = np.zeros((len(DELTAS), len(DELTAS)))
        perc = {c: np.zeros((len(DELTAS), len(DELTAS))) for c in CYPS}
        for c in CYPS:
            col = f"{c}_pIC50_direct_inhibition"
            m = tr[col].notna().to_numpy()
            y = tr.loc[m, col].to_numpy()
            lo = tr.loc[m, col + "_conf_low"].to_numpy()
            hi = tr.loc[m, col + "_conf_high"].to_numpy()
            p = np.asarray(oof[f"FP+DESC+MECH|{c}"] if seed == 0
                           else oofs[f"{seed}|FP+DESC+MECH|{c}"], float)
            f = fold[m]
            W = {d: tilt(y, float(d)) for d in DELTAS}
            for ia, da in enumerate(DELTAS):
                q = fit_apply(p, lo, hi, f, W[da])
                for it, dt in enumerate(DELTAS):
                    v = wnum(q, lo, hi, W[dt]) / wden(y, lo, hi, W[dt])
                    per[ia, it] += v / 4.0
                    perc[c][ia, it] = v
        M[seed] = per
        MC[seed] = perc
    A = np.mean([M[s] for s in SEEDS], axis=0)

    print("=" * 94)
    print("1. Подогнали при одном сдвиге, проверяем при другом (макро ST-RAE, среднее по сидам)")
    print("=" * 94)
    print("   строка --- при каком delta подогнана пара; столбец --- каким delta оказался тест\n")
    print(f"{'подогнано при':>14s} " + " ".join(f"{('тест ' + str(d)):>9s}" for d in DELTAS))
    for ia, da in enumerate(DELTAS):
        mark = "  <- submit.py" if da == 0.0 else ""
        print(f"{da:14.1f} " + " ".join(f"{A[ia, it]:9.4f}" for it in range(len(DELTAS))) + mark)
    print(f"{'оракул':>14s} " + " ".join(f"{A[it, it]:9.4f}" for it in range(len(DELTAS))))

    inb = [i for i, d in enumerate(DELTAS) if BRACKET[0] - 1e-9 <= d <= BRACKET[1] + 1e-9]
    print()
    print("=" * 94)
    print(f"2. Три правила выбора, оценённые по вилке delta от {BRACKET[0]} до {BRACKET[1]}")
    print("=" * 94)
    rules = {
        "подгонять при 0 (сейчас)": int(np.where(DELTAS == 0.0)[0][0]),
        "минимакс по вилке": int(min(range(len(DELTAS)), key=lambda ia: A[ia, inb].max())),
        "среднее по вилке": int(min(range(len(DELTAS)), key=lambda ia: A[ia, inb].mean())),
    }
    print(f"{'правило':26s} {'подгонять при':>14s} {'худший':>8s} {'средний':>8s} "
          f"{'при 0.1':>8s} {'при 0.6':>8s}")
    for name, ia in rules.items():
        print(f"{name:26s} {DELTAS[ia]:14.1f} {A[ia, inb].max():8.4f} {A[ia, inb].mean():8.4f} "
              f"{A[ia, inb[0]]:8.4f} {A[ia, inb[-1]]:8.4f}")
    orc = np.array([A[it, it] for it in inb])
    print(f"{'оракул (недостижим)':26s} {'--':>14s} {orc.max():8.4f} {orc.mean():8.4f} "
          f"{orc[0]:8.4f} {orc[-1]:8.4f}")

    print()
    print("=" * 94)
    print("3. Сколько стоит ошибиться со сдвигом, и держится ли выигрыш по сидам")
    print("=" * 94)
    i0, ib = rules["подгонять при 0 (сейчас)"], rules["среднее по вилке"]
    print(f"Разность «{DELTAS[ib]:.1f} против 0», по сидам (отрицательное = подгонка при "
          f"{DELTAS[ib]:.1f} лучше):\n")
    print(f"{'':10s} " + " ".join(f"{('тест ' + str(d)):>9s}" for d in DELTAS[inb]))
    for s in SEEDS:
        print(f"{('сид ' + str(s)):10s} "
              + " ".join(f"{M[s][ib, it] - M[s][i0, it]:+9.4f}" for it in inb))
    print(f"{'среднее':10s} " + " ".join(f"{A[ib, it] - A[i0, it]:+9.4f}" for it in inb))
    same = [len({int(np.sign(round(M[s][ib, it] - M[s][i0, it], 4))) for s in SEEDS}) == 1
            for it in inb]
    print(f"{'знак 4/4':10s} " + " ".join(f"{str(x):>9s}" for x in same))

    worst0 = A[i0, inb].max()
    worstb = A[ib, inb].max()
    print()
    print("=" * 94)
    print("4. Поферментно: у каждого фермента свой сдвиг, и на 2D6 он отрицательный")
    print("=" * 94)
    print("   Оценка сдвига по verify/k5_shift.py --- сдвиг собственных выходов модели,")
    print("   делённый на затухание аффинной перекалибровки.")
    print()
    print("   ВНИМАНИЕ на чтение последней колонки. Минимум строки при заданном столбце стоит")
    print("   по построению на диагонали, поэтому «лучшее своё» обязано совпасть с оценкой k5:")
    print("   это арифметика, а не независимое подтверждение оценки. Правая пара колонок ---")
    print("   поферментный ОРАКУЛ при условии, что k5 прав, то есть верхняя граница того, что")
    print("   поферментная подгонка могла бы дать, а не то, что она даст.\n")
    K5 = {"CYP1A2": 0.018, "CYP2C9": 0.164, "CYP2D6": -0.256, "CYP3A4": 0.437}
    AC = {c: np.mean([MC[s][c] for s in SEEDS], axis=0) for c in CYPS}
    print(f"{'фермент':8s} {'ядро':>10s} {'единое 0.3':>11s} {'оракул':>8s} "
          f"{'строка оракула':>15s}")
    glob_i = int(min(range(len(DELTAS)), key=lambda ia: A[ia, inb].mean()))
    tot_g = tot_o = 0.0
    for c in CYPS:
        # честная проверка: столбец, отвечающий оценке k5 для этого фермента
        it = int(np.argmin(np.abs(DELTAS - max(K5[c], 0.0))))
        best_own = int(np.argmin(AC[c][:, it]))
        tot_g += AC[c][glob_i, it] / 4.0
        tot_o += AC[c][best_own, it] / 4.0
        print(f"{c:8s} {K5[c]:+10.3f} {AC[c][glob_i, it]:11.4f} {AC[c][best_own, it]:8.4f} "
              f"{DELTAS[best_own]:15.1f}")
    print(f"\n  единое правило (подгонка при {DELTAS[glob_i]:.1f}), если k5 прав:  {tot_g:.4f}")
    print(f"  поферментный оракул, если k5 прав:                {tot_o:.4f}")
    print(f"  верхняя граница выигрыша поферментной подгонки:  {tot_g - tot_o:+.4f}")
    print("""
  Оговорка. Оценки k5 сами имеют разброс, а числа выше посчитаны так, будто он нулевой;
  реальный выигрыш меньше этой границы настолько, насколько k5 неточен.

  Отрицательная оценка на 2D6 --- не ошибка и не противоречие с вилкой +0.1..+0.6. Вилка
  выведена из трёх ферментов, по которым отбирали якоря; 2D6 в отборе не участвовал, и его
  число в reweight.py --- внутренний контроль, а не оценка. Тест обеднён основаниями вдвое,
  а 2D6 --- единственный фермент, у которого основания активнее. Разбор: verify/k7_2d6shift.py.""")

    print()
    print("=" * 94)
    print("5. Поферментная подгонка, оценённая по поферментным же диапазонам истины")
    print("=" * 94)
    print("""   Диапазон --- кластерный бутстрап ядерной оценки (verify/k8_kernel.py). Ядро
   обходит вопрос о коэффициенте передачи целиком: вместо того чтобы делить сдвиг
   предсказаний на неизвестный множитель, оно наклоняет маргиналь метки, пока среднее
   предсказание не сойдётся с наблюдённым на тесте. Инструмент не нужен.
   Правила сравниваются по ХУДШЕМУ случаю внутри диапазона, а не в одной точке, иначе
   выбор правила под оценку и проверка по ней же были бы одним и тем же действием.\n""")
    # ИСПРАВЛЕНО. Прежние диапазоны строились из кластерного бутстрапа сдвига ПРЕДСКАЗАНИЙ
    # и использовались как диапазоны сдвига МЕТОК. Это разные величины: delta связан со
    # сдвигом предсказаний как delta = d_yhat * b / R^2, а множитель b/R^2 равен
    # 3.40 / 2.43 / 4.48 / 1.69. Все четыре диапазона были не в той шкале, и выбранные по
    # ним сдвиги недействительны. Ошибку нашла внешняя проверка, хотя и по неверному следу:
    # там решили, что интервал делили на b, тогда как его не переводили вовсе.
    #
    # Ниже --- оценка через ЯДРО (verify/k8_kernel.py), которая обходится без коэффициента
    # передачи совсем: E_test[yhat] = int E[yhat|y] p_test(y) dy, где ядро E[yhat|y]
    # оценивается изотоникой на OOF, а маргиналь метки наклоняется до совпадения среднего с
    # наблюдённым тестовым. Инструмент не нужен, вопрос "чему равна передача" не ставится.
    # Диапазоны --- её кластерный бутстрап, округлённый наружу до узлов сетки.
    # CYP2D6 --- по СТРАТИФИЦИРОВАННОМУ ядру (verify/k10_strat2d6.py). Объединённое ядро
    # приписывало сдвигу меток и ту часть падения предсказаний, которая идёт от смены
    # состава: на тесте оснований 0.104 против 0.355 в маске 2D6, а модель предсказывает
    # основаниям высоко почти независимо от истинной активности. Ядра страт там расходятся
    # на 0.35-0.58 по всей шкале. После разделения -0.508 вместо -0.917, из них -0.138
    # чистый состав. На трёх остальных ферментах доли различаются на один пункт, и
    # стратификация двигает оценку на 0.003-0.007, то есть ни на что: это контроль
    # конструкции, а не результат.
    # ЧЕСТНЫЕ диапазоны. Прежние были кластерным бутстрапом, в котором ресэмплилась только
    # обучающая сторона, а тестовое среднее держалось фиксированным: интервал ловил
    # неопределённость ядра и выбрасывал неопределённость цели, а цель входит в обращение с
    # множителем 1/наклон, и наклоны здесь 0.22-0.59. После ресэмплинга обеих сторон
    # интервалы шире в 2.3-4.7 раза. Дефект был в verify/k8_kernel.py и унаследован нашим
    # k10_strat2d6.py; исправлены оба.
    TRUTH = {"CYP1A2": (-0.3, 0.4), "CYP2C9": (0.1, 0.7),
             "CYP2D6": (-1.1, 0.1), "CYP3A4": (0.5, 1.0)}
    AC = {c: np.mean([MC[s][c] for s in SEEDS], axis=0) for c in CYPS}
    K5 = {"CYP1A2": 0.038, "CYP2C9": 0.369, "CYP2D6": -0.508, "CYP3A4": 0.742}

    def rng_idx(c):
        a, b = TRUTH[c]
        return [i for i, d in enumerate(DELTAS) if a - 1e-9 <= d <= b + 1e-9]

    i0 = int(np.where(np.isclose(DELTAS, 0.0))[0][0])
    ig = int(np.where(np.isclose(DELTAS, 0.3))[0][0])
    print(f"{'фермент':8s} {'диапазон':>13s} {'k5':>7s} | "
          f"{'при 0':>7s} {'при 0.3':>8s} {'своё лучшее':>12s} {'какое':>7s}")
    tot = {"нуль": 0.0, "глобальное 0.3": 0.0, "поферментное": 0.0}
    pick = {}
    for c in CYPS:
        idx = rng_idx(c)
        worst = {ia: AC[c][ia, idx].max() for ia in range(len(DELTAS))}
        best = min(worst, key=worst.get)
        pick[c] = DELTAS[best]
        tot["нуль"] += worst[i0] / 4.0
        tot["глобальное 0.3"] += worst[ig] / 4.0
        tot["поферментное"] += worst[best] / 4.0
        print(f"{c:8s} [{TRUTH[c][0]:+5.1f},{TRUTH[c][1]:+5.1f}] {K5[c]:+7.3f} | "
              f"{worst[i0]:7.4f} {worst[ig]:8.4f} {worst[best]:12.4f} {DELTAS[best]:+7.1f}")
    print(f"\n{'макро худшего случая':24s} " + "  ".join(
        f"{k} {v:.4f}" for k, v in tot.items()))
    print(f"\nВыигрыш поферментной подгонки над лучшим глобальным правилом: "
          f"{tot['глобальное 0.3'] - tot['поферментное']:+.4f}")
    print(f"Выигрыш над тем, что делается сейчас:                        "
          f"{tot['нуль'] - tot['поферментное']:+.4f}")
    print(f"\nВыбранные сдвиги: " + ", ".join(f"{c[3:]} {pick[c]:+.1f}" for c in CYPS))
    print()
    print("=" * 94)
    print("6. Худший случай или среднее: критерий делает больше работы, чем оценки")
    print("=" * 94)
    print("""   До сих пор правило выбиралось по ХУДШЕМУ случаю внутри диапазона. Это страховка
   от плохого лидерборда. Если же цель --- лучший ожидаемый балл, правило надо выбирать
   по СРЕДНЕМУ по апостериорному распределению delta, а оно у нас есть: 1500 розыгрышей
   на фермент из verify/k10_strat2d6.py. Критерий был унаследован, а не выбран, и два
   критерия расходятся по знаку производной по ширине диапазона: при усреднении широкий
   разброс тянет к общему правилу, при худшем случае --- от него.\n""")
    dpath = _pl.Path(RES + "preds/delta_draws.json")
    if not dpath.exists():
        print("   нет results/preds/delta_draws.json - сначала verify/k10_strat2d6.py")
    else:
        DR = {c: np.asarray(v, float) for c, v in json.load(open(dpath)).items()}

        def mean_over_post(c, ia):
            """Средний ST-RAE по апостериорному delta при подгонке под ia."""
            return float(np.mean(np.interp(DR[c], DELTAS, AC[c][ia, :])))

        print(f"{'фермент':8s} {'P(d>=0)':>8s} | {'ХУДШИЙ СЛУЧАЙ':>22s} | {'СРЕДНЕЕ':>22s}")
        print(f"{'':8s} {'':8s} | {'при 0':>7s} {'своё':>7s} {'какое':>6s} | "
              f"{'при 0':>7s} {'своё':>7s} {'какое':>6s}")
        tot = {"худший 0": 0.0, "худший своё": 0.0, "среднее 0": 0.0, "среднее своё": 0.0}
        picks = {}
        for c in CYPS:
            idx = rng_idx(c)
            w_worst = {ia: AC[c][ia, idx].max() for ia in range(len(DELTAS))}
            bw = min(w_worst, key=w_worst.get)
            m_mean = {ia: mean_over_post(c, ia) for ia in range(len(DELTAS))}
            bm = min(m_mean, key=m_mean.get)
            picks[c] = (DELTAS[bw], DELTAS[bm])
            tot["худший 0"] += w_worst[i0] / 4.0; tot["худший своё"] += w_worst[bw] / 4.0
            tot["среднее 0"] += m_mean[i0] / 4.0; tot["среднее своё"] += m_mean[bm] / 4.0
            print(f"{c:8s} {float((DR[c] >= 0).mean()):8.3f} | "
                  f"{w_worst[i0]:7.4f} {w_worst[bw]:7.4f} {DELTAS[bw]:+6.1f} | "
                  f"{m_mean[i0]:7.4f} {m_mean[bm]:7.4f} {DELTAS[bm]:+6.1f}")
        print(f"\n{'макро':8s} {'':8s} | {tot['худший 0']:7.4f} {tot['худший своё']:7.4f} "
              f"{'':6s} | {tot['среднее 0']:7.4f} {tot['среднее своё']:7.4f}")
        print(f"\nвыигрыш поферментной подгонки: по худшему случаю "
              f"{tot['худший 0'] - tot['худший своё']:+.4f}, "
              f"по среднему {tot['среднее 0'] - tot['среднее своё']:+.4f}")

        print()
        print("=" * 94)
        print("7. Симметричная проверка: не хуже ли выбранное правило, чем бездействие")
        print("=" * 94)
        print("""   Единое delta = 0.3 забраковано тем, что на 2D6 оно ХУЖЕ бездействия. Тем же
   стандартом надо проверить и поферментное правило: на ферментах, чей диапазон
   оседлал ноль, выбранный сдвиг может оказаться таким же дефектом, переехавшим в
   другую ячейку.\n""")
        print(f"{'фермент':8s} | {'ХУДШИЙ СЛУЧАЙ':>28s} | {'СРЕДНЕЕ':>28s}")
        print(f"{'':8s} | {'выбрано':>8s} {'P(хуже)':>8s} {'вердикт':>10s} | "
              f"{'выбрано':>8s} {'P(хуже)':>8s} {'вердикт':>10s}")
        for c in CYPS:
            v_zero = np.interp(DR[c], DELTAS, AC[c][i0, :])
            cells = []
            for pick in picks[c]:
                ib_ = int(np.where(np.isclose(DELTAS, pick))[0][0])
                v = np.interp(DR[c], DELTAS, AC[c][ib_, :])
                pw = float((v > v_zero).mean())
                verdict = "чисто" if pw < 0.25 else ("СПОРНО" if pw < 0.5 else "ХУЖЕ")
                cells.append(f"{pick:+8.1f} {pw:8.3f} {verdict:>10s}")
            print(f"{c:8s} | {cells[0]:>28s} | {cells[1]:>28s}")
        print("""
   Читать так: если правило проигрывает бездействию заметно чаще, чем изредка, оно
   несёт ровно тот дефект, за который отвергнуто единое правило.""")

    print("""
Обратите внимание, что выбранные сдвиги не совпадают с оценками k5 и не обязаны: правило
выбирается по худшему случаю на всём диапазоне, а не по одной точке внутри него. Совпадение
означало бы, что диапазон вырожден, то есть что мы знаем сдвиг точно.""")

    print(f"""
Как это читать. Верхняя строка первой таблицы --- то, что происходит сегодня: пара подогнана
при нулевом сдвиге и встречает тест, сдвинутый на неизвестную величину. Её худший случай по
вилке {worst0:.4f}, у правила «подгонять под середину вилки» --- {worstb:.4f}. Разница между
ними и есть цена того, что мы не знаем delta, --- и её стоит сравнить с 0.005, во столько
обходится весь выбор между бустингом и сетью.

Чего таблица не говорит. Она перевзвешивает наши собственные метки и считает, что связь
предсказания с истиной на тесте та же. Если это не так, ни одна строка не верна, и проверить
это можно только на самом тесте.""")


if __name__ == "__main__":
    main()
