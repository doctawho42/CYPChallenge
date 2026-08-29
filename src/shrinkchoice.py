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
DELTAS = np.round(np.arange(0.0, 0.81, 0.1), 1)
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
    print(f"{'фермент':8s} {'оценка k5':>10s} {'единое 0.3':>11s} {'оракул':>8s} "
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
