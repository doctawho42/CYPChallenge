"""How many compounds carry all three observation channels, and whether the channels are distinct.

The question this gates. A mechanistic proposal on the table separates two latent quantities that
the single pIC50 target confounds: affinity, how tightly a compound binds, and turnover, how
readily the enzyme processes it into something reactive. Direct inhibition sees affinity alone;
time-dependent inhibition sees affinity times turnover times reactivity; Emax sees neither, it
sees what fraction of the enzyme can be inhibited at all, which is mechanism rather than potency.
Three observations and two latents plus a reactivity term from SMARTS alerts is an overdetermined
system, so it would be identifiable - IF the three channels exist on the same molecules in
quantity, and IF they are not collinear.

Both conditions are checked here, because the first is the one people check and the second is the
one that decides. A guess was recorded before running this that the triples would number in the
dozens. They number in the thousands, and guessing rather than counting is a one-command mistake.

The second condition is the real constraint. Direct and TDI pIC50 are measured on the same
compound under two incubation protocols, so they are near-copies of each other: whatever the
second channel adds about turnover lives entirely in the residual between them, and the size of
that residual relative to the signal is the budget the whole construction has to work with.

Reads the three training files. Prints only.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]


def main():
    I = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv").set_index("Molecule_Name")
    T = pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv").set_index("Molecule_Name")
    E = pd.read_csv(D + "cyp-challenge-TRAIN_Emax.csv").set_index("Molecule_Name")
    k = sorted(set(I.index) & set(T.index) & set(E.index))
    print(f"молекул в каждом файле: {len(I)}, {len(T)}, {len(E)}; во всех трёх {len(k)}")
    I, T, E = I.loc[k], T.loc[k], E.loc[k]

    print(f"\n{'фермент':8s} {'альфа':>7s} {'TDI':>7s} {'is_TDI':>7s} {'Emax':>7s} | {'все три':>8s}")
    for c in CYPS:
        a = I[f"{c}_pIC50_direct_inhibition"].notna().to_numpy()
        t = T[f"{c}_pIC50_TDI_condition"].notna().to_numpy()
        e = E[f"{c}_EmaxVsPosCtrl_direct_inhibition"].notna().to_numpy()
        col = f"{c}_is_TDI"
        it = T[col].notna().to_numpy() if col in T.columns else np.zeros(len(k), bool)
        print(f"{c:8s} {a.sum():7d} {t.sum():7d} "
              + (f"{it.sum():7d}" if col in T.columns else f"{'нет':>7s}")
              + f" {e.sum():7d} | {(a & t & e).sum():8d}")

    print(f"\n{'фермент':8s} {'n':>5s} | {'r(a,TDI)':>9s} {'r(a,Emax)':>10s} {'r(TDI,EmaxT)':>12s}"
          f" | {'sd остатка':>10s} {'sd альфа':>8s} {'бюджет':>7s}")
    for c in CYPS:
        a = I[f"{c}_pIC50_direct_inhibition"].to_numpy(float)
        t = T[f"{c}_pIC50_TDI_condition"].to_numpy(float)
        e = E[f"{c}_EmaxVsPosCtrl_direct_inhibition"].to_numpy(float)
        et = E[f"{c}_EmaxVsPosCtrl_TDI_condition"].to_numpy(float)
        m = np.isfinite(a) & np.isfinite(t) & np.isfinite(e) & np.isfinite(et)
        a, t, e, et = a[m], t[m], e[m], et[m]
        d = t - a
        print(f"{c:8s} {int(m.sum()):5d} | {spearmanr(a, t).statistic:9.3f} "
              f"{spearmanr(a, e).statistic:10.3f} {spearmanr(t, et).statistic:12.3f} | "
              f"{d.std():10.3f} {a.std():8.3f} {d.std() / a.std():7.2f}")

    emax_range(E, CYPS)
    print("""
Как читать. Первая таблица снимает возражение про размер: каналов не десятки троек, а весь
размеченный набор, 6524 наблюдения соединение-фермент. Строить есть на чём.

Вторая таблица ставит настоящую границу. Прямой и TDI pIC50 --- это одна молекула в двух
протоколах инкубации, и корреляция 0.91-0.99 говорит, что второй канал почти повторяет первый.
Всё, что он добавляет про оборот, сидит в остатке, а бюджет остатка --- от 0.17 сигнала на
CYP1A2 до 0.45 на CYP2D6. Оценивать латент оборота на 1A2 придётся из шестой части дисперсии.

Третья таблица снимает утверждение, записанное здесь же в первой версии: что Emax --- ось
действительно другая, раз |r| около 0.4. Ранговая корреляция не говорит НИЧЕГО о том, есть ли
там что мерить, и диапазон это показывает. Межквартильный размах 0.02-0.11 при медиане около
-1.0, 99% массы ниже -0.73, выше нуля ни одного соединения из 6524. Каждое соединение в наборе
ингибирует фермент практически полностью, частичных ингибиторов нет, и опереть латент механизма
на почти постоянный канал нельзя.

Это опровергает и предсказание, выведенное из односайтового против многосайтового связывания:
Emax должен был оказаться бимодальным на CYP3A4 и унимодальным на CYP2D6. Бимодальности нет
нигде, а разброс у CYP3A4 самый УЗКИЙ из четырёх, 0.020 против 0.106 у CYP2C9.

Знак связи с Emax на CYP2D6 действительно перевёрнут, +0.770 против -0.394 / -0.407 / -0.462,
но после замера диапазона это факт про порядок внутри полосы шириной 0.05.""")


def emax_range(E, cyps):
    """Диапазон Emax. Печатается рядом с корреляциями, потому что без него они обманывают."""
    print(f"\n{'фермент':8s} {'min':>7s} {'1%':>7s} {'50%':>7s} {'99%':>7s} {'max':>7s} "
          f"{'IQR':>6s} {'доля > -0.5':>11s}")
    for c in cyps:
        v = E[f"{c}_EmaxVsPosCtrl_direct_inhibition"].dropna().to_numpy(float)
        q = np.percentile(v, [1, 50, 99])
        print(f"{c:8s} {v.min():7.2f} {q[0]:7.2f} {q[1]:7.2f} {q[2]:7.2f} {v.max():7.2f} "
              f"{np.subtract(*np.percentile(v, [75, 25])):6.3f} {(v > -0.5).mean():11.3%}")


if __name__ == "__main__":
    main()
