"""Is the Hill slope identifiable per compound, and if so does it carry kinetics or an assay shift.

The algebra, which is the good part of the proposal and is correct. The instrument equation is
`I = E / (1 + 10^{h(pC0 - pi)})` and the screen reports `log2fc = log2(1 - I)`, so with residual
activity `a = 2^log2fc`,

    10^{h(pC0 - pi)} = (E - 1 + a) / (1 - a),   и при E = 1   h = log10(a/(1-a)) / (pC0 - pi)

Three parameters, and two are supplied rather than assumed: `pi` is the label, and `E` is the Emax
column, whose median runs -0.994 to -1.027 across the four enzymes with a narrow spread. **Item 72
read "Emax has no dynamic range" as closing a proposal; it is also the condition that makes this
solvable.** One equation, one unknown, and only on cells where `pi` was measured -- so unlike item
83's concern nothing is extrapolated and no calibration is transported.

Two things are checked here that the proposal did not check, and either could sink it.

**E is not 1 in this repository's own fit.** `verify/g1_calib.py` fits E per enzyme and gets 0.728,
0.621, 0.867 and 0.931. If Emax says compounds reach full inhibition, a fitted E of 0.621 is not
"they only reach 62 per cent" -- it is the two assays disagreeing, with least_squares absorbing the
disagreement into the amplitude. So h is computed both ways and the difference is reported, because
under the fitted E the formula changes shape and `(E - 1 + a)` can go negative, which is itself a
diagnostic: it marks cells the model cannot represent at all.

**The label has its own error and it enters divided by a small number.** The proposal propagated
`log2fc_std_error` only. But `h = u / (pC0 - pi)`, so `dh/dpi = h / (pC0 - pi)`, and the inhibition
table ships `{CYP}_pIC50_direct_inhibition_std`. With |pC0 - pi| filtered at 0.5 and h around 0.5,
that term is about 1.4 times the label's own sd -- the same order as the spread of h being claimed
as signal. Whether h is identifiable turns on this term, not on the screening error.

**Pre-registered, both outcomes, before the run.**

    отношение sd(h) к ПОЛНОЙ ошибке (оба источника) заметно больше 1
        -> поточечный наклон реален, и это второй химический признак после потенции

    отношение около 1 после включения ошибки метки
        -> то, что выглядело кинетикой, было шумом метки, поделённым на (pC0 - pi),
           и предложение закрывается на собственной арифметике

The decisive second stage, run regardless. Fit a **pair** (d, h) per enzyme -- one shift, one slope
-- against all paired cells, and compare the residual to the measurement error. If two numbers per
enzyme suffice, there is no per-compound kinetics and what has been measured is the inter-assay
calibration, which is the previous item's question with a number attached. If the residual stays
far above the measurement error, per-compound kinetics is real.

Reads the inhibition, Emax and screening tables. Writes nothing. ~1 minute.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
PC0 = 4.305
FIT_E = {"CYP1A2": 0.728, "CYP2C9": 0.621, "CYP2D6": 0.867, "CYP3A4": 0.931}
MIN_GAP = 0.5      # |pC0 - pi|, иначе деление на почти ноль


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    em = pd.read_csv(D + "cyp-challenge-TRAIN_Emax.csv").set_index("Molecule_Name")
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")

    print(f"{'фермент':8s} {'Emax медиана':>13s} {'Emax IQR':>10s} {'подогнанное E':>14s}")
    for c in CYPS:
        k = f"{c}_EmaxVsPosCtrl_direct_inhibition"
        v = em[k].dropna()
        print(f"{c:8s} {v.median():13.3f} "
              f"{np.percentile(v,75)-np.percentile(v,25):10.3f} {FIT_E[c]:14.3f}")
    print("\nEmax около -1 значит полное ингибирование при насыщении -> E = 1.")
    print("Подогнанное E ниже единицы тогда не амплитуда, а несогласие двух анализов.\n")

    res = []
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        sub = sc[sc.enzyme == c].set_index("Molecule_Name")
        d = pd.DataFrame({
            "pi": tr[col].to_numpy(float),
            "pi_sd": tr[col + "_std"].to_numpy(float),
            "l2": sub["log2fc_estimate"].reindex(rows.Molecule_Name).to_numpy(float),
            "l2_sd": sub["log2fc_std_error"].reindex(rows.Molecule_Name).to_numpy(float),
        }).dropna()
        a = 2.0 ** d.l2.to_numpy()
        gap = PC0 - d.pi.to_numpy()
        ok = (a > 1e-6) & (a < 1 - 1e-6) & (np.abs(gap) > MIN_GAP)
        a, gap = a[ok], gap[ok]
        l2sd, pisd = d.l2_sd.to_numpy()[ok], d.pi_sd.to_numpy()[ok]

        h1 = np.log10(a / (1 - a)) / gap                       # E = 1
        num = FIT_E[c] - 1 + a
        rep = num > 1e-9                                        # представимо при подогнанном E
        hE = np.full_like(h1, np.nan)
        hE[rep] = np.log10(num[rep] / (1 - a[rep])) / gap[rep]

        # Распространение ошибки: скрининг И метка.
        s_scr = np.abs(np.log(2) / ((1 - a) * np.log(10) * gap)) * l2sd
        s_lab = np.abs(h1 / gap) * pisd
        s_tot = np.sqrt(s_scr ** 2 + s_lab ** 2)

        # Сдвиг метки s: pi -> pi + s, значит gap -> gap - s. Медиана h как функция s
        # разрывна (знаменатель меняет знак при s = gap у части соединений), поэтому
        # решается сеткой, а не делением пополам --- первая версия упиралась в границу.
        u = np.log10(a / (1 - a))
        grid = np.linspace(-2.5, 2.5, 5001)
        with np.errstate(divide="ignore", invalid="ignore"):
            meds = np.array([np.median(u / (gap - sv)) for sv in grid])
        good = np.isfinite(meds)
        if good.any():
            j = int(np.argmin(np.abs(meds[good] - 1.0)))
            shift = float(grid[good][j])
            shift_miss = float(np.abs(meds[good][j] - 1.0))
        else:
            shift, shift_miss = np.nan, np.nan

        res.append(dict(c=c, n=int(ok.sum()), cov=100 * ok.mean(),
                        sd=np.std(h1, ddof=1), med=np.median(h1),
                        iqr=np.percentile(h1, 75) - np.percentile(h1, 25),
                        frac1=100 * (h1 > 1).mean(),
                        e_scr=np.median(s_scr), e_lab=np.median(s_lab),
                        e_tot=np.median(s_tot), shift=shift, miss=shift_miss,
                        l2sd=float(np.median(l2sd)),
                        unrep=100 * (~rep).mean(),
                        med_hE=np.nanmedian(hE), a=a, gap=gap, l2=d.l2.to_numpy()[ok]))

    print("=== 1. наклон при E = 1: разброс против ОБЕИХ ошибок ===")
    print(f"{'фермент':8s} {'n':>5s} {'покрытие':>9s} {'sd(h)':>7s} "
          f"{'ош.скрин':>9s} {'ош.МЕТКИ':>9s} {'полная':>8s} {'отн.скрин':>10s} {'отн.ПОЛНОЕ':>11s}")
    for r in res:
        print(f"{r['c']:8s} {r['n']:5d} {r['cov']:8.1f} % {r['sd']:7.3f} "
              f"{r['e_scr']:9.3f} {r['e_lab']:9.3f} {r['e_tot']:8.3f} "
              f"{r['sd']/r['e_scr']:10.2f} {r['sd']/r['e_tot']:11.2f}")

    print("\n=== 2. что в нём: кинетика или сдвиг между анализами ===")
    print(f"{'фермент':8s} {'медиана h':>10s} {'IQR':>7s} {'доля h>1':>9s} "
          f"{'сдвиг к h=1':>12s} {'промах':>8s} {'непредставимо':>14s} {'медиана h при подогн. E':>24s}")
    for r in res:
        print(f"{r['c']:8s} {r['med']:10.3f} {r['iqr']:7.3f} {r['frac1']:8.1f} % "
              f"{r['shift']:+12.3f} {r['miss']:8.3f} {r['unrep']:12.1f} % {r['med_hE']:24.3f}")

    print("\n=== 3. решающий: хватает ли ПАРЫ (d, h) на фермент ===")
    print(f"{'фермент':8s} {'d':>8s} {'h':>8s} {'ско остатка':>12s} "
          f"{'ош. измерения':>14s} {'отношение':>10s}")
    for r in res:
        a, gap, l2 = r["a"], r["gap"], r["l2"]
        def resid(t):
            dd, hh = t
            I = 1.0 / (1 + 10 ** (hh * (gap + dd)))
            return np.log2(np.clip(1 - I, 1e-6, None)) - l2
        f = least_squares(resid, [0.0, 1.0], bounds=([-3.0, 0.05], [3.0, 6.0]))
        rr = resid(f.x)
        # Остаток и ошибка измерения в одних единицах --- log2fc. Ошибка --- медиана
        # log2fc_std_error по тем же ячейкам.
        me = r["l2sd"]
        print(f"{r['c']:8s} {f.x[0]:+8.3f} {f.x[1]:8.3f} {rr.std():12.3f} "
              f"{me:14.3f} {rr.std() / me:10.2f}")
    print("""
Как читать. Решает столбец «отн.ПОЛНОЕ» в разделе 1, а не «отн.скрин». Отношение,
посчитанное по одной лишь ошибке скрининга, завышено: наклон получается делением на
(pC0 - pi), поэтому ошибка МЕТКИ входит в него умноженной на h/(pC0 - pi), и при пороге
0.5 это множитель около полутора. Если после включения второго источника отношение падает
к единице, то поточечной кинетики нет --- есть шум метки, поделённый на малое число.

Раздел 2 отвечает на вопрос, что внутри. Медианный наклон около 0.2-0.7 физически невозможен
для обратимого конкурентного ингибирования, где он равен единице, поэтому величина вобрала в
себя сдвиг между двумя анализами. Столбец «сдвиг к h=1» --- его размер в единицах pIC50.

Столбец «непредставимо» --- доля ячеек, где при ПОДОГНАННОМ E выражение (E - 1 + a) не
положительно, то есть прибор с такой амплитудой не может дать наблюдённое чтение вовсе.
Большая доля означает, что подогнанное E и колонка Emax описывают разные вещи.

Раздел 3: если пара (d, h) на фермент даёт остаток вровень с ошибкой измерения, поточечной
кинетики нет, зато измерена межанализовая калибровка.""")


if __name__ == "__main__":
    main()
