"""Is the instrument calibration a property of the instrument, or of the population it was fitted on.

The concern, and it is well aimed. `verify/g1_calib.py` fits one pair (E, h) per enzyme on the mask
`m[col].notna() & m[c].notna()` -- compounds that have **both** a curve and a screening reading.
Item 129 established that CYP3A4's training set is two populations: 1805 compounds from the
diversity library, which were screened, and a 530-compound analog campaign, which was not. So on
CYP3A4 the calibration is fitted on 1805 rows and applied to 2335. The other three enzymes have
**zero** labelled compounds without a screening reading, so this is a CYP3A4 question only.

    группа        n   медиана y     IQR   медиана полосы   доля y < pC0
    в скрине   1805       4.203   1.635            0.412         53.6 %
    кампания    530       4.454   1.189            0.315         44.0 %

That is item 83's error class -- a calibration fitted on one population and applied to another --
and on CYP3A4 three things flow through it: the trunk's screening channel, the band model the dead
zone reprojects onto, and the estimate of delta.

**The direct test cannot be run, and the reason is the concern itself.** The calibration residual
is `log2(1-I) - показание`, so it needs a screening reading; the 530 campaign compounds have none.
There is no population on which to measure the misfit, because the missing reading is what defines
the population. Any script claiming to measure it directly is measuring something else.

What is buildable is the question one step back: **how population-dependent is (E, h) at all?** If
splitting the screened set along the very axes on which the campaign differs from it moves the
constants and costs residual when they are transported across the split, then extrapolating to a
population further away than the split is unjustified. If the constants transport within the
screened set, the extrapolation is ordinary and item 129 remains a composition fact with no
consequence for the instrument.

The split axes are chosen to be the ones the campaign actually differs on, not arbitrary ones:
band width (campaign 0.315 against 0.412) and label level (4.454 against 4.203). Each is split at
its median inside the screened set, constants are refitted on each half, and each half's constants
are applied to the other. The excess residual on transport is the quantity of interest.

The last block reports how far the campaign sits from the screened set on those same axes measured
in units of the split just performed, so "further away than the split" is a number rather than an
assertion.

Reads the inhibition table and the screening table. Writes nothing. Seconds.
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


def fit(pi, y):
    """Та же подгонка, что в verify/g1_calib.py, скопированная без изменений."""
    def resid(t):
        E, h = t
        I = E / (1 + 10 ** (h * (PC0 - pi)))
        return np.log2(np.clip(1 - I, 1e-3, None)) - y
    r = least_squares(resid, [1.0, 1.0], bounds=([0.2, 0.2], [1.2, 4.0]))
    return r.x[0], r.x[1]


def resid_of(E, h, pi, y):
    I = E / (1 + 10 ** (h * (PC0 - pi)))
    return np.log2(np.clip(1 - I, 1e-3, None)) - y


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    piv = sc.pivot_table(index="Molecule_Name", columns="enzyme", values="log2fc_estimate")
    m = tr.set_index("Molecule_Name").join(piv)

    print("=== 1. две популяции есть только у CYP3A4 ===")
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        lab = m[col].notna()
        print(f"  {c}: помечено {int(lab.sum()):5d}, из них со скринингом "
              f"{int((lab & m[c].notna()).sum()):5d}, БЕЗ скрининга "
              f"{int((lab & m[c].isna()).sum()):5d}")

    print("\n=== 2. константы на всём скринированном наборе (должны лечь в g1_calib) ===")
    full = {}
    print(f"  {'фермент':8s} {'n':>5s} {'E':>7s} {'h':>7s} {'ско остатка':>12s}")
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        k = m[col].notna() & m[c].notna()
        pi, y = m.loc[k, col].to_numpy(), m.loc[k, c].to_numpy()
        E, h = fit(pi, y)
        full[c] = (E, h, pi, y, k)
        print(f"  {c:8s} {int(k.sum()):5d} {E:7.3f} {h:7.3f} "
              f"{resid_of(E, h, pi, y).std():12.3f}")

    print("\n=== 3. переносится ли калибровка между популяциями ВНУТРИ скрина ===")
    print("  Раскол по осям, на которых кампания отличается от библиотеки.")
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        E0, h0, pi, y, k = full[c]
        band = (m.loc[k, col + "_conf_high"] - m.loc[k, col + "_conf_low"]).to_numpy()
        axes = {"ширина полосы": band, "уровень метки": pi}
        print(f"\n  {c}")
        for nm, v in axes.items():
            lo = v <= np.median(v)
            EA, hA = fit(pi[lo], y[lo])
            EB, hB = fit(pi[~lo], y[~lo])
            # свои константы против чужих, на КАЖДОЙ половине
            sAA = resid_of(EA, hA, pi[lo], y[lo]).std()
            sBA = resid_of(EB, hB, pi[lo], y[lo]).std()
            sBB = resid_of(EB, hB, pi[~lo], y[~lo]).std()
            sAB = resid_of(EA, hA, pi[~lo], y[~lo]).std()
            print(f"    по «{nm}»: E {EA:.3f}/{EB:.3f}  h {hA:.3f}/{hB:.3f}  "
                  f"| ско свои {sAA:.3f}/{sBB:.3f}  чужие {sBA:.3f}/{sAB:.3f}  "
                  f"| перенос стоит {max(sBA - sAA, sAB - sBB):+.3f}")

    print("\n=== 4. насколько кампания дальше, чем этот раскол ===")
    c = "CYP3A4"
    col = f"{c}_pIC50_direct_inhibition"
    lab = m[col].notna()
    ins, out = lab & m[c].notna(), lab & m[c].isna()
    for nm, series in (("уровень метки", m[col]),
                       ("ширина полосы", m[col + "_conf_high"] - m[col + "_conf_low"])):
        a = series[ins].to_numpy(); b = series[out].to_numpy()
        half_lo = a[a <= np.median(a)]; half_hi = a[a > np.median(a)]
        d_split = abs(np.median(half_hi) - np.median(half_lo))
        d_camp = abs(np.median(b) - np.median(a))
        print(f"  {nm:15s} раскол внутри скрина {d_split:.3f}, "
              f"кампания от скрина {d_camp:.3f}  --- в {d_camp/max(d_split,1e-9):.2f} раза")

    print("""
Как читать. Решает строка «перенос стоит» в разделе 3 на CYP3A4. Ноль означает, что
константы --- свойство прибора: их можно подогнать на одной популяции и применять к другой,
и опасение снимается. Заметная величина означает, что они свойство популяции, и тогда
применение констант с 1805 скринированных к 530 из кампании ничем не оправдано, а через них
на CYP3A4 идут скрининговый канал ствола, модель полосы под мёртвой зоной и оценка дельты.

Раздел 4 переводит «дальше» в число: если кампания отстоит от скрина на величину порядка
раскола или больше, то раздел 3 --- нижняя оценка того, что происходит при экстраполяции,
а не её измерение.

Остальные три фермента здесь --- контроль стенда: у них популяция одна, поэтому перенос
между половинами должен быть дешевле, чем у CYP3A4, если дело действительно в составе.""")


if __name__ == "__main__":
    main()
