"""Recover the organisers' rule for choosing which compounds get a dose-response curve.

The puzzle. Curves were run on a third of the screened molecules, and which third depends on the
screening reading in a way that looks different on every enzyme: on CYP2D6 a hard threshold at
the active end, on CYP1A2 a softer one, on CYP2C9 no dependence at all, and on CYP3A4 the
opposite, with curves preferentially on the WEAKER readings. Three shapes, not two, and the
CYP3A4 direction is the one that has to be explained.

The hypothesis, and it makes all four the same rule. A curve is uninformative where the single
point already answers the question: if 49.5 uM gives near-complete inhibition, IC50 lies far
below the range and the curve returns "less than the lowest point". A curve is worth running
where the point sits on the SLOPE. That is exactly the Fisher information of a single reading
about pi, which section 5 of the document derives for another purpose:

    I(pi) = E h ln10 x / (1 + x)^2,   x = 10^{h (pC0 - pi)},   maximal at x = 1, i.e. pi = pC0

The apparent reversal is then an artefact of the axis. On CYP2D6 the screened population is not
saturated, so the slope sits at the active end and the rule looks like "take the actives". On
CYP3A4 the population is deep in saturation - median log2fc of the compounds that DID get curves
is -3.301 in the top tercile - so the slope has moved to the weak end of the ranking and the same
rule looks inverted. On CYP2C9 the whole population sits on the slope and the rule flattens.

The test. Invert every screening reading to pi through the fitted (E, h), compute its Fisher
weight, and plot P(curve was run) against the WEIGHT rather than against the rank. If the rule is
what is claimed, the four curves collapse onto one. Collapse is measured, not eyeballed: the
between-enzyme spread of P at matched weight is compared with the same spread at matched rank.

If it collapses, the experimental design has been recovered from the data - which also settles
what the 11509 unlabelled readings are worth, since by the same rule they are precisely the pool
from which least can be extracted.

Reads the screening file, the inhibition file and rows.csv. Prints only.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D

import numpy as np
import pandas as pd

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
E = {"CYP1A2": 0.728, "CYP2C9": 0.621, "CYP2D6": 0.867, "CYP3A4": 0.931}
H = {"CYP1A2": 1.261, "CYP2C9": 1.112, "CYP2D6": 1.243, "CYP3A4": 1.968}
PC0 = 4.305


def fisher(g, e, h):
    """Fisher weight of one reading about pi, normalised to its own maximum.

    Both saturations send x to 0 or to infinity and the weight to zero, which is the whole
    point, so the non-invertible readings need no special case beyond clamping.
    """
    I = 1.0 - np.power(2.0, g)
    I = np.clip(I, 1e-6, e - 1e-6)
    x = e / I - 1.0
    w = e * h * np.log(10.0) * x / np.square(1.0 + x)
    return w / (e * h * np.log(10.0) * 0.25)


def main():
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())

    dat = {}
    for c in CYPS:
        have = set(tr.loc[tr[f"{c}_pIC50_direct_inhibition"].notna(), "Molecule_Name"])
        s = sc[sc.enzyme == c].copy()
        s["curve"] = s.Molecule_Name.isin(have).astype(float)
        s["w"] = fisher(s.log2fc_estimate.to_numpy(float), E[c], H[c])
        dat[c] = s

    print("P(снята кривая) по децилям скрининга, дециль 1 = самые активные")
    print(f"{'дециль':>7s} " + " ".join(f"{c:>8s}" for c in CYPS))
    Pr = np.zeros((10, 4))
    for c in CYPS:
        s = dat[c]
        q = pd.qcut(s.log2fc_estimate.rank(method="first"), 10, labels=False)
        Pr[:, CYPS.index(c)] = s.groupby(q).curve.mean().to_numpy()
    for i in range(10):
        print(f"{i+1:7d} " + " ".join(f"{Pr[i, j]:8.3f}" for j in range(4)))

    print("\nP(снята кривая) по фишеровскому весу одноточечного показания")
    print(f"{'вес':>12s} " + " ".join(f"{c:>8s}" for c in CYPS) + f" {'n':>7s}")
    edges = np.array([0, .02, .05, .1, .2, .35, .5, .7, .85, .95, 1.001])
    Pw = np.full((len(edges) - 1, 4), np.nan)
    for c in CYPS:
        s = dat[c]
        b = np.digitize(s.w, edges) - 1
        for i in range(len(edges) - 1):
            m = b == i
            if m.sum() >= 25:
                Pw[i, CYPS.index(c)] = s.curve[m].mean()
    for i in range(len(edges) - 1):
        n = sum(int(((np.digitize(dat[c].w, edges) - 1) == i).sum()) for c in CYPS)
        cells = " ".join("     ---" if np.isnan(Pw[i, j]) else f"{Pw[i, j]:8.3f}" for j in range(4))
        print(f"{edges[i]:5.2f}-{edges[i+1]:5.2f} {cells} {n:7d}")

    sd_rank = np.nanmean(np.nanstd(Pr, axis=1))
    sd_w = np.nanmean(np.nanstd(Pw, axis=1))
    print(f"\nразброс между ферментами при равном РАНГЕ: {sd_rank:.3f}")
    print(f"разброс между ферментами при равном ВЕСЕ:  {sd_w:.3f}")
    print(f"схлопывание: {'ДА' if sd_w < 0.7 * sd_rank else 'НЕТ'}, "
          f"отношение {sd_w / sd_rank:.2f}")

    print(f"\n{'фермент':8s} {'медиана log2fc':>15s} {'снятых':>8s} {'несnятых':>9s} "
          f"{'медиана веса снятых':>20s} {'несnятых':>9s}")
    for c in CYPS:
        s = dat[c]
        a, b = s[s.curve == 1], s[s.curve == 0]
        print(f"{c:8s} {'':15s} {a.log2fc_estimate.median():8.3f} "
              f"{b.log2fc_estimate.median():9.3f} {a.w.median():20.3f} {b.w.median():9.3f}")

    print("""
Как читать. Верхняя таблица --- то, что видно глазом: четыре разных рисунка. Нижняя --- то же
самое по фишеровскому весу. Если гипотеза верна, четыре столбца нижней таблицы должны стать
похожи, а разброс между ферментами --- упасть.

Последняя таблица проверяет механизм напрямую: медиана веса у снятых должна быть ВЫШЕ, чем у
несnятых, на всех четырёх ферментах сразу --- включая CYP3A4, где по сырому log2fc всё выглядит
наоборот.""")


if __name__ == "__main__":
    main()
