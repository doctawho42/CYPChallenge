"""Degrade the screening channel in steps and watch the effect move, enzyme by enzyme.

Why this exists. The per-enzyme effect of the joint likelihood is ordered exactly by how
well the screen tracks pIC50 - Spearman -1.000 against 0.936 / 0.896 / 0.862 / 0.828. That
is four points, and at n = 4 a perfect rank correlation has permutation probability 1/12.
Worse, everything differs between enzymes at once: label count, spread, band width, how
well the Hill calibration fits. The screening correlation covaries with all of it, so the
ordering is one bit of evidence and cannot be more.

Injecting noise turns the four points into four curves. Within one enzyme nothing varies
except how much the channel knows, because the compounds, the labels, the bands, the split
and the weight initialisation are all held fixed - so a monotone response inside an enzyme
cannot be explained by anything that distinguishes enzymes.

The dose variable. Noise is added to the STANDARDISED screening target and the result is
divided by sqrt(1 + eta^2), so the target keeps sd 1 and only its correlation with pIC50
falls, by exactly that factor. Without the division the target's sd would grow and lambda
would silently change meaning; "less information" would be confounded with "more weight".
The draw is shared across modes and lambdas, so arms compared at one eta see the same
corrupted channel, and the ladder is nested rather than five independent draws.

Three questions, in increasing order of how much they would settle.

1. Is the response monotone inside each enzyme? A real effect moves smoothly with the dose;
   noise gives a flat or ragged curve. This is the check the four-point ordering could not
   do.

2. Does the size of the movement order the enzymes the way the information story says? If
   the channel's contribution is what noise destroys, the enzymes with the most to lose -
   the ones whose screen is most informative - should move furthest.

3. Do the curves COLLAPSE? This is the strong one. If the effect is a function of channel
   informativeness alone, then plotting all four enzymes at all five doses against the
   effective informativeness rho_e / sqrt(1 + eta^2) should put twenty points on one curve.
   If they do not collapse, enzyme identity matters beyond informativeness and the
   single-factor reading is wrong however well the ordering fits.

Reads results/preds/trunk_{mode}.json for the eta = 0 rung and trunk_{mode}_noise{eta}.json
for the rest. Writes nothing.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import RES, tutorial
tutorial()

import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
ETAS = [0.0, 0.5, 1.0, 2.0, 4.0]
MODES = ["twohead", "calibrated"]
LAM = 3.0

# Measured in src/trunkscore.py: |Spearman| between the screen and pIC50, per enzyme.
RHO_SCR = {"CYP1A2": 0.862, "CYP2C9": 0.896, "CYP2D6": 0.828, "CYP3A4": 0.936}


def load(mode, eta):
    tag = "" if eta == 0 else f"_noise{eta:g}"
    f = _pl.Path(RES + f"preds/trunk_{mode}{tag}.json")
    if not f.exists():
        return None, None
    b = json.load(open(f))
    return pd.DataFrame(b["table"]), b["preds"]


def cell(tab, seed, lam, c):
    r = tab[(tab.seed == seed) & (tab["lambda"] == lam)]
    return float(r[c].iloc[0]) if len(r) else np.nan


def main():
    base = {m: load(m, 0.0) for m in MODES}
    for m in MODES:
        if base[m][0] is None:
            raise SystemExit(f"нет results/preds/trunk_{m}.json - сначала src/trunk.py")
    seeds = sorted(base["twohead"][0].seed.unique())

    print("=" * 96)
    print("0. Контроль: при lambda = 0 канал не входит в потери, шум обязан ничего не менять")
    print("=" * 96)
    t2, p2 = load("twohead", 2.0)
    if t2 is not None and (t2["lambda"] == 0.0).any():
        worst = 0.0
        for s in seeds:
            a = np.asarray(base["twohead"][1][f"twohead|{s}|0.0"], float)
            b = np.asarray(p2[f"twohead|{s}|0.0"], float)
            m = ~np.isnan(a) & ~np.isnan(b)
            worst = max(worst, float(np.max(np.abs(a[m] - b[m]))))
        print(f"  eta = 2 против eta = 0 при lambda = 0, макс |разность| по четырём сидам: "
              f"{worst:.3e}")
        print("  " + ("совпало побитово, как и должно" if worst == 0
                      else "РАСХОЖДЕНИЕ: шум течёт туда, куда не должен"))
    else:
        print("  прогон lambda = 0 при eta = 2 не найден, контроль пропущен")

    print()
    print("=" * 96)
    print(f"1. Лестница: эффект канала при lambda = {LAM} против lambda = 0, среднее по сидам")
    print("=" * 96)
    print("   (отрицательное = канал помогает; eta - порча в ско канала, ранг падает в")
    print("    1/sqrt(1+eta^2) раз)\n")

    eff = {}       # (mode, enzyme, eta) -> list по сидам
    missing = []
    for mode in MODES:
        tb, _ = base[mode]
        for eta in ETAS:
            t, _ = load(mode, eta)
            if t is None:
                missing.append((mode, eta))
                continue
            for c in CYPS:
                eff[(mode, c, eta)] = [cell(t, s, LAM, c) - cell(tb, s, 0.0, c)
                                       for s in seeds]
    if missing:
        print(f"  не найдены прогоны: {missing}\n")
    # Сорванный прогон - это результат, а не помеха: одноголовая рука разошлась численно
    # на отдельных ячейках. Считаем их и называем, а не усредняем молча по остатку.
    dead = [(m, c, e, i) for (m, c, e), v in eff.items()
            for i, x in enumerate(v) if not np.isfinite(x)]
    if dead:
        print(f"  РАЗОШЛИСЬ И ДАЛИ NaN: {len(dead)} прогонов из {len(eff) * len(seeds)}")
        for m, c, e, i in dead:
            print(f"    {m:11s} {c} eta={e:g} сид={seeds[i]}")
        print("  Ниже они исключены из средних; там, где исключать пришлось, стоит звёздочка.\n")

    att = {e: 1.0 / np.sqrt(1.0 + e * e) for e in ETAS}
    for mode in MODES:
        print(f"--- {mode} ---")
        head = f"{'фермент':8s} {'|rho|':>6s}"
        for e in ETAS:
            head += f"  eta={e:<4g}"
        print(head + f" {'сдвиг 0->4':>11s} {'Спирмен':>8s}")
        for c in CYPS:
            row = f"{c:8s} {RHO_SCR[c]:6.3f}"
            vals = []
            for e in ETAS:
                v = eff.get((mode, c, e))
                if v is None:
                    row += f"  {'--':>7s} "; vals.append(np.nan)
                else:
                    mv = float(np.nanmean(v)); star = "*" if not np.all(np.isfinite(v)) else " "
                    row += f" {mv:+7.3f}{star}"; vals.append(mv)
            ok = ~np.isnan(vals)
            mv = vals[int(np.where(ok)[0][-1])] - vals[0] if ok.sum() > 1 else np.nan
            rr = (spearmanr(np.array(ETAS)[ok], np.array(vals)[ok]).statistic
                  if ok.sum() > 2 else np.nan)
            print(row + f" {mv:+11.3f} {rr:+8.3f}")
        print(f"{'затухание ранга':15s} " + "  ".join(f"{att[e]:7.3f}" for e in ETAS))
        print()

    print("=" * 96)
    print("2. Держится ли знак сдвига на каждом сиде по отдельности")
    print("=" * 96)
    top = max(e for e in ETAS if any((m, CYPS[0], e) in eff for m in MODES))
    print(f"   сдвиг = эффект при eta = {top:g} минус эффект при eta = 0, по сидам\n")
    for mode in MODES:
        print(f"--- {mode} ---")
        print(f"{'фермент':8s} " + " ".join(f"{('сид ' + str(s)):>9s}" for s in seeds)
              + f" {'среднее':>9s} {'знак 4/4':>9s}")
        for c in CYPS:
            a, b = eff.get((mode, c, 0.0)), eff.get((mode, c, top))
            if a is None or b is None:
                continue
            d = [x - y for x, y in zip(b, a)]
            fin = [x for x in d if np.isfinite(x)]
            same = len(set(np.sign(fin))) == 1 if fin else False
            print(f"{c:8s} " + " ".join((f"{x:+9.3f}" if np.isfinite(x) else f"{'—':>9s}") for x in d)
                  + f" {np.nanmean(d):+9.3f} {str(same):>9s}")
        print()

    print("=" * 96)
    print("3. Упорядочивает ли информативность величину сдвига")
    print("=" * 96)
    print("   Если шум уничтожает вклад канала, дальше всех уедет фермент, которому есть")
    print("   что терять, то есть с самым информативным скринингом.\n")
    for mode in MODES:
        mv, rho = [], []
        for c in CYPS:
            a, b = eff.get((mode, c, 0.0)), eff.get((mode, c, top))
            if a is None or b is None:
                continue
            mv.append(np.nanmean(b) - np.nanmean(a)); rho.append(RHO_SCR[c])
        if len(mv) == 4:
            r = spearmanr(rho, mv).statistic
            order = " > ".join(c[3:] for c in sorted(
                CYPS, key=lambda c: -(np.nanmean(eff[(mode, c, top)]) - np.nanmean(eff[(mode, c, 0.0)]))))
            print(f"  {mode:12s} Спирмен(|rho|, сдвиг) = {r:+.3f}   порядок сдвига: {order}")
    print("\n  n = 4 снова, но теперь каждая точка сама опирается на кривую, а не на одно")
    print("  измерение, так что случайное совпадение порядка стоит дороже.")

    print()
    print("=" * 96)
    print("4. Схлопываются ли кривые: эффект против эффективной информативности")
    print("=" * 96)
    print("   Если эффект есть функция одной только информативности канала, то все точки")
    print("   ложатся на одну кривую по rho_e / sqrt(1 + eta^2), независимо от фермента.\n")
    for mode in MODES:
        xs, ys, tag = [], [], []
        for c in CYPS:
            for e in ETAS:
                v = eff.get((mode, c, e))
                if v is None:
                    continue
                y_ = float(np.nanmean(v))
                if not np.isfinite(y_):
                    continue
                xs.append(RHO_SCR[c] * att[e]); ys.append(y_); tag.append(c)
        if len(xs) < 8:
            continue
        xs, ys = np.asarray(xs), np.asarray(ys)
        r = spearmanr(xs, ys).statistic
        # Насколько хорошо одна монотонная кривая объясняет все точки сразу.
        fit = IsotonicRegression(out_of_bounds="clip", increasing="auto").fit(xs, ys)
        res = ys - fit.predict(xs)
        share = 1.0 - float(res.var() / ys.var()) if ys.var() > 0 else np.nan
        print(f"--- {mode} ---")
        print(f"  точек {len(xs)}, Спирмен(эффективная информативность, эффект) = {r:+.3f}")
        print(f"  одна монотонная кривая объясняет {share*100:5.1f} % дисперсии эффекта,")
        print(f"  остаток ско {res.std():.4f} при полном ско {ys.std():.4f}")
        by = {}
        for c, rr in zip(tag, res):
            by.setdefault(c, []).append(rr)
        print("  средний остаток по ферментам: "
              + "  ".join(f"{c[3:]} {np.mean(v):+.3f}" for c, v in by.items()))
        print()

    print("""Как это читать. Пункт 1 отвечает на то, ради чего всё затевалось: монотонный отклик
внутри одного фермента нельзя объяснить ничем, что различает ферменты, потому что внутри
фермента не меняется ничего, кроме дозы. Пункт 4 - самый строгий: если остатки вокруг
единой кривой систематически смещены по ферментам, значит информативность канала не
единственное, что определяет эффект, и одного фактора мало, как бы хорошо ни ложился
порядок из четырёх точек.""")


if __name__ == "__main__":
    main()
