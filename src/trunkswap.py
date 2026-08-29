"""Swap the instrument calibration between two enzymes and see whether the damage follows it.

The problem this exists to solve. Four quantities order the four enzymes identically, with
pairwise rank correlation of exactly plus or minus one: how well the screen tracks pIC50, the
mean width of the confidence band, the gain of the calibration curve, and the per-enzyme
damage itself. Nothing measured so far separates them, and noise cannot: degrading the channel
moves informativeness and leaves both the band and the instrument map exactly where they were.
The band was eliminated separately, by rescoring the same predictions with the band equalised
and then removed. Two candidates are left.

Swapping (E, h) between CYP2D6 and CYP3A4 is the one intervention that moves the instrument
map while holding everything else literally fixed - the same compounds, labels, bands, split,
initial weights, batch order, and the same screening readout with the same informativeness.

The fork is clean:

  if the damage follows the MAP, then 2D6 given 3A4's calibration should be hurt much less
  than its usual +0.63, and 3A4 given 2D6's should go from near zero to something large. The
  mechanism is then the gain of g in pIC50 units, and it is one mechanism, not two;

  if the damage stays with the ENZYME, the map is not what orders them, and the remaining
  candidate - how informative that enzyme's screen is - survives alone.

The experiment carries its own internal control, which none of the earlier ones did. CYP1A2
and CYP2C9 keep their own calibrations throughout. Their damage must not move under either
outcome; if it does, something other than the swap changed and the run is void.

A second control is free. At lambda = 0 the screening term is absent from the loss, so
g_of_pi is never called and the swap cannot reach the predictions. Those two arms must agree
bit for bit.

Reads results/preds/trunk_calibrated.json and trunk_calibrated_swap2D6-3A4.json. Writes
nothing. Produce the second with:

    uv run python src/trunk.py --mode calibrated --seeds 0,1,2,3 --lams 0,3.0 \\
        --swap-cal CYP2D6,CYP3A4
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()

import json

import numpy as np
import pandas as pd

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
SWAPPED = ("CYP2D6", "CYP3A4")
PC0 = 4.305
CAL_E = np.array([0.728, 0.621, 0.867, 0.931])
CAL_H = np.array([1.261, 1.112, 1.243, 1.968])


def max_slope(e_, h_):
    """Maximum of |dg/dpi| for g(pi) = log2(1 - E/(1 + 10^(h(pC0 - pi)))).

    Not attained at pC0, which is the easy thing to assume: the log2 outside the Hill term
    moves the peak up-scale by log10(sqrt(1 - E))/h. Solved on a dense grid instead of
    algebraically, because the closed form is easy to get subtly wrong and this costs nothing.
    """
    pi = np.linspace(PC0 - 6, PC0 + 6, 2_000_001)
    inh = e_ / (1.0 + np.power(10.0, h_ * (PC0 - pi)))
    g = np.log2(np.clip(1.0 - inh, 1e-12, None))
    return float(np.max(np.abs(np.gradient(g, pi))))


def load(tag=""):
    f = _pl.Path(RES + f"preds/trunk_calibrated{tag}.json")
    if not f.exists():
        raise SystemExit(f"нет {f} - см. докстринг")
    b = json.load(open(f))
    return pd.DataFrame(b["table"]), b["preds"]


def cell(t, s, lam, c):
    r = t[(t.seed == s) & (t["lambda"] == lam)]
    return float(r[c].iloc[0]) if len(r) else np.nan


def main():
    t0, p0 = load()
    t1, p1 = load("_swap2D6-3A4")
    seeds = sorted(set(t0.seed.unique()) & set(t1.seed.unique()))

    print("=" * 90)
    print("0. Контроль утечки: при lambda = 0 карта прибора не вызывается вовсе")
    print("=" * 90)
    worst = 0.0
    for s in seeds:
        a = np.asarray(p0[f"calibrated|{s}|0.0"], float)
        b = np.asarray(p1[f"calibrated|{s}|0.0"], float)
        m = ~np.isnan(a) & ~np.isnan(b)
        worst = max(worst, float(np.max(np.abs(a[m] - b[m]))))
    print(f"  макс |разность| по {len(seeds)} сидам: {worst:.3e}")
    print("  " + ("совпало побитово, обмен не течёт" if worst == 0
                  else "РАСХОЖДЕНИЕ: обмен достаёт туда, куда не должен, прогон недействителен"))

    print()
    print("=" * 90)
    print("1. Что предсказывает карта прибора")
    print("=" * 90)
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    piv = sc.pivot_table(index="Molecule_Name", columns="enzyme", values="log2fc_estimate")
    print(f"{'фермент':8s} {'макс |dg/dpi|':>14s} {'ско скрина':>11s} {'усиление':>9s}")
    gain = {}
    for e, c in enumerate(CYPS):
        ms = max_slope(CAL_E[e], CAL_H[e])
        sd = float(piv[c].std()) if c in piv else np.nan
        gain[c] = sd / ms
        print(f"{c:8s} {ms:14.3f} {sd:11.3f} {gain[c]:9.3f}")
    print("\nУсиление --- сколько единиц pIC50 ложной тяги даёт одно ско остатка скрининга.")
    print(f"Обмен меняет местами карты {SWAPPED[0]} и {SWAPPED[1]}, то есть меняет местами и")
    print(f"их усиления: {gain[SWAPPED[0]]:.3f} <-> {gain[SWAPPED[1]]:.3f}.")

    print()
    print("=" * 90)
    print("2. Куда пошла порча")
    print("=" * 90)
    # Средним здесь мерить нельзя, и это не абстрактная осторожность: в этом эксперименте
    # один сид уже трижды создавал эффект, которого на остальных трёх нет. Поэтому сдвиг
    # печатается по сидам, а решение принимается по правилу четырёх сидов.
    print(f"{'фермент':8s} {'карта':>7s} " + " ".join(f"{('сид ' + str(s)):>8s}" for s in seeds)
          + f" {'среднее':>8s} {'без сида 0':>11s} {'знак':>6s}  что это")
    dmg, per = {}, {}
    for e, c in enumerate(CYPS):
        a = [cell(t0, s, 3.0, c) - cell(t0, s, 0.0, c) for s in seeds]
        b = [cell(t1, s, 3.0, c) - cell(t1, s, 0.0, c) for s in seeds]
        v = [y - x for x, y in zip(a, b)]
        dmg[c] = (float(np.mean(a)), float(np.mean(b))); per[c] = v
        same = len(set(np.sign(np.round(v, 3)))) == 1
        who = "своя" if c not in SWAPPED else (SWAPPED[1] if c == SWAPPED[0] else SWAPPED[0])[3:]
        role = "КОНТРОЛЬ" if c not in SWAPPED else "обменян"
        print(f"{c:8s} {who:>7s} " + " ".join(f"{x:+8.3f}" for x in v)
              + f" {np.mean(v):+8.3f} {np.mean(v[1:]):+11.3f} {str(same):>6s}  {role}")

    # Порог шума берём по контролям, и тоже без сида 0, иначе он завышен тем же выбросом.
    ctrl = max(abs(np.mean(per[c][1:])) for c in CYPS if c not in SWAPPED)
    u, v_ = SWAPPED
    du, dv = float(np.mean(per[u][1:])), float(np.mean(per[v_][1:]))
    print(f"\nШум перезапуска по контрольным ферментам (без сида 0): {ctrl:.3f}.")
    print("Сдвиг у обменянных надо читать против него, а не против нуля.")

    print()
    print("=" * 90)
    print("3. Что из этого следует")
    print("=" * 90)
    ok_u = len(set(np.sign(np.round(per[u], 3)))) == 1
    ok_v = len(set(np.sign(np.round(per[v_], 3)))) == 1
    moved = abs(du) > 2 * ctrl and abs(dv) > 2 * ctrl and du < 0 < dv and ok_u and ok_v
    stayed = abs(du) <= ctrl and abs(dv) <= ctrl
    if moved:
        print("Порча пошла за картой. Значит поферментный порядок задаёт усиление калибровки")
        print("в единицах pIC50, а не информативность скрининга: механизм один, и он назван.")
    elif stayed:
        print("Порча осталась при ферменте. Значит карта прибора порядок не задаёт, и из двух")
        print("неразличимых величин выживает информативность скрининга.")
    else:
        print("Картина смешанная: сдвиг есть, но не той величины или не того знака, какой")
        print("предсказывает чистый перенос. Ниже числа, разбираться по ним.")
    print(f"\n  {u}: {dmg[u][0]:+.3f} -> {dmg[u][1]:+.3f}  (сдвиг без сида 0 {du:+.3f})")
    print(f"  {v_}: {dmg[v_][0]:+.3f} -> {dmg[v_][1]:+.3f}  (сдвиг без сида 0 {dv:+.3f})")
    print("""
Оговорка, которая остаётся при любом исходе. Ферментов по-прежнему четыре, а обменяны два,
так что это одно сравнение, а не выборка. Зато оно контролируемое: два необменянных фермента
показывают, сколько шума даёт сам перезапуск, и сдвиг у обменянных надо читать против этой
величины, а не против нуля.""")


if __name__ == "__main__":
    main()
