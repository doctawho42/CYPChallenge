"""Four seeds for the per-enzyme offset in the trunk, paired against the coupling without it.

Item 172 measured `--mode calshift` on two seeds and found the offset repairs the scale damage the
coupling causes: at lambda 0.3 the metric improved by 0.161 and 0.069 while rank moved 0.003 and
0.004, inside the floor. Two seeds is not four, and this file's own rule since item 119 is that a
seed-0 or seed-pair result is a hypothesis.

What is compared. `calibrated` and `calshift` differ by exactly one parameter vector -- four
numbers, initialised at zero -- and share weights initialisation, folds, standardisation and the
screening mask, because initialisation is seeded on (seed, fold) and not on the mode. At lambda 0
they are the same model and must agree exactly; that identity is checked here rather than assumed,
and if it fails nothing below can be read.

The test is paired across seeds, because the seeds are the repeated measure and the pairing removes
the split variance that item 165 measured at 0.0036 macro and up to 0.0071 per enzyme.

Two questions, separately answered:

  **сцепление против его отсутствия**   lambda > 0 против lambda 0, внутри calshift
  **сдвиг против его отсутствия**       calshift против calibrated при равной lambda

The first asks whether the screening channel is worth having at all once the offset exists. The
second asks what the offset itself buys. They are different questions and the two-seed item ran
them together.

The fitted offsets are reported across all four seeds, with their spread, because item 172's
reading -- that they carry shrinkage compensation and not only an assay offset -- predicts they
should be **stable across seeds** (shrinkage is a property of the model and the enzyme) rather than
wandering. That is a check the two-seed run could not make.

Reads results/preds/trunk_calibrated.json and trunk_calshift.json. Writes nothing.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import RES

import json

import numpy as np
import pandas as pd
from scipy import stats

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
FLOOR = {"CYP1A2": 0.0061, "CYP2C9": 0.0071, "CYP2D6": 0.0049, "CYP3A4": 0.0033}


def paired(a, b, name):
    """Парный по сидам тест: разность, t, p и сколько знаков совпало."""
    d = np.asarray(a) - np.asarray(b)
    if len(d) < 2 or np.allclose(d, 0):
        return f"{name}: разность ровно ноль ({len(d)} сидов)"
    t, p = stats.ttest_rel(a, b)
    sg = int(np.sum(np.sign(d) == np.sign(np.mean(d))))
    return (f"{name}: {np.mean(d):+.4f}  t = {t:+.2f}  p = {p:.4f}  "
            f"знаков {sg}/{len(d)}")


def main():
    cal = pd.DataFrame(json.load(open(RES + "preds/trunk_calibrated.json"))["table"])
    sh = pd.DataFrame(json.load(open(RES + "preds/trunk_calshift.json"))["table"])
    seeds = sorted(set(cal.seed) & set(sh.seed))
    lams = sorted(set(cal["lambda"]) & set(sh["lambda"]))
    print(f"общих сидов {len(seeds)}: {seeds} | общих lambda: {lams}\n")

    # --- контроль тождества при lambda = 0 ---
    z0 = cal[cal["lambda"] == 0].set_index("seed").reindex(seeds)
    z1 = sh[sh["lambda"] == 0].set_index("seed").reindex(seeds)
    dm = np.abs(z0.MACRO.to_numpy() - z1.MACRO.to_numpy()).max()
    dr = np.abs(z0.MACRO_rho.to_numpy() - z1.MACRO_rho.to_numpy()).max()
    print(f"КОНТРОЛЬ при lambda 0: максимум |разности| пара {dm:.6f}, ранг {dr:.6f}")
    if max(dm, dr) > 1e-6:
        print("  РАСХОЖДЕНИЕ. При lambda 0 модели обязаны совпадать; дальше читать нельзя.")
    else:
        print("  совпало --- сдвиг без скринингового члена не работает, как и задумано.\n")

    print("=== 1. стоит ли канал того, ВНУТРИ calshift (lambda > 0 против lambda 0) ===")
    base = sh[sh["lambda"] == 0].set_index("seed").reindex(seeds)
    for lam in [l for l in lams if l > 0]:
        cur = sh[sh["lambda"] == lam].set_index("seed").reindex(seeds)
        print(f"  lambda {lam}")
        print("    " + paired(cur.MACRO_rho, base.MACRO_rho, "ранг макро"))
        print("    " + paired(cur.MACRO, base.MACRO, "пара макро"))
        for c in CYPS:
            d = cur[f"rho_{c}"].to_numpy() - base[f"rho_{c}"].to_numpy()
            mark = "*" if abs(np.mean(d)) > FLOOR[c] else " "
            print(f"      {c}: {np.mean(d):+.4f}{mark} знаков "
                  f"{int(np.sum(np.sign(d) == np.sign(np.mean(d))))}/{len(d)}")
    print("  * --- больше СОБСТВЕННОГО пола фермента (пункт 165)\n")

    print("=== 2. что покупает САМ сдвиг (calshift против calibrated при равной lambda) ===")
    for lam in [l for l in lams if l > 0]:
        a = sh[sh["lambda"] == lam].set_index("seed").reindex(seeds)
        b = cal[cal["lambda"] == lam].set_index("seed").reindex(seeds)
        print(f"  lambda {lam}")
        print("    " + paired(a.MACRO, b.MACRO, "пара макро (меньше лучше)"))
        print("    " + paired(a.MACRO_rho, b.MACRO_rho, "ранг макро"))
        for c in CYPS:
            print(f"      {c} пара {b[c].mean():.4f} -> {a[c].mean():.4f} "
                  f"({a[c].mean()-b[c].mean():+.4f})")
    print("""
Как читать. Раздел 1 отвечает, нужен ли скрининговый канал вообще, когда сдвиг уже есть.
Раздел 2 --- что стоит сам сдвиг при том же канале; там ожидается движение ПАРЫ при почти
неподвижном ранге, потому что чинится шкала, а не порядок.

Знаки важнее величин: 4/4 при малой разности сильнее, чем 3/4 при большой, потому что пол
пункта 165 намерен именно на посидовом разбросе.""")


if __name__ == "__main__":
    main()
