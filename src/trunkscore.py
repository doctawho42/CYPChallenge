"""Score the lambda_scr response, with a paired bootstrap against the lambda = 0 control.

Follows the house rule from README.md: a difference without an interval is not a result,
comparisons run on one split at a time, and a decision needs the sign to hold across four
seeds. The pairing here is unusually tight - for a given (seed, fold) the two arms share
the split, the weight initialisation and the batch order, so the only thing that differs
between them is whether the screening term was in the loss.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json

import numpy as np
import pandas as pd
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
B = 2000


def truth():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    lo = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    hi = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    return y, lo, hi


def macro(y, lo, hi, p):
    return float(np.mean([
        strae(y[m, e], p[m, e], y_true_upper=hi[m, e], y_true_lower=lo[m, e])
        for e in range(4) for m in [~np.isnan(y[:, e]) & ~np.isnan(p[:, e])]
    ]))


def paired_macro_boot(y, lo, hi, pa, pb, rng):
    """Resample compounds once, then take each enzyme's labelled subset of that resample.

    Resampling the common table rather than each enzyme independently is what preserves
    the correlation between enzymes; four independent bootstraps would understate the
    interval on the macro.
    """
    n = y.shape[0]
    d = np.empty(B)
    for b in range(B):
        i = rng.integers(0, n, n)
        va, vb = [], []
        for e in range(4):
            m = ~np.isnan(y[i, e])
            if m.sum() < 30:
                continue
            j = i[m]
            va.append(strae(y[j, e], pa[j, e], y_true_upper=hi[j, e], y_true_lower=lo[j, e]))
            vb.append(strae(y[j, e], pb[j, e], y_true_upper=hi[j, e], y_true_lower=lo[j, e]))
        d[b] = np.mean(vb) - np.mean(va)
    return d


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inp", default=RES + "preds/trunk_twohead.json")
    a = ap.parse_args()

    blob = json.load(open(a.inp))
    tab = pd.DataFrame(blob["table"])
    # ключ теперь {режим}|{сид}|{lambda}; режим один на файл, поэтому срезаем его
    preds = {"|".join(k.split("|")[1:]): np.asarray(v, float)
             for k, v in blob["preds"].items()}
    y, lo, hi = truth()

    seeds = sorted(tab.seed.unique())
    lams = sorted(tab["lambda"].unique())
    base = min(lams)

    print("=== макро ST-RAE по сидам и lambda (ниже - лучше) ===")
    piv = tab.pivot(index="seed", columns="lambda", values="MACRO")
    print(piv.to_string())
    print(f"\nдля сравнения: бустинг на FP+DESC+MECH даёт 0.7673 на сиде 0")

    print(f"\n=== отклик: разность к lambda={base} на том же сиде ===")
    print("(отрицательное = канал скрининга помогает)")
    rows = []
    for lam in lams:
        if lam == base:
            continue
        d = [piv.loc[s, lam] - piv.loc[s, base] for s in seeds]
        rows.append({"lambda": lam, **{f"сид {s}": round(v, 4) for s, v in zip(seeds, d)},
                     "среднее": round(float(np.mean(d)), 4),
                     "знак одинаков": bool(np.all(np.sign(d) == np.sign(d[0])))})
    print(pd.DataFrame(rows).to_string(index=False))

    print(f"\n=== парный бутстрап по соединениям, {B} ресэмплов, сид разбиения 0 ===")
    rng = np.random.default_rng(7)
    pa = preds[f"0|{base}"]
    for lam in lams:
        if lam == base:
            continue
        d = paired_macro_boot(y, lo, hi, pa, preds[f"0|{lam}"], rng)
        print(f"  lambda {lam:<5} Δмакро {d.mean():+.4f} "
              f"[{np.percentile(d, 2.5):+.4f}, {np.percentile(d, 97.5):+.4f}]  "
              f"P(хуже) {float(np.mean(d > 0)):.3f}")

    print("\n=== по ферментам, среднее по сидам ===")
    out = []
    for lam in lams:
        r = {"lambda": lam}
        for c in CYPS:
            v = [tab[(tab.seed == s) & (tab["lambda"] == lam)][c].iloc[0] for s in seeds]
            r[c] = round(float(np.mean(v)), 4)
        r["MACRO"] = round(float(np.mean([r[c] for c in CYPS])), 4)
        r["rho"] = round(float(np.mean(
            [tab[(tab.seed == s) & (tab["lambda"] == lam)]["MACRO_rho"].iloc[0] for s in seeds])), 3)
        out.append(r)
    print(pd.DataFrame(out).to_string(index=False))

    print("""
Как это читать. Настоящий эффект меняется с lambda гладко и имеет максимум; шум даёт
плоскую или рваную кривую. Знак, одинаковый на четырёх сидах, и интервал бутстрапа,
не накрывающий ноль, - минимум, при котором разность в этом репозитории считается
результатом. Ноль здесь - тоже результат, и записывать его надо так же.""")


if __name__ == "__main__":
    main()
