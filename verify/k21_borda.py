"""Average the members' ranks instead of their values.

The argument. What survives the affine pair is rank (item 77), and the pair is only a
two-parameter monotone family fitted afterwards. So the ensemble's job is to produce a good
ordering, and the scale of any one member is not information -- it is something the pair will
overwrite anyway. Averaging values lets a member with a wider spread pull the joint answer
around in proportion to that spread, which is a weight nobody chose.

Averaging ranks removes that. Each member contributes its ordering and nothing else, the
averaged rank is turned back into pIC50 through the empirical quantile function of the
training labels, and the affine pair is fitted on top exactly as before.

This is a cheap test with a clean falsification: if value-averaging and rank-averaging agree
to within the 0.007 noise floor, then no member's scale was distorting anything and the
question is closed. If they differ, the members' spreads are not commensurate and every
ensemble number in the log was averaging apples with pears.

Reads results/preds/oof_pool_all.json, oof_gp.json, oof_weak.json. Writes nothing.
Four seeds, under a minute -- everything is already computed.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, rankdata
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
SEEDS = (0, 1, 2, 3)
MEMBERS = [("oof_pool_all", "независимо"), ("oof_pool_all", "пул"),
           ("oof_gp", "GP"), ("oof_weak", "гребневая")]


def borda(ps, y_ref):
    """Mean of the members' ranks, mapped back through the label quantile function.

    y_ref is the training label sample whose quantiles the ranks are read against. Using
    the labels rather than any member's own values is what makes the result scale-free.
    """
    r = np.mean([rankdata(p) for p in ps], axis=0)
    q = (rankdata(r) - 0.5) / len(r)
    return np.quantile(y_ref, q)


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    cache = {f: json.load(open(RES + f"preds/{f}.json"))["preds"]
             for f in {m[0] for m in MEMBERS}}

    table = []
    for seed in SEEDS:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for n_use in (2, 3, 4):
            r = {"seed": seed, "членов": n_use}
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                fi, w = fold[m], np.ones(m.sum()) / m.sum()

                ps = [np.asarray(cache[fn][f"{seed}|{arm}|{c}"], float)
                      for fn, arm in MEMBERS[:n_use]]
                mean_p = np.mean(ps, axis=0)
                bord_p = borda(ps, y)

                for tag, p in (("значения", mean_p), ("ранги", bord_p)):
                    q = fit_apply(p, lo, hi, fi, w)
                    r[f"{c} {tag}"] = float(strae(y, q, y_true_upper=hi, y_true_lower=lo))
                    r[f"rho {c} {tag}"] = float(spearmanr(y, p).statistic)
            for tag in ("значения", "ранги"):
                r[f"MACRO {tag}"] = float(np.mean([r[f"{c} {tag}"] for c in CYPS]))
                r[f"MACRO rho {tag}"] = float(np.mean([r[f"rho {c} {tag}"] for c in CYPS]))
            table.append(r)
            print(f"  сид {seed} {n_use} членов: "
                  f"значения {r['MACRO значения']:.4f}/{r['MACRO rho значения']:.4f}  "
                  f"ранги {r['MACRO ранги']:.4f}/{r['MACRO rho ранги']:.4f}", flush=True)

    df = pd.DataFrame(table)
    print()
    print("Среднее по четырём сидам (пара ST-RAE / ранг):")
    print(f"{'членов':>7s} {'значения':>18s} {'ранги':>18s} {'разница пары':>14s} {'разница ранга':>15s}")
    for n_use in (2, 3, 4):
        d = df[df["членов"] == n_use]
        a, b = d["MACRO значения"].mean(), d["MACRO ранги"].mean()
        ra, rb = d["MACRO rho значения"].mean(), d["MACRO rho ранги"].mean()
        print(f"{n_use:7d} {a:9.4f}/{ra:.4f} {b:9.4f}/{rb:.4f} {b-a:+14.4f} {rb-ra:+15.4f}")

    print()
    print("По ферментам, четыре члена:")
    d = df[df["членов"] == 4]
    for c in CYPS:
        print(f"  {c}: пара {d[f'{c} значения'].mean():.4f} -> {d[f'{c} ранги'].mean():.4f}"
              f"   ранг {d[f'rho {c} значения'].mean():.4f} -> {d[f'rho {c} ранги'].mean():.4f}")

    print("""
Как читать. Порог --- шумовой пол 0.007. Разница меньше него не является результатом.

Ранг у двух схем может отличаться только за счёт того, что усреднение значений и усреднение
рангов дают разный ПОРЯДОК: сама по себе замена шкалы ранг не меняет. Так что ненулевая
разница в rho --- это прямое доказательство, что шкалы членов несоизмеримы, а нулевая
закрывает вопрос.""")


if __name__ == "__main__":
    main()
