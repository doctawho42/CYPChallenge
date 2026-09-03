"""Can the trunk's screening channel be moved into the ensemble.

The state item 118 left. `src/trunk.py`'s screening channel is worth -0.0264 after the affine
pair with the sign holding on all four seeds (item 79), which makes it the largest surviving
intervention in this log. But the trunk carrying it only ties with the boosting -- 0.7149
against 0.7155, sign alternating, p = 0.72 -- so the vehicle is shelved and the channel is
stranded inside it.

Ties are exactly the situation item 100 was about. The ridge is worse than the boosting alone
and improves the ensemble anyway, because it errs smoothly where the boosting errs locally. A
model that *ties* with the boosting while coming from an entirely different family is a better
prior for that than the ridge was, and it costs nothing to check: every prediction involved is
already on disk.

Two arms, four seeds:

    current four   per-enzyme boosting, pooled boosting, GP on descriptors, ridge on descriptors
    plus trunk     the same four with the two-head trunk at lambda = 3
    control        the same four with the trunk at lambda = 0 -- same architecture, same
                   parameter count, same initial weights, screening head receiving no
                   gradient. Without this arm a gain could be the model family rather than
                   the channel, and item 5 would be credited with something it did not do

and the trunk alone as the reference, so a null result can be read as "nothing to add" rather
than "the member is broken".

Clipping is not optional and is not a thumb on the scale. Item 79 records that seed 0 contains a
compound the trunk predicts at -360, and `src/trunkdose.py` clips to the enzyme's label range
plus or minus two units before any post-processing for that reason. The same clip is used here,
so the comparison is the one the document already makes.

Reads oof_pool_all.json, oof_gp.json, oof_weak.json, trunk_twohead.json. Writes nothing. ~4 min.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr, ttest_1samp
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
SEEDS = (0, 1, 2, 3)
MEMBERS = [("oof_pool_all", "независимо"), ("oof_pool_all", "пул"),
           ("oof_gp", "GP"), ("oof_weak", "гребневая")]
LAMS = ("3.0", "0.0")   # с каналом и без; вторая --- контроль на семейство


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    cache = {f: json.load(open(RES + f"preds/{f}.json"))["preds"]
             for f in {m[0] for m in MEMBERS}}
    T = json.load(open(RES + "preds/trunk_twohead.json"))
    T = T.get("preds", T)

    table = []
    for seed in SEEDS:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        TK = {L: np.asarray(T[f"twohead|{seed}|{L}"], float) for L in LAMS}
        r = {"seed": seed}
        for tag in ("четыре", "пять", "пять lam0", "ствол один", "ствол lam0"):
            pr, rk = [], []
            for e, c in enumerate(CYPS):
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                fi, u = fold[m], np.ones(m.sum()) / m.sum()
                # Тот же клип, что в src/trunkdose.py: диапазон меток фермента +-2.
                t3 = np.clip(TK["3.0"][m, e], y.min() - 2.0, y.max() + 2.0)
                t0 = np.clip(TK["0.0"][m, e], y.min() - 2.0, y.max() + 2.0)
                ps = [np.asarray(cache[fn][f"{seed}|{arm}|{c}"], float) for fn, arm in MEMBERS]
                p = {"четыре": np.mean(ps, axis=0),
                     "пять": np.mean(ps + [t3], axis=0),
                     "пять lam0": np.mean(ps + [t0], axis=0),
                     "ствол один": t3,
                     "ствол lam0": t0}[tag]
                q = fit_apply(p, lo, hi, fi, u)
                pr.append(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)))
                rk.append(float(spearmanr(y, p).statistic))
                r[f"{c}|{tag}"] = rk[-1]
            r[f"пара {tag}"] = float(np.mean(pr))
            r[f"ранг {tag}"] = float(np.mean(rk))
        table.append(r)
        print(f"  сид {seed}: четыре {r['пара четыре']:.4f}/{r['ранг четыре']:.4f}  "
              f"пять {r['пара пять']:.4f}/{r['ранг пять']:.4f}  "
              f"пять-lam0 {r['пара пять lam0']:.4f}/{r['ранг пять lam0']:.4f}", flush=True)

    df = pd.DataFrame(table)
    dp = df["пара пять"] - df["пара четыре"]
    dr = df["ранг пять"] - df["ранг четыре"]
    print(f"\n{'':14s} {'пара':>9s} {'ранг':>9s}")
    for tag in ("четыре", "пять", "пять lam0", "ствол один", "ствол lam0"):
        print(f"{tag:14s} {df[f'пара {tag}'].mean():9.4f} {df[f'ранг {tag}'].mean():9.4f}")
    tp = ttest_1samp(dp, 0.0); tr_ = ttest_1samp(dr, 0.0)
    print(f"\nдобавление ствола: пара {dp.mean():+.4f} (t={tp.statistic:+.2f}, p={tp.pvalue:.3f}, "
          f"знаков {int((dp < 0).sum())}/4)")
    print(f"                   ранг {dr.mean():+.4f} (t={tr_.statistic:+.2f}, p={tr_.pvalue:.3f}, "
          f"знаков {int((dr > 0).sum())}/4)")
    d0 = df["пара пять lam0"] - df["пара четыре"]
    tt0 = ttest_1samp(d0, 0.0)
    print(f"контроль без канала: пара {d0.mean():+.4f} (t={tt0.statistic:+.2f}, "
          f"p={tt0.pvalue:.3f}, знаков {int((d0 < 0).sum())}/4)")
    print(f"вклад самого канала (пять минус пять-lam0): "
          f"{(df['пара пять'] - df['пара пять lam0']).mean():+.4f}")
    print("\nПо ферментам, прирост ранга от пятого члена:")
    for c in CYPS:
        d = df[f"{c}|пять"] - df[f"{c}|четыре"]
        print(f"  {c}: {d.mean():+.4f}  (знаков {int((d > 0).sum())}/4)")

    print("""
Как читать. Порог 0.007 и знак на четырёх сидах из четырёх; ни то ни другое поодиночке не
хватает при n = 4.

Контрольная рука решает вопрос атрибуции, и без неё вывод был бы неверен. Ствол при lam = 0
--- та же архитектура, то же число параметров, те же начальные веса, но скрининговая голова не
получает градиента. Если он даёт ансамблю столько же, то помогает НЕ канал, а просто другое
семейство моделей, и пункт 5 к выигрышу отношения не имеет. Если меньше --- разница между двумя
руками и есть вклад канала, измеренный внутри ансамбля.""")


if __name__ == "__main__":
    main()
