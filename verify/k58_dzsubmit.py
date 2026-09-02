"""Does the dead-zone pass inside src/submit.py reproduce item 164's +0.0167?

Item 164 measured the pass by composing SAVED predictions -- `oof_pool_all`, `oof_gp`,
`oof_weak` -- in `verify/k46_five.py`, and got +0.0167 of rank on the submitted five-member
configuration over four seeds. `src/submit.py` does not read those files; it recomputes every
member from `data/feats.npz`. So the pass now living in `submit.py` is a second implementation
of a measured thing, and the only honest way to trust it is to reproduce the number.

That is the whole point of this file. Item 198 is the reason it exists: there I wrote a
booster from scratch instead of reusing the debugged one next door, reproduced exactly the
defect its comment warned about, and the resulting table looked entirely plausible for a day.
A new implementation of a measured quantity is not trustworthy until it reproduces it.

What is compared, on the same folds and the same seed:

    ансамбль5, без прохода      must land near item 120's 0.6063 of rank
    ансамбль5, проход в четырёх must land near item 164's 0.6230, i.e. about +0.0167

The trunk is not reprojected in either arm, exactly as in item 164, so the comparison is the
one that file made rather than a different one that happens to be nearby.

Two things this check cannot do. It runs one seed, and the per-seed spread of macro rank is
0.016 (f3), so a single seed cannot confirm +0.0167 -- it can only show that the pass fires
and moves rank the right way by roughly the right amount. And it recomputes members rather
than reading them, so exact agreement with item 164 is not expected; a difference of a few
thousandths is the noise between two implementations of the same estimator, a difference of
0.01 or a sign flip is a defect.

Reads data/feats.npz, data/rows.csv, results/preds/trunk_twohead.json. Writes nothing.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

import submit as SB
from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = SB.CYPS


def score(P, y, LO, HI, mask, fold):
    r = {}
    for e, c in enumerate(CYPS):
        m = mask[:, e]
        yy, lo, hi, fi = y[m, e], LO[m, e], HI[m, e], fold[m]
        p = P[e]
        q = fit_apply(p, lo, hi, fi, np.ones(len(yy)) / len(yy))
        r[f"{c} пара"] = float(strae(yy, q, y_true_upper=hi, y_true_lower=lo))
        r[f"{c} rho"] = float(spearmanr(yy, p).statistic)
    for t in ("пара", "rho"):
        r[f"MACRO {t}"] = float(np.mean([r[f"{c} {t}"] for c in CYPS]))
    return r


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--mode", default="ансамбль5")
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)
    fold, _ = butina_folds(list(rows.SMILES), seed=a.seed)

    print(f"сид {a.seed}, режим {a.mode}\n", flush=True)
    t0 = time.time()
    parts = SB.oof_members(X, y, mask, fold, a.mode)
    print(f"члены посчитаны за {time.time()-t0:.0f} с: "
          + ", ".join(k for k, _ in parts), flush=True)

    plain = [np.mean([P[e] for _, P in parts], axis=0) for e in range(len(CYPS))]
    t0 = time.time()
    dzp = SB.dz_pass(parts, X, mask, fold, (LO, HI))
    print(f"проход мёртвой зоны за {time.time()-t0:.0f} с", flush=True)
    dz = [np.mean([P[e] for _, P in dzp], axis=0) for e in range(len(CYPS))]

    rows_out = []
    for nm, P in (("без прохода", plain), ("проход в четырёх", dz)):
        r = score(P, y, LO, HI, mask, fold)
        r["рука"] = nm
        rows_out.append(r)

    df = pd.DataFrame(rows_out).set_index("рука")
    cols = ["MACRO пара", "MACRO rho"] + [f"{c} rho" for c in CYPS]
    print()
    print(df[cols].round(4).to_string())
    d_rho = df.loc["проход в четырёх", "MACRO rho"] - df.loc["без прохода", "MACRO rho"]
    d_pair = df.loc["проход в четырёх", "MACRO пара"] - df.loc["без прохода", "MACRO пара"]
    print(f"\nприрост ранга {d_rho:+.4f}, пары {d_pair:+.4f}")
    print(f"пункт 164 на четырёх сидах: ранг +0.0167, пара -0.0171 "
          f"(0.6063 -> 0.6230, 0.6824 -> 0.6653)")
    print("""
Как читать. Совпадения до четвёртого знака здесь быть не должно: пункт 164 складывал
СОХРАНЁННЫЕ предсказания из ablate/ablpool/ablgp, а этот файл пересчитывает члены кодом
подачи, и это два разных прогона одного оценщика. Ожидается тот же знак и порядок величины.
Знак наоборот или ноль --- дефект прохода, а не шум. Один сид, разброс макро по сидам 0.016.""")


if __name__ == "__main__":
    main()
