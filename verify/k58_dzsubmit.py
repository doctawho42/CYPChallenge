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

import json

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
    ap.add_argument("--cache", default="", help="куда сложить/откуда взять предсказания членов")
    ap.add_argument("--trunk-dead", default="",
                    help="файл ствола, обученного против спроецированной мишени "
                         "(src/trunk.py --dead). Даёт третью руку: проход во ВСЕХ пяти.")
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

    def load_cache():
        pp = a.cache + ".parts"
        if a.cache and not _pl.Path(a.cache).exists() and _pl.Path(pp).exists():
            C = json.load(open(pp))
            if C.get("seed") == a.seed and C.get("mode") == a.mode:
                print(f"члены взяты из {pp}, проход пересчитывается", flush=True)
                return ([(k, [np.asarray(v, float) for v in P]) for k, P in C["parts"]], None)
        if not a.cache or not _pl.Path(a.cache).exists():
            return None
        C = json.load(open(a.cache))
        if C.get("seed") != a.seed or C.get("mode") != a.mode:
            return None
        return ([(k, [np.asarray(v, float) for v in P]) for k, P in C["parts"]],
                [(k, [np.asarray(v, float) for v in P]) for k, P in C["dzp"]])

    got = load_cache()
    if got is not None and got[1] is not None:
        parts, dzp = got
        print(f"члены и проход взяты из {a.cache}", flush=True)
    elif got is not None:
        parts = got[0]
        t0 = time.time()
        dzp = SB.dz_pass(parts, X, mask, fold, (LO, HI))
        print(f"проход мёртвой зоны за {time.time()-t0:.0f} с", flush=True)
        if a.cache:
            json.dump({"seed": a.seed, "mode": a.mode,
                       "parts": [(k, [v.tolist() for v in P]) for k, P in parts],
                       "dzp": [(k, [v.tolist() for v in P]) for k, P in dzp]},
                      open(a.cache, "w"))
    else:
        t0 = time.time()
        parts = SB.oof_members(X, y, mask, fold, a.mode, seed=a.seed)
        print(f"члены посчитаны за {time.time()-t0:.0f} с: "
              + ", ".join(k for k, _ in parts), flush=True)
        if a.cache:
            # Члены складываются ДО прохода: он стоит вдвое дороже их, и падение в нём
            # не должно уносить с собой сорок минут уже посчитанного. Ровно это и
            # случилось, когда сторож ствола остановил сиды 1-3 после стадии членов.
            json.dump({"seed": a.seed, "mode": a.mode, "частично": True,
                       "parts": [(k, [v.tolist() for v in P]) for k, P in parts],
                       "dzp": []}, open(a.cache + ".parts", "w"))
        t0 = time.time()
        dzp = SB.dz_pass(parts, X, mask, fold, (LO, HI))
        print(f"проход мёртвой зоны за {time.time()-t0:.0f} с", flush=True)
        if a.cache:
            json.dump({"seed": a.seed, "mode": a.mode,
                       "parts": [(k, [v.tolist() for v in P]) for k, P in parts],
                       "dzp": [(k, [v.tolist() for v in P]) for k, P in dzp]},
                      open(a.cache, "w"))
            print(f"сложено в {a.cache}", flush=True)

    plain = [np.mean([P[e] for _, P in parts], axis=0) for e in range(len(CYPS))]
    dz = [np.mean([P[e] for _, P in dzp], axis=0) for e in range(len(CYPS))]

    arms = [("без прохода", plain), ("проход в четырёх", dz)]

    if a.trunk_dead:
        # Ствол, обученный против спроецированной мишени, подставляется НА МЕСТО обычного.
        # Обрезка применяется, потому что её применяет src/submit.py, и потому что проход
        # выброс не убрал, а переселил: на CYP2D6 минимум -360.26 стал максимумом 76.82.
        J = json.load(open(a.trunk_dead))["preds"]
        key = f"{SB.TRUNK_MODE}|{a.seed}|{SB.TRUNK_LAM}"
        if key not in J:
            raise SystemExit(f"нет ключа {key} в {a.trunk_dead}")
        A = np.asarray(J[key], float)
        td = [SB._trunk_clip(A[mask[:, e], e], y[mask[:, e], e]) for e in range(len(CYPS))]
        five = []
        for e in range(len(CYPS)):
            other = [P[e] for k, P in dzp if k != "ствол"]
            five.append(np.mean(other + [td[e]], axis=0))
        arms.append(("проход во всех пяти", five))

    rows_out = []
    for nm, P in arms:
        r = score(P, y, LO, HI, mask, fold)
        r["рука"] = nm
        rows_out.append(r)

    df = pd.DataFrame(rows_out).set_index("рука")
    cols = ["MACRO пара", "MACRO rho"] + [f"{c} rho" for c in CYPS]
    print()
    print(df[cols].round(4).to_string())
    base = df.loc["без прохода"]
    for nm in df.index:
        if nm == "без прохода":
            continue
        print(f"\n{nm}: прирост ранга {df.loc[nm, 'MACRO rho'] - base['MACRO rho']:+.4f}, "
              f"пары {df.loc[nm, 'MACRO пара'] - base['MACRO пара']:+.4f}")
    print(f"пункт 164 на четырёх сидах: ранг +0.0167, пара -0.0171 "
          f"(0.6063 -> 0.6230, 0.6824 -> 0.6653)")
    print("""
Как читать. Совпадения до четвёртого знака здесь быть не должно: пункт 164 складывал
СОХРАНЁННЫЕ предсказания из ablate/ablpool/ablgp, а этот файл пересчитывает члены кодом
подачи, и это два разных прогона одного оценщика. Ожидается тот же знак и порядок величины.
Знак наоборот или ноль --- дефект прохода, а не шум. Один сид, разброс макро по сидам 0.016.""")


if __name__ == "__main__":
    main()
