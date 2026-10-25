"""Keep one submission's ORDERING and another's SCALE: map the source's per-enzyme ranks onto the
reference's sorted values.

Why this exists. On 25 October the board scored `activity_submission_stack374.csv` at MA-ST-RAE
0.7715 and rank 135 with **MA-Spearman 0.7212, the best this project has ever had** --- against
`best4`'s 0.5357 and 0.7062. The ordering improved by +0.0150, more than the +0.0099 measured out of
fold, and the scale collapsed: MA-R² fell from 0.4946 to 0.1454. That is item 308's diagnosis word
for word, "the ordering is better than we predicted and the scale is far worse", and item 366's
mechanism --- ST-RAE scores BAND CONTAINMENT against credible intervals whose median half-widths are
0.16 to 0.26, so where the level sits matters more than how good the model is.

The levels say the same thing. `best4` predicts CYP2D6 at a mean of 3.725, a full log unit BELOW the
training mean of 4.784, and scores well; `stack374` puts it at 4.659 and does not. The challenge's
own About page explains why: the test set is an analogue expansion of the top 25 hits of CYP1A2,
CYP2C9 and CYP3A4 only, so for any one isoform two thirds of the test were chosen for a DIFFERENT
isoform and are largely inactive on it. The test distribution sits below the training set's, which
was built from compounds promoted to dose-response BECAUSE they were hits.

A rank map is a monotone per-enzyme transform of our own predictions onto our own earlier
predictions' empirical distribution. No test label is involved and none may be. It is the same
family as the `--delta`, `--spread` and affine-pair calibration this repository already ships, with
more freedom: it reproduces the reference's distribution exactly rather than matching two moments.

    uv run python src/rankmap.py --source results/submission/activity_submission_stack374.csv \
                                 --reference results/submission/activity_submission_best4.csv \
                                 --out results/submission/activity_submission_rankmap.csv

FOUR CONTROLS, all of which must succeed, and two of them are the point of the file:

  * per-enzyme Spearman between output and SOURCE is exactly 1.000000 --- the ordering, which is
    the thing worth keeping, survives untouched;
  * the output's sorted values equal the REFERENCE's sorted values bit for bit --- the scale, which
    is the thing that scored, is reproduced and not approximated;
  * molecule names match the test table element for element, and the same comparison on the names
    rolled by one must FAIL;
  * the output differs from BOTH inputs, or the map did nothing and the run is a no-op dressed up
    as a result.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from submit import CYPS

MIN_STD = 0.01   # the Space's own floor: config.py MIN_PREDICTION_STD


def rank_map(src, ref):
    """src's order, ref's values. Ties in src are broken by the stable argsort, so the map is
    deterministic; `rank_map` is the identity when src and ref are the same vector."""
    if len(src) != len(ref):
        raise SystemExit("длины не совпали: %d против %d" % (len(src), len(ref)))
    out = np.empty(len(src), float)
    out[np.argsort(src, kind="stable")] = np.sort(ref)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--source", required=True, help="файл, чей ПОРЯДОК берём")
    ap.add_argument("--reference", required=True, help="файл, чью ШКАЛУ берём")
    ap.add_argument("--out", required=True)
    a = ap.parse_args(argv)

    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    S = pd.read_csv(a.source)
    R = pd.read_csv(a.reference)
    for nm, df in (("источник", S), ("эталон", R)):
        if len(df) != len(te):
            raise SystemExit("%s: %d строк против %d в тестовой таблице" % (nm, len(df), len(te)))
        n = df.Molecule_Name.to_numpy(object)
        t = te.Molecule_Name.to_numpy(object)
        ok, rolled = bool((n == t).all()), bool((np.roll(n, 1) == t).all())
        if not ok or rolled:
            raise SystemExit("%s: порядок молекул (совпал=%s, сдвинутый=%s)" % (nm, ok, rolled))
        print("  %-9s %s, молекулы сверены, сдвинутый контроль False" % (nm, df.shape), flush=True)

    out = S.copy()
    print("\n  %-8s %10s %10s %10s %10s %12s" % ("фермент", "ср источ.", "ср эталон",
                                                 "sd источ.", "sd эталон", "ро с источ."))
    for c in CYPS:
        col = "%s_pIC50_direct_inhibition" % c
        s = S[col].to_numpy(float)
        r = R[col].to_numpy(float)
        if not (np.isfinite(s).all() and np.isfinite(r).all()):
            raise SystemExit("%s: не-конечные значения на входе" % c)
        v = rank_map(s, r)
        rho = float(spearmanr(s, v).statistic)
        if abs(rho - 1.0) > 1e-12:
            raise SystemExit("%s: порядок НЕ сохранён, ро = %.9f" % (c, rho))
        d = float(np.max(np.abs(np.sort(v) - np.sort(r))))
        if d != 0.0:
            raise SystemExit("%s: шкала не воспроизведена, max|Δ| отсортированных = %.2e" % (c, d))
        if v.std(ddof=1) < MIN_STD:
            raise SystemExit("%s: разброс %.4f ниже порога площадки %.2f" % (c, v.std(ddof=1),
                                                                             MIN_STD))
        out[col] = v
        print("  %-8s %10.3f %10.3f %10.3f %10.3f %12.6f"
              % (c, s.mean(), r.mean(), s.std(ddof=1), r.std(ddof=1), rho))
    print("\n  КОНТРОЛЬ порядка: ро с источником 1.000000 на всех четырёх")
    print("  КОНТРОЛЬ шкалы: отсортированные значения совпали с эталоном побитово")

    cols = ["%s_pIC50_direct_inhibition" % c for c in CYPS]
    for nm, df in (("источника", S), ("эталона", R)):
        if np.allclose(out[cols].to_numpy(float), df[cols].to_numpy(float)):
            raise SystemExit("результат совпал с %s — отображение ничего не сделало" % nm)
    print("  КОНТРОЛЬ непустоты: результат отличается и от источника, и от эталона")

    out.to_csv(a.out, index=False)
    print("\n  записано %s" % a.out)

    from validation.activity_validation import validate_activity_submission
    ok, errs = validate_activity_submission(a.out, expected_ids=set(te.Molecule_Name))
    print("  валидатор организаторов: %s" % ("принято" if ok else "ОТКЛОНЕНО %s" % errs))
    # The clause has to be shown live, or "accepted" proves nothing about the id check.
    bad, berrs = validate_activity_submission(a.out, expected_ids={"нет-такой"})
    print("  КОНТРОЛЬ валидатора (заведомо неверные id, ожидается отказ): %s"
          % ("ОТКАЗ как и ожидалось" if not bad else "ПРИНЯЛ — проверка не работает"))
    return 0 if ok and not bad else 2


if __name__ == "__main__":
    raise SystemExit(main())
