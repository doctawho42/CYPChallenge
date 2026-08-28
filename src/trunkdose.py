"""The lambda_scr dose curve, what the damage is made of, and at which delta any of it holds.

Four questions that a macro number at one lambda cannot answer, and one that says whether
the whole trunk-versus-boosting comparison was made at the right operating point.

1. Is the control exact? The two arms differ only inside `if lam > 0`, and the weight seed
   is hash((seed, fold)) with neither lambda nor mode in it, so at lambda = 0 the arms
   should be identical bit for bit, not merely close. That is a stronger statement than
   "the rows agree", and it is checkable, so check it rather than assert it.

2. What is the damage made of? Two stories predict opposite orderings between two enzymes,
   which is what makes them separable at all:
     saturation      - the fixed map g is blind where its derivative vanishes, so damage
                       follows the blind fraction, 3A4 .44 > 1A2 .32 > 2D6 .21 > 2C9 .09;
     bias absorption - one latent has nowhere to put a misfit of the population (E, h), so
                       pi absorbs it, and damage follows how badly that fit misses.
   2D6 and 3A4 rank opposite ways under the two, so the per-enzyme loss at lambda = 3
   decides it. The lambda slope corroborates: under saturation it should be steepest on
   3A4 and flattest on 2C9.

3. Is the post-isotonic gap between the arms noise? Measured against the RAW spread across
   seeds it looks like a third of it. But both arms have been recalibrated by then, and
   recalibration collapses the spread; the gap has to be measured against the spread of
   the thing it is a gap in.

4. Does any of it survive the tilt? Every trunk number in the document is computed at
   delta = 0, and src/reweight.py puts the test's label marginal at delta +0.1 to +0.6.
   The boosting reference moves under that tilt, so a comparison made at delta = 0 is a
   comparison at an operating point we have argued is the wrong one. There is a
   substantive expectation: the trunk's advantage is a rank advantage, while shrinkage's
   advantage over raw boosting is purely one of scale, and scale advantages are the first
   thing a shift of the mean takes away.

Reads results/preds/trunk_twohead.json and trunk_calibrated.json from src/trunk.py, and
results/preds/oof.json for the boosting reference. Writes nothing.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
# tilt() and wstrae() are definitions two scripts have to agree on, exactly as the split
# is, so they are imported rather than pasted. reweight.py guards its main().
from reweight import tilt, wstrae

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
LAMS = [0.0, 0.3, 1.0, 3.0]
OFFGRID = np.round(np.arange(-0.2, 1.61, 0.05), 2)
LAMGRID = np.round(np.arange(0.20, 1.01, 0.02), 2)

# Candidate per-enzyme orderings, all measured elsewhere, quoted here only to be correlated
# against the damage. Each is a hypothesis about what the damage is made of.
#
# BLIND_G is the honest saturation proxy: the share of compounds where |dg/dpi| falls below
# a tenth of its maximum. g is the map the loss actually compares against the screen, so its
# derivative is what governs how much the residual can say about pi.
#
# BLIND_I is the same fraction computed on |dI/dpi|, the derivative of the inhibited
# fraction rather than of log2(1 - I). It is the version quoted as 0.32/0.09/0.21/0.44 in
# earlier drafts. Kept because it is what the document says, and because the two order the
# enzymes identically - which is why the wrong derivative never showed up.
#
# EMAX_SPREAD was proposed as the competing "bias absorption" proxy, on the grounds that it
# ranks 2D6 above 3A4 and so points the opposite way from saturation. It does - but it is
# q95-q5 of the MEASURED per-compound Emax column in data/cyp-challenge-TRAIN_Emax.csv, a
# different quantity from the calibration E in a different file on a different sign
# convention. The spread of the actual fitted E is 30-60x smaller.
#
# HILL_RESID is the proxy that does measure misfit of a single fixed (E, h): the residual sd
# printed by verify/g1_calib.py. Note it correlates with EMAX_SPREAD at Spearman -1.000 -
# the two "misfit" proxies are perfectly reversed, and only one of them can be the right
# one. Under HILL_RESID, bias absorption predicts the same 3A4-first ordering as saturation,
# so the per-enzyme ordering stops separating the two stories at all.
BLIND_G = {"CYP1A2": 0.218, "CYP2C9": 0.068, "CYP2D6": 0.137, "CYP3A4": 0.397}
BLIND_I = {"CYP1A2": 0.316, "CYP2C9": 0.089, "CYP2D6": 0.212, "CYP3A4": 0.440}
EMAX_SPREAD = {"CYP1A2": 0.159, "CYP2C9": 0.240, "CYP2D6": 0.091, "CYP3A4": 0.067}
HILL_RESID = {"CYP1A2": 0.231, "CYP2C9": 0.185, "CYP2D6": 0.497, "CYP3A4": 0.612}
RHO_SCR = {"CYP1A2": 0.862, "CYP2C9": 0.896, "CYP2D6": 0.828, "CYP3A4": 0.936}
HYPS = [("слепая доля |dg/dpi| (насыщение)", BLIND_G),
        ("слепая доля |dI/dpi| (как в тексте)", BLIND_I),
        ("разброс измеренного Emax", EMAX_SPREAD),
        ("остаточное ско Хилла (промах E,h)", HILL_RESID),
        ("|rho| скрининг~pIC50", RHO_SCR)]
BLIND = BLIND_G


def truth():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    lo = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    hi = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    return list(rows.SMILES), y, lo, hi


def iso_oof(p, y, fold):
    """Out-of-fold isotonic: each fold gets a map fitted on the other four.

    Monotone, so within one fold's map the ranking is preserved exactly apart from ties,
    which isotonic does create by collapsing ranges to constants. Across the glued vector
    five different maps are in play, so rank is preserved only approximately - that is
    measured in block 5 rather than waved at.
    """
    o = np.full_like(p, np.nan)
    for f in np.unique(fold):
        te, trn = fold == f, fold != f
        if te.sum() == 0 or trn.sum() < 10:
            continue
        ir = IsotonicRegression(out_of_bounds="clip").fit(p[trn], y[trn])
        o[te] = ir.predict(p[te])
    return o


def fit_affine_oof(p, y, lo, hi, fold):
    """Offset and lambda fitted on the training folds, applied to the held-out fold.

    The family c + L(p - c) is identically a + b p with b = L and a = c(1 - L), so fixing
    the offset slices a two-parameter family along an arbitrary line. Both are searched.
    """
    o = np.full_like(p, np.nan)
    for f in np.unique(fold):
        te, trn = fold == f, fold != f
        if te.sum() == 0 or trn.sum() < 10:
            continue
        mu = p[trn].mean()
        best, bv = None, np.inf
        for off in OFFGRID:
            c = mu + off
            for L in LAMGRID:
                v = strae(y[trn], c + L * (p[trn] - c),
                          y_true_upper=hi[trn], y_true_lower=lo[trn])
                if v < bv:
                    bv, best = v, (c, L)
        c, L = best
        o[te] = c + L * (p[te] - c)
    return o


def per_enzyme(y, lo, hi, P):
    """ST-RAE per enzyme over each enzyme's own labelled cells."""
    out = []
    for e in range(4):
        m = ~np.isnan(y[:, e]) & ~np.isnan(P[:, e])
        out.append(float(strae(y[m, e], P[m, e], y_true_upper=hi[m, e], y_true_lower=lo[m, e])))
    return out


def main():
    smiles, y, lo, hi = truth()
    blob = {}
    for mode in ("twohead", "calibrated"):
        f = _pl.Path(RES + f"preds/trunk_{mode}.json")
        if not f.exists():
            raise SystemExit(f"нет {f} - сначала: uv run python src/trunk.py --mode {mode} "
                             "--seeds 0,1,2,3 --lams 0,0.3,1.0,3.0")
        blob[mode] = json.load(open(f))
    P = {k: np.asarray(v, float) for m in blob for k, v in blob[m]["preds"].items()}
    seeds = sorted({int(k.split("|")[1]) for k in P})
    folds = {s: butina_folds(smiles, seed=s)[0] for s in seeds}

    print("=" * 92)
    print("1. Точен ли контроль: совпадают ли руки при lambda = 0")
    print("=" * 92)
    print("Ветка режима стоит внутри `if lam > 0`, а зерно весов - hash((seed, fold)) без")
    print("lambda и без режима. Значит при lambda = 0 руки обязаны совпасть побитово.\n")
    worst = 0.0
    for s in seeds:
        a, b = P[f"twohead|{s}|0.0"], P[f"calibrated|{s}|0.0"]
        m = ~np.isnan(a) & ~np.isnan(b)
        d = float(np.max(np.abs(a[m] - b[m])))
        worst = max(worst, d)
        print(f"  зерно {s}: макс |разность| = {d:.3e}   "
              f"{'побитово' if d == 0 else 'РАСХОЖДЕНИЕ'}")
    print("\n" + ("Контроль точен по построению, а не по совпадению." if worst == 0
                   else "ВНИМАНИЕ: руки при lambda=0 разошлись, сравнение не контролируемо."))

    print()
    print("=" * 92)
    print("2. Дозовая кривая: макро ST-RAE, среднее по зёрнам (ниже - лучше)")
    print("=" * 92)
    curve = {}
    for mode in ("twohead", "calibrated"):
        row = {}
        for lam in LAMS:
            v = [np.mean(per_enzyme(y, lo, hi, P[f"{mode}|{s}|{lam}"])) for s in seeds]
            row[lam] = (float(np.mean(v)), float(np.std(v, ddof=1)))
        curve[mode] = row
    print(f"{'lambda':>7s} | {'двухголовая':>22s} | {'калиброванная':>22s} | {'во сколько раз':>14s}")
    print(f"{'':>7s} | {'макро':>9s} {'потеря':>7s} {'sd':>4s} | "
          f"{'макро':>9s} {'потеря':>7s} {'sd':>4s} | {'хуже':>14s}")
    b0, c0 = curve["twohead"][0.0][0], curve["calibrated"][0.0][0]
    for lam in LAMS:
        bm, bs = curve["twohead"][lam]
        cm, cs = curve["calibrated"][lam]
        db, dc = bm - b0, cm - c0
        r = f"{dc/db:13.1f}x" if abs(db) > 1e-9 else f"{'-':>14s}"
        print(f"{lam:7.1f} | {bm:9.4f} {db:+7.4f} {bs:4.3f} | "
              f"{cm:9.4f} {dc:+7.4f} {cs:4.3f} | {r}")

    print()
    print("=" * 92)
    print("3. Из чего сделана порча: потеря по ферментам при lambda = 3 против lambda = 0")
    print("=" * 92)
    dmg = {}
    for mode in ("twohead", "calibrated"):
        base = np.mean([per_enzyme(y, lo, hi, P[f"{mode}|{s}|0.0"]) for s in seeds], 0)
        top = np.mean([per_enzyme(y, lo, hi, P[f"{mode}|{s}|3.0"]) for s in seeds], 0)
        dmg[mode] = top - base
    print(f"{'фермент':8s} {'слепая доля':>12s} {'двухгол.':>10s} {'калибр.':>9s} {'отношение':>10s}")
    for e, c in enumerate(CYPS):
        r = dmg["calibrated"][e] / dmg["twohead"][e] if abs(dmg["twohead"][e]) > 1e-9 else np.nan
        print(f"{c:8s} {BLIND[c]:12.2f} {dmg['twohead'][e]:+10.4f} "
              f"{dmg['calibrated'][e]:+9.4f} {r:10.1f}")
    order = [CYPS[i] for i in np.argsort(-dmg["calibrated"])]
    print(f"\nпорядок порчи калиброванной руки, от худшего: {' > '.join(order)}")
    print(f"\n{'гипотеза':38s} {'порядок ферментов':34s} {'Спирмен с порчей':>17s}")
    for name, h in HYPS:
        o = " > ".join(c[3:] for c in sorted(CYPS, key=lambda c: -h[c]))
        r = spearmanr([h[c] for c in CYPS], dmg["calibrated"]).statistic
        print(f"{name:38s} {o:34s} {r:+17.3f}")
    print("\nn = 4: перестановочная вероятность |Спирмена| = 1 равна 1/12 = 0.083 для")
    print("двустороннего чтения. Ни одна строка сама по себе не доказывает механизм; смысл")
    print("таблицы в том, какие гипотезы порядок ИСКЛЮЧАЕТ.")

    print()
    print("=" * 92)
    print("4. Крутизна по lambda, по ферментам (насыщение: круче всего 3A4, площе всего 2C9)")
    print("=" * 92)
    print(f"{'фермент':8s} " + " ".join(f"{('l=' + str(l)):>9s}" for l in LAMS) + f" {'наклон':>9s}")
    for e, c in enumerate(CYPS):
        v = [float(np.mean([per_enzyme(y, lo, hi, P[f"calibrated|{s}|{l}"])[e] for s in seeds]))
             for l in LAMS]
        sl = float(np.polyfit(LAMS, v, 1)[0])
        print(f"{c:8s} " + " ".join(f"{x:9.4f}" for x in v) + f" {sl:+9.4f}")

    print()
    print("=" * 92)
    print("5. Разрыв после изотоники: против какой линейки его мерить")
    print("=" * 92)
    iso = {}
    for mode in ("twohead", "calibrated"):
        for lam in (0.0, 3.0):
            for s in seeds:
                p = P[f"{mode}|{s}|{lam}"]
                q = np.full_like(p, np.nan)
                for e in range(4):
                    m = ~np.isnan(y[:, e]) & ~np.isnan(p[:, e])
                    q[m, e] = iso_oof(p[m, e], y[m, e], folds[s][m])
                iso[(mode, lam, s)] = float(np.mean(per_enzyme(y, lo, hi, q)))
    print(f"{'рука':28s} " + " ".join(f"{('зерно ' + str(s)):>9s}" for s in seeds)
          + f" {'среднее':>9s} {'sd':>7s} {'размах':>7s}")
    lab = {("twohead", 0.0): "двухголовая l=0 + изо",
           ("twohead", 3.0): "двухголовая l=3 + изо",
           ("calibrated", 3.0): "калиброванная l=3 + изо"}
    for k, name in lab.items():
        v = [iso[(k[0], k[1], s)] for s in seeds]
        print(f"{name:28s} " + " ".join(f"{x:9.4f}" for x in v)
              + f" {np.mean(v):9.4f} {np.std(v, ddof=1):7.4f} {max(v)-min(v):7.4f}")
    gap = [iso[("calibrated", 3.0, s)] - iso[("twohead", 3.0, s)] for s in seeds]
    sd_post = float(np.std([iso[("twohead", 3.0, s)] for s in seeds], ddof=1))
    rng_post = float(np.ptp([iso[("twohead", 3.0, s)] for s in seeds]))
    sd_raw = curve["twohead"][3.0][1]
    print(f"\nразрыв калиброванная - двухголовая при l=3, по зёрнам: "
          + " ".join(f"{g:+.4f}" for g in gap))
    print(f"знак одинаков на всех {len(seeds)} зёрнах: "
          f"{bool(np.all(np.sign(gap) == np.sign(gap[0])))}, среднее {np.mean(gap):+.4f}")
    print(f"против СЫРОГО разброса   sd {sd_raw:.4f}: {abs(np.mean(gap))/sd_raw:.2f} sd")
    print(f"против разброса ПОСЛЕ изо sd {sd_post:.4f}, размах {rng_post:.4f}: "
          f"{abs(np.mean(gap))/sd_post:.2f} sd, {abs(np.mean(gap))/rng_post:.2f} размаха")
    print("Обе руки к этому моменту перекалиброваны, поэтому линейка - вторая, а не первая.")

    print()
    print("=" * 92)
    print("6. Под наклоном: при каком delta сравнение ствола с бустингом ещё то же самое")
    print("=" * 92)
    oof = json.load(open(RES + "preds/oof.json"))
    f0 = folds[0]
    arms = {}
    for e, c in enumerate(CYPS):
        m = ~np.isnan(y[:, e])
        g = np.asarray(oof[f"FP+DESC+MECH|{c}"], float)
        arms.setdefault("бустинг сырой", {})[c] = g
        arms.setdefault("бустинг + аффинная пара", {})[c] = fit_affine_oof(
            g, y[m, e], lo[m, e], hi[m, e], f0[m])
        for mode, lam, name in (("twohead", 0.0, "ствол l=0 + изо"),
                                ("twohead", 3.0, "ствол двухгол. l=3 + изо"),
                                ("calibrated", 3.0, "ствол калибр. l=3 + изо")):
            p = P[f"{mode}|0|{lam}"][m, e]
            arms.setdefault(name, {})[c] = iso_oof(p, y[m, e], f0[m])
        # Like for like. Above, the boosting gets the affine pair and the trunk gets
        # isotonic - two different post-processings, so the comparison confounds model with
        # post-processing. These two rows give the trunk exactly what the boosting got.
        for mode, lam, name in (("twohead", 3.0, "ствол двухгол. l=3 + пара"),
                                ("calibrated", 3.0, "ствол калибр. l=3 + пара")):
            p = P[f"{mode}|0|{lam}"][m, e]
            arms.setdefault(name, {})[c] = fit_affine_oof(
                p, y[m, e], lo[m, e], hi[m, e], f0[m])
    deltas = np.round(np.arange(0.0, 0.61, 0.1), 1)
    print(f"{'вариант':26s} " + " ".join(f"{('d=' + str(d)):>8s}" for d in deltas))
    tabs = {}
    for name, byc in arms.items():
        row = []
        for d in deltas:
            v = []
            for e, c in enumerate(CYPS):
                m = ~np.isnan(y[:, e])
                w = tilt(y[m, e], float(d))
                v.append(wstrae(y[m, e], byc[c], lo[m, e], hi[m, e], w))
            row.append(float(np.mean(v)))
        tabs[name] = row
        print(f"{name:26s} " + " ".join(f"{x:8.4f}" for x in row))
    print()
    a, b = "бустинг + аффинная пара", "ствол двухгол. l=3 + изо"
    print(f"{b + ' минус ' + a:26.26s} " + " ".join(
        f"{tabs[b][i] - tabs[a][i]:+8.4f}" for i in range(len(deltas))))
    flip = [i for i in range(len(deltas)) if (tabs[b][i] - tabs[a][i]) < 0]
    if flip:
        print(f"\nПорядок переворачивается начиная с delta = {deltas[flip[0]]}: до него лучше")
        print("бустинг со сдвинутым центром, после - ствол. Правдоподобный диапазон delta")
        print("по src/reweight.py - от +0.1 до +0.6, так что точка переворота лежит внутри")
        print("него, а не за ним.")
    else:
        print("\nПорядок не переворачивается ни при одном delta из проверенного диапазона.")
    print("""
Оговорка, без которой шестой блок читается сильнее, чем следует. Наклон правит маргиналь
МЕТОК и предполагает, что p(y | yhat) на тесте та же; постобработка здесь подогнана при
delta = 0 и под наклоном только оценивается, что и есть вопрос переноса, а не оракула.
Геометрию сходства наклон не трогает вовсе.""")


if __name__ == "__main__":
    main()
