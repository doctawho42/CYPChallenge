"""How much ORDERING do these labels even permit? Converted into a ceiling on Spearman, the
published assay error gives 0.978 of macro against 0.7062 for us and 0.7971 for the best ordering
on the board --- and the deficit is LARGEST where the labels are cleanest.

WHY THIS QUESTION AND WHY NOW. The goal is first place; item 353 decomposed the gap and found it is
ordering rather than placement, 38.5 per cent of it on CYP2D6; item 357 then put placement's entire
in-sample headroom at 0.0386 of macro ST-RAE. Every lever since has been an attempt to buy ordering
and eight have been closed (items 358--361). Before buying more it is worth knowing whether the
ordering exists to be bought. This needs no GPU and nothing that is not already on disk.

WHAT IS ALREADY KNOWN AND MUST NOT BE RE-CLAIMED. **Item 221 reached the conclusion first, in the
error-magnitude metric:** with the assay floor at 0.04 of pIC50 against a model RMSE of 0.735, the
instrument explains 0.3 per cent of our error variance and 99.7 per cent is epistemic. `src/ceiling.py`
turns the same column into reliability --- 1 - E[sigma^2]/Var(y), the R^2 ceiling --- and gets 0.898
to 0.948. So "the labels are not the wall" is item 221's finding, not this one.

**What has never been computed here is the ceiling on SPEARMAN, and it is not that number.** Rank
correlation attenuates according to how the noise is DISTRIBUTED over the compounds, not according
to its mean variance, so it cannot be read off R^2 and has to be simulated. Spearman is also the
metric the gap actually lives in: item 322's within-entrant slope converts rho into ST-RAE at -1.62,
and item 353's decomposition is in rho. The ceiling computed here is E rho(latent truth, one
measured panel) over 200 redraws, with the latent sample being the label distribution rescaled to
the variance left after noise.

    фермент   надёжность   потолок A   повтор B   гомоск. A |   наш ро   лучший   наш %   лучш %
    CYP1A2       0.9481      0.9892     0.9792     0.9460   |   0.7513   0.8561     76 %    87 %
    CYP2C9       0.9050      0.9665     0.9357     0.9359   |   0.8231   0.8836     85 %    91 %
    CYP2D6       0.9335      0.9786     0.9607     0.9231   |   0.4359   0.5745     45 %    59 %
    CYP3A4       0.8978      0.9783     0.9592     0.9392   |   0.8146   0.8741     83 %    89 %
    МАКРО                    0.9782     0.9587              |   0.7062   0.7971     72 %    81 %

**The board sits at 81 per cent of what the labels permit and we sit at 72 per cent.** Column B is
the stricter reference --- what a SECOND run of the same assay would score against the first --- and
even that is 0.9587, far above anyone on the board. CYP2D6, carrying 38.5 per cent of our gap, is
the extreme case at 45 per cent of its own ceiling.

WHY THE RANK CEILING EXCEEDS THE VARIANCE ONE, and it is not an artefact. The `гомоск.` column holds
the mean noise variance fixed and makes it constant over compounds: 0.946 / 0.936 / 0.923 / 0.939,
clearly below the measured 0.989 / 0.967 / 0.979 / 0.978. The published error is strongly
heteroscedastic --- median sigma 0.07 to 0.14 against an rms of 0.23 to 0.35 --- so most labels are
far more precise than the variance decomposition suggests and a rank metric barely notices the few
that are not. Simulating rather than using the closed form is what makes that visible: the Gaussian
(6/pi)arcsin(r/2) under-reads the ceiling by 0.018 to 0.036.

WHAT `_std` IS, STATED BEFORE IT IS USED. Item 114 measured that the band width is a DETERMINISTIC
function of the label --- isotonic R^2 0.928 to 0.970, reproduced here as control 7 to three
decimals --- so this column is not per-compound measured uncertainty; it is a formula, and item 220
says so in those words. What item 220 did establish is its SCALE, by a chain closing entirely on the
organisers' side: well noise 0.161 log2fc propagated through a twelve-point fit lands on their own
published se, and our band is 3.8 times that, i.e. on the reproducibility rather than the precision
scale. For a ceiling on rank that is enough, because what the simulation needs is the noise scale as
a function of potency, which is exactly what a label-dependent width formula supplies. It is not
enough for any per-compound claim, and none is made.

SENSITIVITY, IN BOTH DIRECTIONS, since item 220 fixes which direction is live:

    kappa      0.27     1.00     1.25     1.50     2.00     3.00     4.00
    макро    0.9984   0.9782   0.9647   0.9478   0.9063   0.8128   0.7251

kappa = 1 is the published width, already the generous reading on item 220's finding; kappa = 0.27 is
the same band read on the precision scale; and item 221's label-INDEPENDENT instrument floor of 0.04
pIC50 gives 0.9979, which is the arm that holds nothing of item 114's formula at all. **The band's
scale would have to be understated by a factor of 3.2 for the best ordering on the board to already
be at the ceiling, and 4.2 for ours to be** --- on a quantity item 220 measured as 3.8 times too
WIDE, not too narrow. Section 4 pushes the same way: item 354 bounded the test's label spread below
at 1.346 / 0.976 / 1.359 / 1.116 against training's 1.030 / 0.782 / 0.916 / 1.093, and a wider truth
at the same noise is a higher ceiling, 0.9861 macro.

THE PART THAT IS NEW RATHER THAN A RESTATEMENT, and it is the actionable half. Item 114's formula
puts the width on the WEAK compounds: corr(sigma, y) is -0.74 to -0.90 and the weak half of each
enzyme carries 95 per cent of the whole noise variance. So the potent half is measured almost
exactly, and each half has to be scored against its own ceiling:

    фермент  | сильные: наше ро  потолок  доля | слабые: наше ро  потолок  доля
    CYP1A2   |          0.2969   0.9772   30 % |         0.4705   0.9507   49 %
    CYP2C9   |          0.4783   0.9564   50 % |         0.5407   0.8042   67 %
    CYP2D6   |          0.3660   0.9845   37 % |         0.2571   0.8594   30 %
    CYP3A4   |          0.5883   0.9878   60 % |         0.6734   0.8001   84 %
    среднее                               44 %                            58 %

**We take 44 per cent of the ceiling where the labels are clean and 58 per cent where they are
noisy.** The ordering deficit is concentrated in the half of the data with almost no measurement
error in it. This is NOT item 115's mistake in a new wrapper: item 115 measured where the ST-RAE
NUMERATOR sits and found it U-shaped, 45.5 per cent potent and 27.6 per cent weak, and nothing here
contradicts that. The claim is about ordering within a half against that half's own rank ceiling,
which is a different quantity and points the other way. (These are out-of-fold numbers --- macro rho
0.6434 against 0.7062 on the board, the test being an easier sample, item 357 --- so each share is a
lower bound on the board-side share. The comparison BETWEEN halves is unaffected, both sides being
out of fold.)

WHAT THIS DOES NOT SAY, and the distinction is the whole value of the item. It bounds what the
LABELS permit, not what STRUCTURE determines. This project has measured the second kind twice and
found it low --- items 134/135 on the screening gap and item 229 on the TDI shift, where MCC near
0.3 is close to what the label allows from structure. A ceiling of 0.978 says a model that recovered
the latent truth would score it; it does not promise the truth is recoverable from a SMILES. So the
reading is narrow and still worth having: **the remaining gap to first place is not a measurement
limit, and the claim that it is can now be refused with a number.** The +0.1019 of macro rho item
353 needs is 37 per cent of the +0.2719 the labels leave unclaimed, and three quarters of what is
unclaimed sits on CYP1A2 (+0.238) and CYP2D6 (+0.543) --- the two enzymes where the potent half is
also furthest below its own ceiling.

    uv run python verify/k138_noiseceiling.py

Reads the inhibition table, data/rows.csv, results/preds/oof_submitted.json and the logs of k128 and
k131. Writes results/logs/k138_noiseceiling.json. Runs in about a minute.

Controls, all seven able to fail and four of them against numbers committed elsewhere:
 1. (hi-lo)/`_std` must be 2*1.96 on all four enzymes, which is what licenses the band reading;
 2. the labels must carry item 354's per-enzyme training sd to 1e-9;
 3. at sigma == 0 both ceilings must be exactly 1.0 to 1e-12;
 4. known answer --- homoscedastic Gaussian data at reliability 0.80 must return the closed form
    (6/pi)arcsin(r/2) within 0.01 for both ceilings;
 5. the ceiling must fall strictly as kappa rises;
 6. the k128 log must carry item 353's committed board row;
 7. the label-dependence of the width must reproduce item 114's table within 0.01.
Exit: 0 ran, 2 a control failed, 3 an input is missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import json

import numpy as np
import pandas as pd
from scipy.stats import rankdata, spearmanr
from sklearn.isotonic import IsotonicRegression

from submit import CYPS

LOG = RES + "logs/k138_noiseceiling.json"
K128 = RES + "logs/k128_boardgap.json"
K131 = RES + "logs/k131_scaleshift.json"

REPEATS = 200           # ceilings
REPEATS_FAST = 60       # the kappa sweep and the bisection
Z95 = 1.959963984540054  # the band is 2*Z95*_std; control 1 checks that it is
KAPPA = (0.27, 1.0, 1.25, 1.5, 2.0, 3.0, 4.0)
FLOOR221 = 0.04      # item 221: the instrument noise floor in pIC50, label-independent
# Item 114's committed table: the band width is a function of the label, not a measurement.
ITEM114 = {"CYP1A2": (-0.885, 0.963), "CYP2C9": (-0.899, 0.928),
           "CYP2D6": (-0.558, 0.955), "CYP3A4": (-0.928, 0.970)}
N_BOARD = 375           # the public half the live board scores
# Item 354's committed per-enzyme training sd, which identifies the labels being read.
ITEM354_SD = {"CYP1A2": 1.0301112449331375, "CYP2C9": 0.7819900948636752,
              "CYP2D6": 0.9157893015336759, "CYP3A4": 1.0929851992409918}
# Item 353's committed board row, which identifies the k128 log being read.
ITEM353 = {"us_strae": 0.5357, "us_rho": 0.7062, "first_strae": 0.3706, "n_board": 273}
SLOPE_WITHIN = -1.62    # item 322: ST-RAE per unit rho WITHIN an entrant's own trajectory


def rho(a, b):
    """Spearman, as Pearson of average ranks -- same tie handling as scipy.spearmanr."""
    ra, rb = rankdata(a), rankdata(b)
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    return float(ra @ rb / np.sqrt((ra @ ra) * (rb @ rb)))


def ceilings(t, sig, Z1, Z2, kappa=1.0, take=None):
    """E rho(truth, one noisy panel) and E rho(two independent noisy panels).

    `t` is the latent sample, `sig` its per-compound assay sd, `Z1`/`Z2` pre-drawn standard
    normals of shape (repeats, n) so that a kappa sweep reuses the same draws and the curve
    is monotone up to nothing but kappa itself."""
    A, B = [], []
    for r in range(len(Z1)):
        e1 = kappa * sig * Z1[r]
        e2 = kappa * sig * Z2[r]
        k = slice(None) if take is None else take[r]
        A.append(rho(t[k], (t + e1)[k]))
        B.append(rho((t + e1)[k], (t + e2)[k]))
    return float(np.mean(A)), float(np.mean(B)), float(np.std(A))


def gauss_spearman(r):
    """Spearman of a bivariate normal with Pearson correlation r -- the closed-form reference."""
    return float(6.0 / np.pi * np.arcsin(r / 2.0))


def latent(y, mean_var, var_target=None):
    """A latent sample with the label distribution's SHAPE and the variance left after noise."""
    v = y.var(ddof=1) if var_target is None else var_target
    c = np.sqrt(max(v - mean_var, 0.0)) / y.std(ddof=1)
    return y.mean() + c * (y - y.mean()), max(0.0, 1.0 - mean_var / v)


def main():
    f = D + "cyp-challenge-TRAIN_inhibition.csv"
    if not _pl.Path(f).exists():
        print("нет %s --- make setup" % f, file=_sys.stderr)
        return 3
    for p in (K128, K131):
        if not _pl.Path(p).exists():
            print("нет %s --- сначала его скрипт" % p, file=_sys.stderr)
            return 3
    b128 = json.load(open(K128))
    b131 = json.load(open(K131))

    # ---- control 6: the k128 log must carry item 353's committed board row ----
    got = {"us_strae": b128["us"]["strae"], "us_rho": b128["us"]["rho"],
           "first_strae": b128["first"]["strae"], "n_board": b128["n_board"]}
    if any(abs(got[k] - v) > 1e-9 for k, v in ITEM353.items()):
        print("КОНТРОЛЬ 6 ПРОВАЛЕН: k128 несёт %r, а не числа пункта 353 %r" % (got, ITEM353),
              file=_sys.stderr)
        return 2
    print("контроль 6: k128 несёт ST-RAE %.4f, ро %.4f, первое место %.4f, n=%d пункта 353"
          % (got["us_strae"], got["us_rho"], got["first_strae"], got["n_board"]))

    tr = pd.read_csv(f)
    rng = np.random.default_rng(0)
    out = {"repeats": REPEATS, "kappa_grid": list(KAPPA), "per_enzyme": {},
           "snapshot": b128["snapshot"]}

    # ---- control 4: known answer, homoscedastic Gaussian at a reliability we set ----
    n_s, rel_s = 2000, 0.80
    ts = rng.standard_normal(n_s) * np.sqrt(rel_s)
    ss = np.full(n_s, np.sqrt(1.0 - rel_s))
    Zs1, Zs2 = rng.standard_normal((40, n_s)), rng.standard_normal((40, n_s))
    As, Bs, _ = ceilings(ts, ss, Zs1, Zs2)
    wantA, wantB = gauss_spearman(np.sqrt(rel_s)), gauss_spearman(rel_s)
    if abs(As - wantA) > 0.01 or abs(Bs - wantB) > 0.01:
        print("КОНТРОЛЬ 4 ПРОВАЛЕН: на гауссовых данных при надёжности %.2f получилось "
              "A=%.4f (ждали %.4f), B=%.4f (ждали %.4f)" % (rel_s, As, wantA, Bs, wantB),
              file=_sys.stderr)
        return 2
    print("контроль 4: известный ответ --- гауссовы данные с надёжностью %.2f дают "
          "A %.4f против %.4f и B %.4f против %.4f (закрытая форма (6/pi)arcsin(r/2))"
          % (rel_s, As, wantA, Bs, wantB))

    # ---- control 3: with no noise at all both ceilings must be exactly one ----
    A0, B0, _ = ceilings(ts, np.zeros(n_s), Zs1[:5], Zs2[:5])
    if abs(A0 - 1.0) > 1e-12 or abs(B0 - 1.0) > 1e-12:
        print("КОНТРОЛЬ 3 ПРОВАЛЕН: при нулевом шуме потолки %.12f и %.12f" % (A0, B0),
              file=_sys.stderr)
        return 2
    print("контроль 3: при нулевом шуме оба потолка ровно 1.0 (|Δ| < 1e-12)")

    # ---- 1. the noise, exactly as the organisers publish it -------------------------
    print("\n1. ШУМ, КАК ЕГО ПУБЛИКУЮТ. `{CYP}_pIC50_direct_inhibition_std` --- погрешность фита")
    print("   кривой на то же соединение, по одной на каждую метку, ни одного пропуска.")
    print("   %-8s %6s %9s %9s %9s %9s %8s %9s"
          % ("фермент", "n", "sd метки", "мед. сиг", "rms сиг", "макс сиг", "кор(с,y)",
             "надёжн."))
    dat = {}
    for c in CYPS:
        col = "%s_pIC50_direct_inhibition" % c
        d = tr[tr[col].notna()]
        y = d[col].to_numpy(float)
        sig = d[col + "_std"].to_numpy(float)
        lo = d[col + "_conf_low"].to_numpy(float)
        hi = d[col + "_conf_high"].to_numpy(float)
        if np.isnan(sig).any():
            print("пропуск в _std на %s" % c, file=_sys.stderr)
            return 3
        # ---- control 2: these must be item 354's labels ----
        if abs(y.std() - ITEM354_SD[c]) > 1e-9:
            print("КОНТРОЛЬ 2 ПРОВАЛЕН: sd %s = %.12f, у пункта 354 %.12f"
                  % (c, y.std(), ITEM354_SD[c]), file=_sys.stderr)
            return 2
        mean_var = float(np.mean(sig ** 2))
        t, rel = latent(y, mean_var)
        dat[c] = {"y": y, "sig": sig, "t": t, "rel": rel, "mean_var": mean_var,
                  "width_ratio": float(np.median((hi - lo) / sig))}
        print("   %-8s %6d %9.4f %9.4f %9.4f %9.4f %8.3f %9.4f"
              % (c, len(y), y.std(), np.median(sig), np.sqrt(mean_var), sig.max(),
                 np.corrcoef(sig, y)[0, 1], rel))
        out["per_enzyme"][c] = {"n": len(y), "sd_y": float(y.std()),
                                "median_sigma": float(np.median(sig)),
                                "rms_sigma": float(np.sqrt(mean_var)),
                                "corr_sigma_y": float(np.corrcoef(sig, y)[0, 1]),
                                "reliability": rel, "width_over_std": dat[c]["width_ratio"]}

    # ---- control 7: reproduce item 114's table -- the width is a FUNCTION of the label ----
    print("   контроль 7, пункт 114: ширина полосы --- детерминированная функция метки, а не")
    print("   измерение на соединение. Воспроизводим его таблицу:")
    print("   %-8s %14s %14s %12s %12s"
          % ("фермент", "ро(y,шир)", "п.114", "R2 изотоники", "п.114"))
    for c in CYPS:
        col = "%s_pIC50_direct_inhibition" % c
        d = tr[tr[col].notna()]
        y = d[col].to_numpy(float)
        w = (d[col + "_conf_high"].to_numpy(float) - d[col + "_conf_low"].to_numpy(float))
        rs = float(spearmanr(y, w).statistic)
        fit = IsotonicRegression(increasing=False).fit_transform(y, w)
        r2 = 1.0 - float(np.sum((w - fit) ** 2) / np.sum((w - w.mean()) ** 2))
        want_r, want_r2 = ITEM114[c]
        if abs(rs - want_r) > 0.01 or abs(r2 - want_r2) > 0.01:
            print("КОНТРОЛЬ 7 ПРОВАЛЕН: %s даёт ро %.3f, R2 %.3f против пункта 114 %.3f, %.3f"
                  % (c, rs, r2, want_r, want_r2), file=_sys.stderr)
            return 2
        out["per_enzyme"][c].update(rho_y_width=rs, isotonic_r2=r2)
        print("   %-8s %14.3f %14.3f %12.3f %12.3f" % (c, rs, want_r, r2, want_r2))
    print("   так что ниже `_std` --- это ШКАЛА шума как функция потентности, проверенная")
    print("   пунктом 220 против независимой цепочки, а не информация о каждом соединении.")
    print("   Раздел 3 поэтому варьирует и шкалу, и саму зависимость от метки.")

    # ---- control 1: the band is 2*1.96*_std, so reading it as a 95 % CI is a FACT ----
    bad = {c: v["width_ratio"] for c, v in dat.items() if abs(v["width_ratio"] - 2 * Z95) > 0.05}
    if bad:
        print("КОНТРОЛЬ 1 ПРОВАЛЕН: (hi-lo)/_std = %r, а не %.4f" % (bad, 2 * Z95),
              file=_sys.stderr)
        return 2
    print("   контроль 1: (hi-lo)/_std = " +
          ", ".join("%s %.3f" % (c, dat[c]["width_ratio"]) for c in CYPS) +
          " против 2*1.96 = %.3f --- полоса ЕСТЬ 95 %% ДИ этого же sd, не допущение" % (2 * Z95))

    # ---- 2. the ceiling on Spearman -------------------------------------------------
    print("\n2. ПОТОЛОК НА SPEARMAN. A --- идеальная модель латентной истины против одной")
    print("   измеренной панели; B --- что набрала бы ВТОРАЯ такая же панель. Гомоск. ---")
    print("   тот же средний шум, но постоянный; гаусс --- (6/pi)arcsin(r/2) при той надёжности.")
    print("   %-8s %8s %9s %9s %9s %9s | %8s %8s %7s %7s"
          % ("фермент", "надёжн.", "потолок A", "потолок B", "гомоск.A", "гаусс A",
             "наш ро", "лучший", "наш %", "лучш %"))
    pe = b128["per_enzyme"]
    cA, cB, ours, best = [], [], [], []
    for c in CYPS:
        y, sig, t = dat[c]["y"], dat[c]["sig"], dat[c]["t"]
        n = len(y)
        Z1 = rng.standard_normal((REPEATS, n))
        Z2 = rng.standard_normal((REPEATS, n))
        A, B, sdA = ceilings(t, sig, Z1, Z2)
        Ah, Bh, _ = ceilings(t, np.full(n, np.sqrt(dat[c]["mean_var"])), Z1, Z2)
        gA = gauss_spearman(np.sqrt(dat[c]["rel"]))
        # the board scores 375 compounds, so the same ceiling on a 375-row subsample
        take = [rng.choice(n, size=min(N_BOARD, n), replace=False) for _ in range(REPEATS)]
        A375, B375, _ = ceilings(t, sig, Z1, Z2, take=take)
        o, bst = pe[c]["ours_rho"], pe[c]["best_rho"]
        cA.append(A); cB.append(B); ours.append(o); best.append(bst)
        out["per_enzyme"][c].update(ceiling_A=A, ceiling_B=B, ceiling_A_sd=sdA,
                                    ceiling_A_homo=Ah, ceiling_B_homo=Bh, gauss_A=gA,
                                    ceiling_A_n375=A375, ceiling_B_n375=B375,
                                    ours_rho=o, best_rho=bst)
        print("   %-8s %8.4f %9.4f %9.4f %9.4f %9.4f | %8.4f %8.4f %6.0f%% %6.0f%%"
              % (c, dat[c]["rel"], A, B, Ah, gA, o, bst, 100 * o / A, 100 * bst / A))
    mA, mB = float(np.mean(cA)), float(np.mean(cB))
    mo, mb = float(np.mean(ours)), float(np.mean(best))
    print("   %-8s %8s %9.4f %9.4f %9s %9s | %8.4f %8.4f %6.0f%% %6.0f%%"
          % ("МАКРО", "", mA, mB, "", "", mo, mb, 100 * mo / mA, 100 * mb / mA))
    print("   контроль: на подвыборке n=375, как на живой доске, потолок A = " +
          ", ".join("%.4f" % out["per_enzyme"][c]["ceiling_A_n375"] for c in CYPS))

    # ---- 3. how far off the ceiling has to be for first place to be AT it -----------
    print("\n3. ЧУВСТВИТЕЛЬНОСТЬ ПО ШКАЛЕ, в обе стороны. Пункт 220 измерил, что наша полоса живёт")
    print("   на шкале ВОСПРОИЗВОДИМОСТИ, в 3.8 раза выше собственной se подгонки организаторов")
    print("   (0.090 против 0.024 pIC50), так что kappa=1 --- уже щедрое чтение, а kappa=0.27 ---")
    print("   чтение по шкале точности. kappa --- множитель на sigma, те же случайные draw'ы:")
    print("   %-7s %9s %9s %9s %9s %9s" % ("kappa", "1A2", "2C9", "2D6", "3A4", "МАКРО"))
    sweep = {}
    Zf = {c: (rng.standard_normal((REPEATS_FAST, len(dat[c]["y"]))),
              rng.standard_normal((REPEATS_FAST, len(dat[c]["y"])))) for c in CYPS}

    def macro_at(k):
        v = [ceilings(dat[c]["t"], dat[c]["sig"], Zf[c][0], Zf[c][1], kappa=k)[0] for c in CYPS]
        return float(np.mean(v)), v

    for k in KAPPA:
        m, v = macro_at(k)
        sweep[k] = {"per_enzyme": v, "macro": m}
        print("   %-7.2f %9.4f %9.4f %9.4f %9.4f %9.4f" % (k, v[0], v[1], v[2], v[3], m))
    # ---- control 5: the ceiling must fall strictly as the noise is scaled up ----
    ms = [sweep[k]["macro"] for k in KAPPA]
    if any(ms[i + 1] >= ms[i] for i in range(len(ms) - 1)):
        print("КОНТРОЛЬ 5 ПРОВАЛЕН: потолок не монотонен по kappa: %r" % ms, file=_sys.stderr)
        return 2
    print("   контроль 5: потолок строго падает по kappa, %.4f -> %.4f" % (ms[0], ms[-1]))

    # the label-independent arm: item 221's instrument floor, constant across compounds
    flo = [ceilings(dat[c]["t"], np.full(len(dat[c]["y"]), FLOOR221),
                    Zf[c][0], Zf[c][1])[0] for c in CYPS]
    print("   %-7s %9.4f %9.4f %9.4f %9.4f %9.4f  <-- постоянная sigma = %.2f pIC50, пункт 221"
          % ("пол", flo[0], flo[1], flo[2], flo[3], float(np.mean(flo)), FLOOR221))
    out["floor221"] = {"sigma": FLOOR221, "per_enzyme": flo, "macro": float(np.mean(flo))}

    kstar = {}
    for lab, target in (("наш", mo), ("лучший на доске", mb)):
        lo_, hi_ = 1.0, 40.0
        if macro_at(hi_)[0] > target:
            kstar[lab] = float("inf")
        else:
            for _ in range(28):
                mid = 0.5 * (lo_ + hi_)
                if macro_at(mid)[0] > target:
                    lo_ = mid
                else:
                    hi_ = mid
            kstar[lab] = 0.5 * (lo_ + hi_)
        print("   чтобы потолок опустился до %s макро-ро %.4f, шкала полосы должна быть занижена "
              "в %s раз" % (lab, target, ("%.1f" % kstar[lab]) if np.isfinite(kstar[lab])
                            else ">40"))

    # ---- 4. the test's wider labels raise the ceiling, not lower it ------------------
    print("\n4. НА ТЕСТЕ ПОТОЛОК ВЫШЕ. Пункт 354 ограничил sd меток теста снизу через")
    print("   MAE/sqrt(1-R2); при том же шуме более широкая истина --- более высокая надёжность:")
    print("   %-8s %9s %9s %9s %9s" % ("фермент", "sd обуч.", "sd тест >=", "надёжн.", "потолок A"))
    tA = []
    for c in CYPS:
        bnd = b131["spread"][c]["bound"]
        t2, rel2 = latent(dat[c]["y"], dat[c]["mean_var"], var_target=bnd ** 2)
        n = len(dat[c]["y"])
        A2, _, _ = ceilings(t2, dat[c]["sig"], Zf[c][0], Zf[c][1])
        tA.append(A2)
        out["per_enzyme"][c].update(test_sd_bound=bnd, reliability_test=rel2, ceiling_A_test=A2)
        print("   %-8s %9.4f %9.4f %9.4f %9.4f" % (c, dat[c]["sd_y"] if False else
                                                   dat[c]["y"].std(), bnd, rel2, A2))
    print("   %-8s %9s %9s %9s %9.4f" % ("МАКРО", "", "", "", float(np.mean(tA))))

    # ---- 5. the arithmetic the goal actually needs -----------------------------------
    need = b128["rho_needed_within"]
    avail = mA - mo
    print("\n5. АРИФМЕТИКА ЦЕЛИ. Пункту 353 нужно +%.4f макро-ро, чтобы закрыть отрыв до первого"
          % need)
    print("   места по наклону пункта 322 (%.2f ST-RAE на единицу ро). Потолок оставляет"
          % SLOPE_WITHIN)
    print("   незанятыми +%.4f. Требуется %.0f %% того, что метки вообще позволяют."
          % (avail, 100 * need / avail))
    print("   по ферментам, сколько ро не занято: " +
          ", ".join("%s +%.3f" % (c, out["per_enzyme"][c]["ceiling_A"] -
                                  out["per_enzyme"][c]["ours_rho"]) for c in CYPS))

    # ---- 6. WHERE the noise is, and it is not spread evenly over the ranking ---------
    print("\n6. ГДЕ ШУМ. кор(сигма, y) выше --- от -0.74 до -0.90: погрешность фита почти")
    print("   полностью на СЛАБЫХ соединениях. Значит потолок по половинам метки разный, и")
    print("   сравнивать надо подобное с подобным --- наше ро внутри той же половины:")
    print("   %-8s | %8s %9s %7s | %8s %9s %7s | %7s"
          % ("фермент", "сил. ро", "потолок", "доля", "слаб. ро", "потолок", "доля", "дисп.%"))
    oofp = _pl.Path(RES + "preds/oof_submitted.json")
    if not oofp.exists():
        print("нет %s --- раздел 6 пропущен" % oofp, file=_sys.stderr)
    else:
        rows = pd.read_csv(D + "rows.csv")
        al = (tr.set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
        P = [np.asarray(v, float) for v in json.load(open(oofp))["P"]]
        halves = {}
        for e, c in enumerate(CYPS):
            col = "%s_pIC50_direct_inhibition" % c
            m = al[col].notna().to_numpy()
            ya = al[col].to_numpy(float)[m]
            sa = al[col + "_std"].to_numpy(float)[m]
            pa = P[e]
            if len(pa) != len(ya):
                print("раздел 6: %s --- %d предсказаний против %d меток, пропуск"
                      % (c, len(pa), len(ya)), file=_sys.stderr)
                continue
            med = float(np.median(ya))
            cut = {"сильные": ya >= med, "слабые": ya < med}
            row, rec = [], {}
            for lab, k in cut.items():
                mv = float(np.mean(sa[k] ** 2))
                th, relh = latent(ya[k], mv)
                Zh1 = rng.standard_normal((REPEATS_FAST, int(k.sum())))
                Zh2 = rng.standard_normal((REPEATS_FAST, int(k.sum())))
                Ah2, _, _ = ceilings(th, sa[k], Zh1, Zh2)
                oh = rho(ya[k], pa[k])
                rec[lab] = {"n": int(k.sum()), "reliability": relh, "ceiling_A": Ah2,
                            "ours_rho": oh, "rms_sigma": float(np.sqrt(mv))}
                row += [oh, Ah2, oh / Ah2]
            share = float(np.sum(sa[cut["слабые"]] ** 2) / np.sum(sa ** 2))
            rec["weak_var_share"] = share
            rec["ours_rho_full"] = rho(ya, pa)
            halves[c] = rec
            print("   %-8s | %8.4f %9.4f %6.0f%% | %8.4f %9.4f %6.0f%% | %6.0f%%"
                  % (c, row[0], row[1], 100 * row[2], row[3], row[4], 100 * row[5],
                     100 * share))
        out["halves"] = halves
        if halves:
            sp = float(np.mean([h["сильные"]["ours_rho"] / h["сильные"]["ceiling_A"]
                                for h in halves.values()]))
            wp = float(np.mean([h["слабые"]["ours_rho"] / h["слабые"]["ceiling_A"]
                                for h in halves.values()]))
            print("   в среднем мы берём %.0f %% потолка на СИЛЬНОЙ половине и %.0f %% на слабой;"
                  % (100 * sp, 100 * wp))
            print("   слабая половина несёт %.0f %% всей дисперсии шума."
                  % (100 * np.mean([h["weak_var_share"] for h in halves.values()])))
            print("   справка, наше ро ВНЕ ФОЛДА на всей выборке: " +
                  ", ".join("%s %.4f" % (c, halves[c]["ours_rho_full"]) for c in halves) +
                  ", макро %.4f --- против %.4f на доске, то есть тест --- более лёгкая"
                  % (float(np.mean([h["ours_rho_full"] for h in halves.values()])), mo))
            print("   выборка (пункт 357), и доли выше поэтому --- оценка снизу на нашу долю")
            print("   потолка на доске; СРАВНЕНИЕ половин друг с другом от этого не страдает.")
            out["macro_halves"] = {"share_strong": sp, "share_weak": wp}

    out["macro"] = {"ceiling_A": mA, "ceiling_B": mB, "ours_rho": mo, "best_rho": mb,
                    "ceiling_A_test": float(np.mean(tA)),
                    "share_ours": mo / mA, "share_best": mb / mA,
                    "rho_needed_within": need, "rho_available": avail,
                    "need_over_available": need / avail}
    out["kappa_sweep"] = {str(k): v for k, v in sweep.items()}
    out["kappa_star"] = {k: (None if not np.isfinite(v) else v) for k, v in kstar.items()}
    out["known_answer"] = {"reliability": rel_s, "A": As, "A_closed": wantA,
                           "B": Bs, "B_closed": wantB}
    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
    print("\nзапись: %s" % LOG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
