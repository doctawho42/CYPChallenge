"""Does the layer disagreement grow with lipophilicity, and does it exist at all under the fitted E.

Two questions, and the second is prior to the first.

**The hypothesis under test.** Observed pIC50 is a property of the molecule *and* the assay: free
concentration at the enzyme is below nominal because of non-specific binding to microsomal protein
and lipid, and the free fraction is a known decreasing function of lipophilicity (Hallifax-Houston).
Lipophilic compounds should therefore appear weaker than they are, and the distortion should grow
with logP. Combined with Cheng-Prusoff, which contributes a constant per-enzyme offset, this
predicts a **positive** correlation between the shift `d` and logP.

**The confounder that must be removed, and it was named by the proposer.** `d` is defined as the
shift that returns the Hill slope to one, and it contains `-pIC50` by construction, while pIC50
itself correlates with logP. So the raw correlation is guaranteed to be positive whether or not the
hypothesis is true. The control is to remove pIC50 from `d` isotonically and correlate the residual.

**The prior question this file must ask first.** Item 170 measured that `d` is computed under
`E = 1`, and that under this repository's own fitted E the median Hill slope is 1.102, 0.949, 1.000
and 0.641 -- essentially one on three enzymes out of four. If the slope is already one, the shift
that returns it to one is zero, and there is no disagreement to explain, correlate or model. So `d`
is recomputed here under both amplitudes before anything is regressed against logP. If it collapses
under the fitted E, then the quantity whose lipophilicity dependence is being tested is an artefact
of the amplitude assumption, and both the hypothesis and its refutation are about nothing.

Pre-registered, all three outcomes:

    d выживает при подогнанном E и растёт с logP после контроля -> связывание подтверждено
    d выживает, но с logP после контроля не связан или связан обратно -> расхождение реально,
        механизм не микросомальное связывание
    d схлопывается при подогнанном E -> расхождения нет, обсуждать нечего

Reads the inhibition, Emax and screening tables plus data/feats.npz for MolLogP. Writes nothing.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
PC0 = 4.305
FIT_E = {"CYP1A2": 0.728, "CYP2C9": 0.621, "CYP2D6": 0.867, "CYP3A4": 0.931}
MIN_GAP = 0.5


def shift_to_unit(u, gap):
    """Сдвиг s, при котором медиана u/(gap - s) равна единице. Сеткой: медиана разрывна."""
    grid = np.linspace(-2.5, 2.5, 5001)
    with np.errstate(divide="ignore", invalid="ignore"):
        meds = np.array([np.median(u / (gap - s)) for s in grid])
    ok = np.isfinite(meds)
    if not ok.any():
        return np.nan, np.nan
    j = int(np.argmin(np.abs(meds[ok] - 1.0)))
    return float(grid[ok][j]), float(abs(meds[ok][j] - 1.0))


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    z = np.load(D + "feats.npz")
    dn = [l.strip() for l in open(D + "desc_names.csv")]
    logp_all = z["DESC"][:, dn.index("MolLogP")].astype(float)

    print("=== 0. существует ли сдвиг при ПОДОГНАННОМ E (пункт 170) ===")
    print(f"{'фермент':8s} {'n':>5s} {'мед h при E=1':>14s} {'сдвиг при E=1':>14s} "
          f"{'мед h при E подогн':>19s} {'сдвиг при E подогн':>19s} {'промах':>8s}")
    keep = {}
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        sub = sc[sc.enzyme == c].set_index("Molecule_Name")["log2fc_estimate"]
        pi = tr[col].to_numpy(float)
        l2 = sub.reindex(rows.Molecule_Name).to_numpy(float)
        a = 2.0 ** l2
        gap = PC0 - pi
        ok = np.isfinite(a) & np.isfinite(pi) & (a > 1e-6) & (a < 1 - 1e-6) & (np.abs(gap) > MIN_GAP)
        E = FIT_E[c]
        numE = E - 1 + a
        okE = ok & (numE > 1e-9)
        u1 = np.log10(a[ok] / (1 - a[ok]))
        uE = np.log10(numE[okE] / (1 - a[okE]))
        s1, m1 = shift_to_unit(u1, gap[ok])
        sE, mE = shift_to_unit(uE, gap[okE])
        h1 = u1 / gap[ok]
        hE = uE / gap[okE]
        print(f"{c:8s} {int(ok.sum()):5d} {np.median(h1):14.3f} {s1:+14.3f} "
              f"{np.median(hE):19.3f} {sE:+19.3f} {mE:8.3f}")
        keep[c] = dict(ok=ok, okE=okE, h1=h1, hE=hE, gap=gap, u1=u1, uE=uE, pi=pi)

    print("\n=== 1. сырая связь сдвига с липофильностью, и контроль на потенцию ===")
    print("  d определён поточечно как s_i = gap_i - u_i/1 = gap_i - u_i, то есть сдвиг,")
    print("  при котором ЭТО соединение даёт наклон единицу.")
    print(f"\n{'фермент':8s} {'rho(d,logP)':>12s} {'rho(pIC50,logP)':>16s} "
          f"{'rho(остаток,logP)':>18s} {'d при logP<2':>13s} {'d при logP>4':>13s}")
    for c in CYPS:
        k = keep[c]
        ok = k["ok"]
        lp = logp_all[ok]
        pi = k["pi"][ok]
        d = k["gap"][ok] - k["u1"]          # s, при котором h_i = 1
        raw = spearmanr(d, lp).statistic
        rpp = spearmanr(pi, lp).statistic
        # Контроль: убрать потенцию изотоникой, коррелировать остаток.
        # increasing="auto" ОБЯЗАТЕЛЕН: d убывает по pIC50 (содержит -pi), а умолчание
        # increasing=True подгоняет неубывающую функцию к убывающей связи и даёт
        # почти константу --- остаток тогда равен d со сдвигом, и корреляция не меняется.
        iso = IsotonicRegression(increasing="auto", out_of_bounds="clip").fit(pi, d)
        res = d - iso.predict(pi)
        ctl = spearmanr(res, lp).statistic
        lo, hi = lp < 2, lp > 4
        print(f"{c:8s} {raw:+12.3f} {rpp:+16.3f} {ctl:+18.3f} "
              f"{np.median(d[lo]) if lo.sum() > 20 else np.nan:13.3f} "
              f"{np.median(d[hi]) if hi.sum() > 20 else np.nan:13.3f}")

    print("\n=== 2. то же при ПОДОГНАННОМ E ===")
    print(f"{'фермент':8s} {'rho(d,logP)':>12s} {'rho(остаток,logP)':>18s} {'медиана d':>11s}")
    for c in CYPS:
        k = keep[c]
        okE = k["okE"]
        lp = logp_all[okE]
        pi = k["pi"][okE]
        d = k["gap"][okE] - k["uE"]
        iso = IsotonicRegression(increasing="auto", out_of_bounds="clip").fit(pi, d)
        res = d - iso.predict(pi)
        print(f"{c:8s} {spearmanr(d, lp).statistic:+12.3f} "
              f"{spearmanr(res, lp).statistic:+18.3f} {np.median(d):11.3f}")

    print("""
Как читать. Раздел 0 --- предварительный и решающий. Столбец «сдвиг при E подогн» рядом со
столбцом «мед h при E подогн»: если медианный наклон уже единица, то сдвигать нечего, и
величина, чью связь с липофильностью проверяют, существует только при допущении E = 1.

Раздел 1 воспроизводит проверку и её контроль. Сырая корреляция обязана быть положительной
по построению --- d содержит -pIC50, а pIC50 растёт с logP, и второй столбец показывает
насколько. Решает третий: знак и величина ПОСЛЕ снятия потенции. Отрицательный знак означает,
что при равной потенции более липофильные расходятся СЛАБЕЕ --- противоположно неспецифическому
связыванию, которое предсказывало обратное.

Раздел 2 повторяет то же в параметризации, где прибор физически осмыслен.""")


if __name__ == "__main__":
    main()
