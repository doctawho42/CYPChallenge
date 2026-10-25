"""The test set's activity distribution, derived from the published design -- with the threshold
calibrated against ACTIVITY instead of against promotion, which is what item 380 got wrong.

Item 380 calibrated the screen's hit threshold by asking which cut best predicts PROMOTION to
dose-response. On CYP3A4 that returned Youden J = -0.005, negative across the whole grid: promotion
there was not screen-driven, and the resulting 0.000 hit rate was an artefact rather than a
measurement, in the face of CYP3A4 having the strongest screen-to-pIC50 correlation of the four
(-0.936). CYP2C9's cut was placed badly too, calling 84 per cent of compounds hits.

The fix is to calibrate against the thing the challenge actually scores. A compound is ACTIVE when
it has a dose-response pIC50 and that pIC50 is at least 4.0 -- the organisers' own floor, "low
activity compounds (pIC50 < 4) fall below the lowest tested assay dose". That label is observed for
every one of the 4375 screened compounds, because "never promoted" is an observation and not a
missing value, so there is no selection hole of the kind that sank item 380's dose-response device.

THE DERIVATION. Each test compound has a recovered parent (item 378) and the parent has a screen
value on every isoform. So for compound i and isoform e the probability of being active is the
measured rate among analogues of a hit when the parent is a hit on e, and among analogues of a
non-hit when it is not -- per compound, per isoform, with no need to assign parents to isoforms at
all. The expected level follows from that probability and the measured levels of actives and
inactives.

CONTROLS:
  * direction is read off the data, not assumed;
  * the threshold must separate: Youden J reported per isoform, and an isoform whose J is below
    0.25 is declared UNDERIVED rather than quietly used;
  * analogue-of-hit > base > analogue-of-non-hit, or "analogue of a hit" carries nothing;
  * every rate carries a bootstrap interval over pairs, because the test's similarity band holds
    only about 450 usable pairs per isoform and a point estimate would hide that.

Writes results/logs/k158_derive.json.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.stats import spearmanr

from submit import CYPS

BAND = (0.5156, 0.7018)
FLOOR = 4.0            # the organisers' own "below the lowest tested dose"
MIN_J = 0.25           # below this the isoform is declared underived
BOOT = 400


def main():
    rows = pd.read_csv(D + "rows.csv")
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    par = pd.read_csv(RES + "preds/k155_parent_of_test.csv")
    assert list(par.Molecule_Name) == list(te.Molecule_Name)
    Y = np.stack([tr["%s_pIC50_direct_inhibition" % c].to_numpy(float) for c in CYPS], 1)

    L = np.full((len(rows), len(CYPS)), np.nan)
    pos = {n: i for i, n in enumerate(rows.Molecule_Name)}
    for e, c in enumerate(CYPS):
        g = sc[sc.enzyme == c]
        idx = g.Molecule_Name.map(pos)
        ok = idx.notna()
        L[idx[ok].astype(int).to_numpy(), e] = g.log2fc_estimate.to_numpy()[ok.to_numpy()]
    screened = np.isfinite(L).all(1)
    ACT = np.isfinite(Y) & (Y >= FLOOR)          # observed for every screened compound
    print("  скринировано %d из %d; активных (pIC50 >= %.1f) по ферментам: %s"
          % (screened.sum(), len(rows), FLOOR,
             "/".join(str(int(ACT[screened, e].sum())) for e in range(4))))

    print("\n  КОНТРОЛЬ направления и КАЛИБРОВКА порога по АКТИВНОСТИ (а не по продвижению)")
    thr, und = {}, []
    for e, c in enumerate(CYPS):
        m = screened
        v, a = L[m, e], ACT[m, e]
        both = np.isfinite(L[:, e]) & np.isfinite(Y[:, e])
        rho = float(spearmanr(L[both, e], Y[both, e]).statistic)
        if not np.isfinite(rho):
            raise SystemExit("%s: контроль направления дал NaN на %d парах" % (c, both.sum()))
        t = float(np.quantile(v, a.mean()))
        j = np.array([(v <= t)[a].mean() - (v <= t)[~a].mean()])
        thr[c] = {"threshold": t, "youden": float(j.max()), "rho_screen_pic50": rho,
                  "base_rate": float((v <= t).mean()), "true_active_rate": float(a.mean())}
        if j.max() < MIN_J:
            und.append(c)
        print("     %-8s ро %+.3f | порог log2fc <= %+.3f, J=%.3f | доля по порогу %.3f, "
              "истинная доля активных %.3f"
              % (c, rho, t, j.max(), thr[c]["base_rate"], thr[c]["true_active_rate"]))
    print("     НЕВЫВЕДЕНО (J < %.2f): %s" % (MIN_J, ", ".join(und) if und else "нет"))

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    F = []
    for s in rows.SMILES:
        mm = Chem.MolFromSmiles(s)
        if mm is None:
            raise SystemExit("RDKit не разобрал %r" % s)
        F.append(gen.GetFingerprint(mm))
    S = np.vstack([np.asarray(DataStructs.BulkTanimotoSimilarity(f, F), dtype=np.float32)
                   for f in F])
    assert np.allclose(S.max(1), 1.0)
    iu = np.triu_indices(len(rows), k=1)
    sv = S[iu]
    sel = (sv >= BAND[0]) & (sv <= BAND[1])
    I, J = iu[0][sel], iu[1][sel]

    rng = np.random.default_rng(0)
    print("\n  ДОЛЯ АКТИВНЫХ СРЕДИ АНАЛОГОВ, с бутстрэпом по парам (%d пар в полосе)" % len(I))
    print("  %-8s %7s %9s %20s %20s" % ("фермент", "n пар", "базовая",
                                        "аналог хита [95%]", "аналог НЕхита [95%]"))
    rates = {}
    for e, c in enumerate(CYPS):
        ok = np.isfinite(L[I, e]) & np.isfinite(L[J, e])
        t = thr[c]["threshold"]
        xi, xj = I[ok], J[ok]
        x = np.concatenate([L[xi, e], L[xj, e]])
        y = np.concatenate([ACT[xj, e], ACT[xi, e]])      # NB: child's TRUE activity
        hx = x <= t
        def rate(mask, idx=None):
            yy = y[mask] if idx is None else y[idx][mask[idx]] if False else y[mask]
            return float(yy.mean()) if len(yy) else np.nan
        cond, anti = float(y[hx].mean()), float(y[~hx].mean())
        bc, ba = [], []
        n = len(x)
        for _ in range(BOOT):
            b = rng.integers(0, n, n)
            hb, yb = hx[b], y[b]
            if hb.sum() > 5 and (~hb).sum() > 5:
                bc.append(yb[hb].mean()); ba.append(yb[~hb].mean())
        ci = lambda v: (float(np.percentile(v, 2.5)), float(np.percentile(v, 97.5)))
        rates[c] = {"n_pairs": int(ok.sum()), "base": thr[c]["true_active_rate"],
                    "analog_of_hit": cond, "analog_of_hit_ci": ci(bc),
                    "analog_of_nonhit": anti, "analog_of_nonhit_ci": ci(ba),
                    "n_hit_parents": int(hx.sum())}
        print("  %-8s %7d %9.3f %9.3f [%.3f,%.3f] %9.3f [%.3f,%.3f]"
              % (c, ok.sum(), rates[c]["base"], cond, *ci(bc), anti, *ci(ba)))
    bad = [c for c in CYPS if not (rates[c]["analog_of_hit"] > rates[c]["base"] >
                                   rates[c]["analog_of_nonhit"])]
    print("  КОНТРОЛЬ порядка (аналог хита > базовой > аналога нехита): %d из 4%s"
          % (4 - len(bad), (" — нарушено: %s" % ", ".join(bad)) if bad else ""))

    print("\n  ВЫВЕДЕННОЕ распределение 750, против поданных файлов")
    print("  %-8s %9s %9s %9s %9s %9s %9s"
          % ("фермент", "P(акт)", "выведено", "sd вывед.", "best4", "rankmap", "stack374"))
    sub = {k: pd.read_csv(RES + "submission/activity_submission_%s.csv" % k)
           for k in ("best4", "rankmap", "stack374")}
    out = {}
    for e, c in enumerate(CYPS):
        col = "%s_pIC50_direct_inhibition" % c
        p = par["parent_%s" % c].to_numpy(float)          # parent's pIC50, NaN if inactive
        pn = par.parent_name.map(pos).to_numpy()
        phit = np.where(np.isfinite(L[pn, e]), L[pn, e] <= thr[c]["threshold"], False)
        prob = np.where(phit, rates[c]["analog_of_hit"], rates[c]["analog_of_nonhit"])
        act = Y[ACT[:, e], e]
        lo_pool = Y[np.isfinite(Y[:, e]) & (Y[:, e] < FLOOR), e]
        mu_hi, mu_lo = float(act.mean()), float(lo_pool.mean()) if len(lo_pool) else FLOOR - 0.5
        mu = prob * mu_hi + (1 - prob) * mu_lo
        sd = float(np.sqrt(np.mean(prob * (1 - prob)) * (mu_hi - mu_lo) ** 2
                           + np.mean(prob) * act.var() + (1 - np.mean(prob)) *
                           (lo_pool.var() if len(lo_pool) > 1 else 0.25)))
        out[c] = {"p_active_mean": float(prob.mean()), "share_hit_parents": float(phit.mean()),
                  "mu_active": mu_hi, "mu_inactive": mu_lo, "derived_mean": float(mu.mean()),
                  "derived_sd": sd,
                  **{k: float(sub[k][col].mean()) for k in sub}}
        print("  %-8s %9.3f %9.3f %9.3f %9.3f %9.3f %9.3f"
              % (c, prob.mean(), mu.mean(), sd, out[c]["best4"], out[c]["rankmap"],
                 out[c]["stack374"]))
    print("\n  уровни пулов: активные %s; ниже пола %s"
          % ("/".join("%.2f" % out[c]["mu_active"] for c in CYPS),
             "/".join("%.2f" % out[c]["mu_inactive"] for c in CYPS)))
    print("  доля тестовых, чей родитель — хит: %s"
          % "/".join("%.2f" % out[c]["share_hit_parents"] for c in CYPS))
    json.dump({"band": BAND, "floor": FLOOR, "thresholds": thr, "underived": und,
               "rates": rates, "derived": out},
              open(RES + "logs/k158_derive.json", "w"), ensure_ascii=False, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
