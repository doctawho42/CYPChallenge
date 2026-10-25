"""What fraction of a hit's close analogues is itself a hit? Answered on the PRIMARY SCREEN, which
is the only readout in this project that is free of promotion bias.

Item 379 measured the parent-to-child law on dose-response labels and the answer contradicted the
board: it put CYP1A2 at 5.399 and CYP2D6 at 5.355, higher still than the file that had just lost 71
places. The defect is selection. A compound carries a pIC50 only if it was promoted to dose-response
for that isoform --- 1412, 1285, 1493 and 2335 of 4905 --- and promotion happened BECAUSE the
compound screened active. Conditioning on an active parent inside a set where everything already
passed a screen measures hit-to-hit pairs, while the test's children are catalogue compounds bought
for their similarity and never screened at all.

The single-concentration screen has no such hole: `log2fc_estimate` exists for all 4376 screened
compounds whether or not they were ever promoted, so "analogue of a hit" can be sampled without
conditioning on the outcome.

THE QUANTITY. Over pairs at the test's own similarity (0.5156--0.7018, item 378's 10th--90th
percentile of parent similarity), the probability that a compound is a hit GIVEN that its neighbour
is one. Everything the test's level depends on follows from it: the 750 are ten analogues each of
about 75 parents, so the expected share of actives is that conditional rate, and the rest sit at or
below the assay floor.

THREE CONTROLS:

  * the screen's direction is established from the data, not assumed: `log2fc_estimate` must
    correlate with pIC50 on the compounds that have both, and the sign of that correlation decides
    which tail is "active";
  * the hit threshold on the screen is calibrated against PROMOTION, not chosen -- the cut that
    best separates compounds that received a pIC50 from those that did not;
  * the conditional rate must beat the base rate, and the SAME statistic computed for analogues of
    NON-hits must not, or "analogue of a hit" carries no information and the derivation is void.

Writes results/logs/k157_hitrate.json.
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


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    Y = np.stack([tr["%s_pIC50_direct_inhibition" % c].to_numpy(float) for c in CYPS], 1)

    # screen value per (compound, enzyme), aligned to rows.csv; NaN where not screened
    L = np.full((len(rows), len(CYPS)), np.nan)
    pos = {n: i for i, n in enumerate(rows.Molecule_Name)}
    for e, c in enumerate(CYPS):
        g = sc[sc.enzyme == c]
        idx = g.Molecule_Name.map(pos)
        ok = idx.notna()
        L[idx[ok].astype(int).to_numpy(), e] = g.log2fc_estimate.to_numpy()[ok.to_numpy()]
    print("  скрининг выровнен: по ферментам есть у %s из %d"
          % ("/".join(str(int(np.isfinite(L[:, e]).sum())) for e in range(4)), len(rows)))

    print("\n  КОНТРОЛЬ направления: связь скрининга с pIC50 там, где есть оба")
    sign = {}
    for e, c in enumerate(CYPS):
        m = np.isfinite(L[:, e]) & np.isfinite(Y[:, e])
        r = float(spearmanr(L[m, e], Y[m, e]).statistic)
        sign[c] = r
        print("     %-8s n=%4d  Spearman(log2fc, pIC50) = %+.3f  -> активнее = %s"
              % (c, m.sum(), r, "НИЖЕ log2fc" if r < 0 else "ВЫШЕ log2fc"))
    if not all(v < -0.1 for v in sign.values()):
        raise SystemExit("ОТКАЗ: направление скрининга не одинаково по ферментам: %s" % sign)

    print("\n  КОНТРОЛЬ порога: калибровка по ПРОДВИЖЕНИЮ в дозовые кривые")
    thr = {}
    for e, c in enumerate(CYPS):
        m = np.isfinite(L[:, e])
        promoted = np.isfinite(Y[:, e])[m]
        v = L[m, e]
        grid = np.quantile(v, np.linspace(0.01, 0.99, 99))
        # Youden's J on "log2fc <= t predicts promotion"
        j = [( (v <= t)[promoted].mean() - (v <= t)[~promoted].mean() ) for t in grid]
        t = float(grid[int(np.argmax(j))])
        thr[c] = {"threshold": t, "youden": float(np.max(j)),
                  "base_rate": float((v <= t).mean()),
                  "promoted_share": float(promoted.mean())}
        print("     %-8s порог log2fc <= %+.3f (J=%.3f); доля «хитов» %.3f, продвинутых %.3f"
              % (c, t, np.max(j), thr[c]["base_rate"], thr[c]["promoted_share"]))

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    F = []
    for s in rows.SMILES:
        mm = Chem.MolFromSmiles(s)
        if mm is None:
            raise SystemExit("RDKit не разобрал %r" % s)
        F.append(gen.GetFingerprint(mm))
    S = np.vstack([np.asarray(DataStructs.BulkTanimotoSimilarity(f, F), dtype=np.float32)
                   for f in F])
    assert np.allclose(S.max(1), 1.0), "КОНТРОЛЬ: максимум с собой не 1"
    iu = np.triu_indices(len(rows), k=1)
    sv = S[iu]
    sel = (sv >= BAND[0]) & (sv <= BAND[1])
    I, J = iu[0][sel], iu[1][sel]
    print("\n  пар в полосе %.4f..%.4f: %d" % (*BAND, len(I)))

    print("\n  ДОЛЯ ХИТОВ СРЕДИ АНАЛОГОВ, на скрининге (без смещения продвижения)")
    print("  %-8s %7s %10s %12s %12s %10s"
          % ("фермент", "n пар", "базовая", "аналог хита", "аналог НЕхита", "подъём"))
    out = {}
    for e, c in enumerate(CYPS):
        t = thr[c]["threshold"]
        ok = np.isfinite(L[I, e]) & np.isfinite(L[J, e])
        a, b = L[I[ok], e], L[J[ok], e]
        x = np.concatenate([a, b])            # "parent"
        y = np.concatenate([b, a])            # "child"
        hit_x, hit_y = x <= t, y <= t
        base = float(thr[c]["base_rate"])
        cond = float(hit_y[hit_x].mean()) if hit_x.sum() else np.nan
        anti = float(hit_y[~hit_x].mean()) if (~hit_x).sum() else np.nan
        out[c] = {"n_pairs": int(ok.sum()), "base_rate": base, "analog_of_hit": cond,
                  "analog_of_nonhit": anti, "lift": cond / base if base else np.nan,
                  "n_hit_parents": int(hit_x.sum())}
        print("  %-8s %7d %10.3f %12.3f %12.3f %10.2f"
              % (c, ok.sum(), base, cond, anti, out[c]["lift"]))
    bad = [c for c in CYPS if not (out[c]["analog_of_hit"] > out[c]["base_rate"] >
                                   out[c]["analog_of_nonhit"])]
    print("  КОНТРОЛЬ: аналог хита > базовой > аналога нехита на %d из 4%s"
          % (4 - len(bad), (" (нарушено: %s)" % ", ".join(bad)) if bad else ""))

    print("\n  СЛЕДСТВИЕ для 750: ожидаемая доля активных и медиана pIC50 активных")
    for e, c in enumerate(CYPS):
        act = Y[np.isfinite(Y[:, e]), e]
        out[c]["active_pIC50_median"] = float(np.median(act))
        print("     %-8s ожидаемая доля активных %.3f; у активных медиана pIC50 %.3f"
              % (c, out[c]["analog_of_hit"], out[c]["active_pIC50_median"]))
    json.dump({"band": BAND, "direction": sign, "thresholds": thr, "rates": out},
              open(RES + "logs/k157_hitrate.json", "w"), ensure_ascii=False, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
