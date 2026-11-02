"""How much does a parent's activity tell us about its analogue's? Measured on ALL labelled pairs at
the test's own similarity, then applied to the 75 recovered parents to DERIVE the test set's
activity distribution.

Item 377 cost 71 places by putting the level in the wrong place, and item 378 recovered why the
level is derivable at all: the test is about 75 series of ten catalogue analogues around known
training compounds --- children of hits, not a fresh chemical sample. So the question that sets the
level is empirical and answerable without a single test label: at the similarity the test actually
sits at, how far does a child regress from its parent?

THE DEVICE, and a first version of this file got it wrong. Restricting to NEAREST-NEIGHBOUR pairs
left 131 labelled pairs on CYP1A2, which is too few to say anything. The question is not about
nearest neighbours; it is about any two compounds at Tanimoto ~0.59, so every in-band pair counts
and there are tens of thousands of them.

TWO READINGS, because the second does not assume the first's shape:

  наклон   a least-squares line of child on parent over all in-band pairs. Interpretable: how much
           of a parent's potency survives into an analogue, and what the residual spread is.
  условное the non-parametric version and the one the derivation uses: take in-band pairs whose
           PARENT is as active as the real parents are, and read off the distribution of the other
           member directly. No linearity, no normality.

TWO CONTROLS:

  * SHUFFLING the pairing must destroy the slope, or a positive slope is just both members being
    drawn from one labelled distribution;
  * the conditional mean must differ from the unconditional one, or "conditioning on an active
    parent" is doing nothing and the test is simply the training distribution.

Reads results/preds/k155_parent_of_test.csv. Writes results/logs/k156_childlaw.json.
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

from submit import CYPS

BAND = (0.5156, 0.7018)      # the test's own 10th-90th percentile of parent similarity, from k155


def main():
    rows = pd.read_csv(D + "rows.csv")
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    par = pd.read_csv(RES + "preds/k155_parent_of_test.csv")
    assert list(par.Molecule_Name) == list(te.Molecule_Name), "родители не в порядке теста"
    Y = np.stack([tr["%s_pIC50_direct_inhibition" % c].to_numpy(float) for c in CYPS], 1)

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    F = []
    for s in rows.SMILES:
        m = Chem.MolFromSmiles(s)
        if m is None:
            raise SystemExit("RDKit не разобрал %r" % s)
        F.append(gen.GetFingerprint(m))
    S = np.vstack([np.asarray(DataStructs.BulkTanimotoSimilarity(f, F), dtype=np.float32)
                   for f in F])
    assert np.allclose(S.max(1), 1.0), "КОНТРОЛЬ: максимум с собой не 1"
    iu = np.triu_indices(len(rows), k=1)
    sv = S[iu]
    sel = (sv >= BAND[0]) & (sv <= BAND[1])
    I, J = iu[0][sel], iu[1][sel]
    print("  все пары в полосе теста %.4f..%.4f: %d (из %d возможных)"
          % (*BAND, len(I), len(sv)))

    rng = np.random.default_rng(0)
    res, cond, out = {}, {}, {}
    print("\n  НАКЛОН ребёнка по родителю, обе стороны каждой пары")
    print("  %-8s %8s %8s %8s %9s %9s %12s"
          % ("фермент", "n пар", "наклон", "R^2", "ост. sd", "sd метки", "перемеш."))
    for e, c in enumerate(CYPS):
        ok = np.isfinite(Y[I, e]) & np.isfinite(Y[J, e])
        # both directions, so the fit is symmetric in the pair
        x = np.concatenate([Y[I[ok], e], Y[J[ok], e]])
        y = np.concatenate([Y[J[ok], e], Y[I[ok], e]])
        b, a = np.polyfit(x, y, 1)
        r2 = 1.0 - np.var(y - (a + b * x)) / np.var(y)
        sd = float(np.std(y - (a + b * x), ddof=2))
        bs = [float(np.polyfit(rng.permutation(x), y, 1)[0]) for _ in range(30)]
        res[c] = {"n_pairs": int(ok.sum()), "slope": float(b), "intercept": float(a),
                  "r2": float(r2), "resid_sd": sd, "label_sd": float(np.nanstd(Y[:, e])),
                  "shuffled_slope_mean": float(np.mean(bs)),
                  "shuffled_slope_sd": float(np.std(bs, ddof=1))}
        print("  %-8s %8d %8.3f %8.3f %9.3f %9.3f %12.4f"
              % (c, ok.sum(), b, r2, sd, res[c]["label_sd"], np.mean(bs)))
    bad = [c for c in CYPS if res[c]["slope"] <= 3 * res[c]["shuffled_slope_sd"]]
    print("  КОНТРОЛЬ перемешивания: наклон переживает его на %d ферментах из 4%s"
          % (4 - len(bad), (" (не пережил: %s)" % ", ".join(bad)) if bad else ""))

    print("\n  УСЛОВНОЕ распределение: пары, чей РОДИТЕЛЬ не слабее реальных родителей")
    print("  %-8s %9s %8s %9s %9s %10s %10s"
          % ("фермент", "порог род.", "n пар", "ср. реб.", "sd реб.", "ср. безусл.", "сдвиг"))
    for e, c in enumerate(CYPS):
        p = par["parent_%s" % c].to_numpy(float)
        thr = float(np.nanpercentile(p, 25))       # as active as the lower quartile of real parents
        ok = np.isfinite(Y[I, e]) & np.isfinite(Y[J, e])
        x = np.concatenate([Y[I[ok], e], Y[J[ok], e]])
        y = np.concatenate([Y[J[ok], e], Y[I[ok], e]])
        m = x >= thr
        uncond = float(np.nanmean(Y[:, e]))
        cond[c] = {"parent_threshold": thr, "n": int(m.sum()), "child_mean": float(y[m].mean()),
                   "child_sd": float(y[m].std(ddof=1)), "uncond_mean": uncond,
                   "shift": float(y[m].mean() - uncond)}
        print("  %-8s %9.3f %8d %9.3f %9.3f %10.3f %10.3f"
              % (c, thr, m.sum(), y[m].mean(), y[m].std(ddof=1), uncond, cond[c]["shift"]))
    flat = [c for c in CYPS if abs(cond[c]["shift"]) < 0.05]
    print("  КОНТРОЛЬ: условное отличается от безусловного на %d ферментах из 4%s"
          % (4 - len(flat), (" (почти не отличается: %s)" % ", ".join(flat)) if flat else ""))

    print("\n  ВЫВЕДЕННОЕ распределение теста против того, что подавали")
    print("  %-8s %10s %10s %10s %10s %10s"
          % ("фермент", "выведено", "sd вывед.", "best4", "stack374", "метка обуч."))
    sub = {k: pd.read_csv(RES + "submission/activity_submission_%s.csv" % k)
           for k in ("best4", "stack374")}
    for e, c in enumerate(CYPS):
        col = "%s_pIC50_direct_inhibition" % c
        out[c] = {"mean": cond[c]["child_mean"], "sd": cond[c]["child_sd"],
                  "best4": float(sub["best4"][col].mean()),
                  "stack374": float(sub["stack374"][col].mean()),
                  "train": float(np.nanmean(Y[:, e]))}
        print("  %-8s %10.3f %10.3f %10.3f %10.3f %10.3f"
              % (c, out[c]["mean"], out[c]["sd"], out[c]["best4"], out[c]["stack374"],
                 out[c]["train"]))
    json.dump({"band": BAND, "slope_fit": res, "conditional": cond, "derived": out},
              open(RES + "logs/k156_childlaw.json", "w"), ensure_ascii=False, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
