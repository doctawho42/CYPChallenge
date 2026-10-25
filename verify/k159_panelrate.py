"""Does an analogue of an active compound tend to be active? Measured on the NCGC panel, which is
the only unselected readout that covers CYP2C9 and CYP3A4.

Item 381 derived the test set's activity distribution for CYP1A2 and CYP2D6 and REFUSED it for
CYP2C9 and CYP3A4, because our single-concentration screen does not separate actives on those two
(Youden 0.139 and 0.151). Item 382 then closed the question against our own files: every other
readout we hold --- pIC50, Emax, the TDI curves --- exists only for compounds PROMOTED to
dose-response, and promotion was driven by the very activity we want to predict.

The NCGC panel is the one source without that hole, and `src/fetch1851.py` says why in its own
words: "one laboratory, one protocol, one campaign, the whole library screened, so there is no
publication filter to correct for". Inactives are reported rather than missing --- `censored` is
True for a compound that was tested and gave no fitted AC50 --- which is exactly what a rate over
actives and inactives needs and exactly what ChEMBL cannot give (items 363 and 364 measured
ChEMBL's selection reaching the ORDERING, not merely the level).

THIS IS A FIFTH ROLE FOR THE PANEL, not a reopening of item 356. That item measured it as a source
of TRAINING ROWS -- +0.0113 at member level, -0.0092 as a seventh member at sign 0/2 -- and closed
it there. Nothing is trained here. The panel is used as an instrument to measure one structural
rate, and the rate is a property of chemistry rather than of our model.

THE CONTROL THAT DECIDES WHETHER ANY OF THIS TRANSFERS, and it is the reason the two "easy"
isoforms are measured at all. On CYP1A2 and CYP2D6 item 381 already has an answer from our OWN
screen: a conditional rate of 0.521 against a base of 0.270, and 0.671 against 0.312, i.e. LIFTS of
1.93 and 2.15. A different laboratory, a different library and a different protocol must reproduce
those lifts before its answer on CYP2C9 and CYP3A4 may be believed. Absolute rates are not expected
to carry across assays -- sensitivity differs -- so the lift is the quantity compared.

TWO MORE CONTROLS:
  * self-similarity is exactly 1 on a sample, which proves the fingerprint matrix before any rate
    is read off it;
  * SHUFFLING the activity flags within an enzyme must drive the lift to 1, or the lift is an
    artefact of the pairing rather than a fact about analogues.

Writes results/logs/k159_panelrate.json.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES

import json

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from submit import CYPS

RDLogger.DisableLog("rdApp.*")
BAND = (0.5156, 0.7018)          # the test's own 10th-90th percentile of parent similarity (k155)
BLOCK = 512
BOOT = 400
# Item 381's own numbers are read from its log rather than copied here, so the two sources are
# compared on the SAME statistic and a change there cannot silently desync this file.
DERIVED_OK = ["CYP1A2", "CYP2D6"]        # the two item 381 could derive from our screen


def main():
    p = pd.read_csv(D + "ncgc/panel.csv", low_memory=False)
    p = p[p.enzyme.isin(CYPS)]
    act = (p.assign(active=~p.censored.astype(bool))
             .pivot_table(index="smiles", columns="enzyme", values="active", aggfunc="max"))
    smis = list(act.index)
    # pivot_table leaves object dtype (True/False/NaN); float makes it 1.0/0.0/NaN.
    A = np.stack([act[c].astype(float).to_numpy() for c in CYPS], 1)
    print("  панель: %d структур, по ферментам измерено %s"
          % (len(smis), "/".join(str(int(np.isfinite(A[:, e]).sum())) for e in range(4))))
    print("  доля активных: %s"
          % "  ".join("%s %.3f" % (c, np.nanmean(A[:, e])) for e, c in enumerate(CYPS)))

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    F, keep = [], []
    for i, s in enumerate(smis):
        m = Chem.MolFromSmiles(s)
        if m is not None:
            F.append(gen.GetFingerprint(m)); keep.append(i)
    A = A[keep]
    print("  отпечатки построены для %d из %d структур" % (len(F), len(smis)))
    chk = np.asarray(DataStructs.BulkTanimotoSimilarity(F[0], F[:200]))
    assert abs(chk[0] - 1.0) < 1e-9, "КОНТРОЛЬ: сходство с собой не 1"
    print("  КОНТРОЛЬ: сходство с собой = %.6f" % chk[0])

    I, J = [], []
    n = len(F)
    for a in range(0, n, BLOCK):
        b = min(a + BLOCK, n)
        S = np.vstack([np.asarray(DataStructs.BulkTanimotoSimilarity(F[i], F), dtype=np.float32)
                       for i in range(a, b)])
        S[:, :a] = -1.0                                   # upper triangle only
        for r in range(b - a):
            S[r, : a + r + 1] = -1.0
        ii, jj = np.where((S >= BAND[0]) & (S <= BAND[1]))
        I.append(ii + a); J.append(jj)
    I, J = np.concatenate(I), np.concatenate(J)
    print("  пар в полосе %.4f..%.4f: %d" % (*BAND, len(I)))

    rng = np.random.default_rng(0)
    out = {}
    print("\n  %-8s %9s %9s %11s %11s %7s %9s"
          % ("фермент", "n пар", "базов.пар", "аналог акт.", "аналог неакт.", "отнош.", "перемеш."))
    for e, c in enumerate(CYPS):
        ok = np.isfinite(A[I, e]) & np.isfinite(A[J, e])
        x = np.concatenate([A[I[ok], e], A[J[ok], e]]).astype(bool)
        y = np.concatenate([A[J[ok], e], A[I[ok], e]]).astype(bool)
        base = float(y.mean())                 # within the PAIRED population, not the library
        cond, anti = float(y[x].mean()), float(y[~x].mean())
        bs = []
        for _ in range(BOOT):
            k = rng.integers(0, len(x), len(x))
            if x[k].sum() > 5 and (~x[k]).sum() > 5:
                bs.append(y[k][x[k]].mean() / max(y[k][~x[k]].mean(), 1e-9))
        col = A[:, e].copy()
        fin = np.isfinite(col)
        col[fin] = rng.permutation(col[fin])
        xs = np.concatenate([col[I[ok]], col[J[ok]]]).astype(bool)
        ys = np.concatenate([col[J[ok]], col[I[ok]]]).astype(bool)
        sh = float(ys[xs].mean() / max(ys[~xs].mean(), 1e-9))
        out[c] = {"n_pairs": int(ok.sum()), "base": base, "cond": cond, "anti": anti,
                  "ratio": cond / max(anti, 1e-9),
                  "ratio_ci": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                  "ratio_shuffled": sh}
        print("  %-8s %9d %9.3f %11.3f %11.3f %7.2f %9.2f"
              % (c, ok.sum(), base, cond, anti, out[c]["ratio"], sh))

    k158 = json.load(open(RES + "logs/k158_derive.json"))["rates"]
    print("\n  КОНТРОЛЬ ПЕРЕНОСА: отношение акт./неакт., панель против НАШЕГО скрининга (п. 381)")
    print("  %-8s %13s %20s %12s" % ("фермент", "наш скрининг", "панель [95%]", "вердикт"))
    agree = []
    for c in CYPS:
        lo, hi = out[c]["ratio_ci"]
        if c in DERIVED_OK:
            r = k158[c]["analog_of_hit"] / max(k158[c]["analog_of_nonhit"], 1e-9)
            ok = lo <= r <= hi
            agree.append(ok)
            print("  %-8s %13.2f %8.2f [%.2f,%.2f] %12s"
                  % (c, r, out[c]["ratio"], lo, hi, "СОШЛОСЬ" if ok else "РАЗОШЛОСЬ"))
        else:
            print("  %-8s %13s %8.2f [%.2f,%.2f] %12s"
                  % (c, "невыведен", out[c]["ratio"], lo, hi, "— нечем сверить"))
    print("\n  перенос подтверждён на %d из 2 сверяемых ферментов" % sum(agree))
    bad = [c for c in CYPS if abs(out[c]["ratio_shuffled"] - 1.0) > 0.15]
    print("  КОНТРОЛЬ перемешивания (отношение должно уйти к 1): нарушено на %s"
          % (", ".join(bad) if bad else "ни одном"))
    json.dump({"band": BAND, "rates": out, "transfer_ok": sum(agree)},
              open(RES + "logs/k159_panelrate.json", "w"), ensure_ascii=False, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
