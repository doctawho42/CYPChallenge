"""The similarity geometry of the three regimes, recomputed rather than recalled.

The question this answers: our validation holds out whole Butina clusters, which measures
EXTRAPOLATION to new scaffolds. The blind test was drawn by the organisers and may be much closer
to our training set than that. If so, every "this does not help" in the journal was measured on a
harder problem than the one the board scores.

Fingerprint and radius are cypsplit's own, so the numbers are commensurable with the split's 0.35
threshold. Self-similarity is excluded in the out-of-fold regimes and must be: a compound compared
against its own fold would score 1.000 against itself.

CONTROL THAT MUST SUCCEED and exercises the suspect clause: every training compound's maximum
similarity to the FULL training set INCLUDING itself must be exactly 1.000, which proves the matrix
and the argmax are right; the out-of-fold figures then differ only by the exclusion.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES
import numpy as np, pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from rdkit import DataStructs
from cypsplit import butina_folds

rows = pd.read_csv(D + "rows.csv")
te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)

def fps(smis):
    out = []
    for s in smis:
        m = Chem.MolFromSmiles(s)
        out.append(gen.GetFingerprint(m) if m is not None else None)
    bad = sum(1 for x in out if x is None)
    if bad:
        raise SystemExit("RDKit не разобрал %d структур" % bad)
    return out

FTR = fps(list(rows.SMILES))
FTE = fps(list(te.SMILES))
print("  обучение %d, тест %d, отпечаток Morgan r=2 2048 (как в cypsplit)" % (len(FTR), len(FTE)))

# full 4905 x 4905, as a list of rows to keep memory sane
S = np.vstack([np.asarray(DataStructs.BulkTanimotoSimilarity(f, FTR), dtype=np.float32)
               for f in FTR])
self_max = S.max(1)
assert np.allclose(self_max, 1.0), "КОНТРОЛЬ: не все максимумы с собой равны 1 (%.4f)" % self_max.min()
print("  КОНТРОЛЬ: максимум с собой = 1.000 для всех %d строк" % len(FTR))

STE = np.vstack([np.asarray(DataStructs.BulkTanimotoSimilarity(f, FTR), dtype=np.float32)
                 for f in FTE])

def report(name, nn):
    q = np.percentile(nn, [10, 25, 50, 75, 90])
    print("  %-28s медиана %.4f   >0.5 %5.1f%%  >0.6 %5.1f%%  >0.7 %5.1f%%  >0.8 %5.1f%%"
          % (name, np.median(nn), 100*(nn > .5).mean(), 100*(nn > .6).mean(),
             100*(nn > .7).mean(), 100*(nn > .8).mean()))
    return dict(median=float(np.median(nn)), q10=float(q[0]), q25=float(q[1]),
                q75=float(q[3]), q90=float(q[4]),
                **{"gt%d" % int(t*10): float((nn > t).mean()) for t in (.5, .6, .7, .8)})

out = {}
print("\n  БЛИЖАЙШИЙ СОСЕД В ОБУЧАЮЩЕМ МНОЖЕСТВЕ, по режимам")
out["тест -> всё обучение"] = report("тест -> всё обучение (4905)", STE.max(1))

fold, ncl = butina_folds(list(rows.SMILES))
nn_but = np.array([S[i][fold != fold[i]].max() for i in range(len(FTR))])
out["Butina 0.35, вне фолда"] = report("Butina 0.35, вне фолда", nn_but)

rng = np.random.default_rng(0)
rf = rng.integers(0, 5, len(FTR))
nn_rnd = np.array([S[i][rf != rf[i]].max() for i in range(len(FTR))])
out["случайные 5 фолдов"] = report("случайные 5 фолдов", nn_rnd)

# leave-one-out: the densest regime our own data can express
nn_loo = np.array([np.partition(S[i], -2)[-2] for i in range(len(FTR))])
out["исключая по одному"] = report("исключая по одному (потолок)", nn_loo)

print("\n  кластеров Butina: %d на %d соединений" % (ncl, len(FTR)))
print("  РАЗРЫВ: тест ближе к обучению, чем наш фолд, на %.4f медианы и в %.1f раза по доле >0.7"
      % (out["тест -> всё обучение"]["median"] - out["Butina 0.35, вне фолда"]["median"],
         out["тест -> всё обучение"]["gt7"] / max(out["Butina 0.35, вне фолда"]["gt7"], 1e-9)))
import json
json.dump(out, open(RES + "logs/k152_geometry.json", "w"), ensure_ascii=False, indent=1)
np.save(RES + "preds/k152_nn_test.npy", STE.max(1))
np.save(RES + "preds/k152_nn_butina.npy", nn_but)
