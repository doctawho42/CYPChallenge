"""Is the low effective sample size in item 113 the estimator's fault or the data's.

k25 reweighted the out-of-fold rows to the test set's similarity distribution and got
ESS/n between 20.7 and 39.3 per cent, low enough that two of four enzymes fall under the
one-third rule the file states for itself.  Before trying to fix that with a smoother
density ratio, it is worth knowing whether there is anything to fix.

The bound is the chi-square divergence between the two similarity distributions.  For
self-normalised importance weights,

    ESS / n  ->  1 / (1 + chi2(p_test || p_oof))

as n grows, whatever estimator produces the weights.  If the binned ratio already sits near
that bound, the weights are as good as weights can be here and the low ESS is a property of
how little the two distributions overlap -- which is, after all, the whole finding.  Only a
gap between the two numbers would be an argument for a better estimator.
"""
import sys as _sys
import pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
import numpy as np, pandas as pd
from cyppaths import D
from cypsplit import butina_folds
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
RDLogger.DisableLog("rdApp.*")

BINS = np.array([0.0, .25, .30, .35, .40, .45, .50, .55, .60, .70, 1.01])
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]

rows = pd.read_csv(D + "rows.csv")
tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
        .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]
fte = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in te.SMILES if Chem.MolFromSmiles(s)]
fold, _ = butina_folds(list(rows.SMILES))

s_oof = np.zeros(len(rows))
for f in range(5):
    t_, tn = np.where(fold == f)[0], np.where(fold != f)[0]
    ref = [fps[j] for j in tn]
    for i in t_:
        s_oof[i] = max(DataStructs.BulkTanimotoSimilarity(fps[i], ref))
rng = np.random.default_rng(0); n_sub = int(round(0.8 * len(fps)))
s_te = np.zeros(len(fte))
for _ in range(5):
    ref = [fps[j] for j in rng.choice(len(fps), n_sub, replace=False)]
    for t_, f_ in enumerate(fte):
        s_te[t_] += max(DataStructs.BulkTanimotoSimilarity(f_, ref)) / 5

h_o = np.histogram(s_oof, BINS)[0] / len(s_oof)
h_t = np.histogram(s_te, BINS)[0] / len(s_te)

# chi2(p_test || p_oof) = sum_b p_t^2 / p_o - 1, по тем же бинам.
ok = h_o > 0
chi2 = float((h_t[ok] ** 2 / h_o[ok]).sum() - 1.0)
print(f"chi2(тест || OOF) по бинам: {chi2:.3f}")
print(f"потолок ESS/n = 1/(1+chi2) = {100/(1+chi2):.1f} %")
lost = float(h_t[~ok].sum()) if (~ok).any() else 0.0
if lost:
    print(f"  (в бинах без опоры лежит {100*lost:.1f} % теста --- там вес не определён)")

for clip in (4.0, 8.0, 20.0, np.inf):
    ratio = np.where(ok, h_t / np.maximum(h_o, 1e-12), 1.0)
    ratio = np.minimum(ratio, clip)
    w_all = ratio[np.clip(np.digitize(s_oof, BINS) - 1, 0, len(ratio) - 1)]
    line = []
    for c in CYPS:
        m = tr[f"{c}_pIC50_direct_inhibition"].notna().to_numpy()
        w = w_all[m]
        line.append(100 * w.sum() ** 2 / (w * w).sum() / m.sum())
    tag = "без обрезки" if not np.isfinite(clip) else f"обрезка {clip:g}"
    print(f"  {tag:14s} ESS/n по ферментам: " + " ".join(f"{v:5.1f} %" for v in line))

print("""
Как читать. Если достигнутый ESS близок к потолку, оценщик ни при чём и низкий ESS ---
это и есть измеряемое явление: тест лежит там, где нашей проверки почти нет. Тогда
дисклеймер в пункте 113 остаётся, а искать более гладкий оценщик бессмысленно.""")
