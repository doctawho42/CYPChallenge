"""Why pooling's advantage grows exactly where the test set sits.

k22 measured the largest unexplained number in this log: pooling beats the per-enzyme model by
+0.0367 of rank in the layer where the nearest training neighbour is closer than 0.55, against
+0.0067 overall, and it is actively harmful in the far layer. The sign flip asks for a mechanism
and this file proposes one that is cheap to falsify.

The proposal. The label matrix is sparse -- CYP1A2 has 1412 labels out of 4905 rows, CYP2C9 has
1285 -- so a molecule's nearest structural neighbour usually carries a label for some *other*
enzyme. The per-enzyme model cannot see that row at all; it is not in its training set. The
pooled model can, through the enzyme indicator. So the pooled model's extra reach is exactly the
neighbours that are near but unlabelled for the enzyme being predicted, and that reach is worth
most where neighbours are near in the first place.

The count that suggested it. Pairs of *co-labelled* molecules above Tanimoto 0.55 number 58, 48,
62 and 996 across the four enzymes -- essentially none, except on CYP3A4 -- while twelve per cent
of rows have some neighbour that near. The neighbours exist and are mostly labelled elsewhere.

The measurement. For every row compute two nearest-neighbour similarities across the fold: one
over all training molecules, one over training molecules labelled for this enzyme. Their
difference is how much of the neighbourhood is invisible to the per-enzyme model. If the
proposal is right, pooling's advantage is concentrated in the rows where that difference is
large, and it should be near zero where the difference is zero.

The falsification is sharp and it is not guaranteed to pass: if pooling's advantage is flat in
the gap, then whatever pooling is doing it is not borrowing neighbours, and the sign flip needs
a different explanation.

Reads results/preds/oof_pool_all.json. Writes nothing. ~5 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES

import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from cypsplit import butina_folds

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
SEEDS = (0, 1, 2, 3)
GAPS = [(-0.01, 0.02), (0.02, 0.08), (0.08, 0.18), (0.18, 1.01)]
GLAB = ["~ 0", "0.02-0.08", "0.08-0.18", "> 0.18"]


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    P = json.load(open(RES + "preds/oof_pool_all.json"))["preds"]
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]

    share, acc = [], {}
    for seed in SEEDS:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for c in CYPS:
            col = f"{c}_pIC50_direct_inhibition"
            m = tr[col].notna().to_numpy()
            y = tr.loc[m, col].to_numpy()
            a = np.asarray(P[f"{seed}|независимо|{c}"], float)
            b = np.asarray(P[f"{seed}|пул|{c}"], float)
            idx = np.where(m)[0]

            nn_all = np.zeros(len(idx))
            nn_lab = np.zeros(len(idx))
            for f in range(5):
                te = np.where(fold[idx] == f)[0]
                if len(te) == 0:
                    continue
                pool_all = [fps[j] for j in np.where(fold != f)[0]]
                pool_lab = [fps[j] for j in idx[fold[idx] != f]]
                for t in te:
                    i = idx[t]
                    nn_all[t] = max(DataStructs.BulkTanimotoSimilarity(fps[i], pool_all))
                    nn_lab[t] = max(DataStructs.BulkTanimotoSimilarity(fps[i], pool_lab))
            gap = nn_all - nn_lab
            share.append((c, float(np.mean(gap > 0.02)), float(np.median(gap)),
                          float(np.median(nn_all)), float(np.median(nn_lab))))
            for k, (g0, g1) in enumerate(GAPS):
                s = (gap >= g0) & (gap < g1)
                if s.sum() < 30:
                    continue
                acc.setdefault(k, []).append(
                    (float(spearmanr(y[s], b[s]).statistic)
                     - float(spearmanr(y[s], a[s]).statistic), int(s.sum())))
        print(f"  сид {seed} готов", flush=True)

    print("\nСколько соседства не видно поферментной модели (сид 0):")
    print(f"{'фермент':8s} {'меток':>7s} {'медиана nn все':>15s} {'медиана nn с меткой':>20s} "
          f"{'медианный зазор':>16s} {'доля зазор>0.02':>17s}")
    for c in CYPS:
        r = [x for x in share if x[0] == c][0]
        n = int(tr[f"{c}_pIC50_direct_inhibition"].notna().sum())
        print(f"{c:8s} {n:7d} {r[3]:15.3f} {r[4]:20.3f} {r[2]:16.3f} {100*r[1]:16.1f} %")

    print("\nПреимущество пула по величине зазора (ранг пула минус ранг поферментной):")
    print(f"{'зазор':>12s} {'n':>7s} {'преимущество пула':>20s}")
    for k in sorted(acc):
        v = np.array([x[0] for x in acc[k]])
        n = sum(x[1] for x in acc[k]) // len(SEEDS)
        print(f"{GLAB[k]:>12s} {n:7d} {v.mean():+20.4f}")

    print("""
Как читать. Механизм предсказывает монотонный рост слева направо и величину около нуля в
первой строке: там, где вся окрестность уже размечена по этому ферменту, пулу нечего
одалживать, и он обязан не выигрывать.

Плоская таблица опровергает механизм целиком --- и это самый полезный исход из двух, потому
что тогда преимущество пула в ближнем слое остаётся необъяснённым, а необъяснённый эффект
такого размера опаснее объяснённого: мы не знаем, перенесётся ли он на тест.""")


if __name__ == "__main__":
    main()
