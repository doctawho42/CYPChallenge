"""Every model's rank, measured in the regime the test set actually sits in.

The gap this closes. Item 97 measured that cross-validation runs at a nearest-neighbour
similarity of 0.435 while the test set sits at 0.587. Every architectural choice in the log
was therefore made in a harder regime than the one it will be scored in. That was found for
the Gaussian process, where the contribution collapses from -0.0161 in the far layer to
-0.0004 in the near one, and there is no reason it should be the only case.

So: recompute the rank of every saved model inside each similarity layer, and see whether
the ordering of the models is the same in the layer the test lives in as it is overall.

The immediate use is item 100's rejection of kNN. It scored rho 0.4816 against about 0.56
for everything else, and was closed on that. But a nearest-neighbour method is precisely the
method whose accuracy depends on how near the neighbours are, and it was judged at 0.435.
If it recovers at > 0.55 the rejection was delivered in the wrong regime.

What is compared. Rank inside a layer, not across layers: the layers have different label
spreads and their correlations are not comparable to each other. Comparing models within one
layer is sound, and that is the only comparison drawn.

Reads results/preds/oof_pool_all.json, oof_gp.json, oof_weak.json. Writes nothing.
Four seeds, about 4 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
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
LAYERS = [(0.0, 0.35), (0.35, 0.45), (0.45, 0.55), (0.55, 1.01)]
LAB = ["< 0.35", "0.35-0.45", "0.45-0.55", "> 0.55"]
SEEDS = (0, 1, 2, 3)

SRC = {"бустинг поферментно": ("oof_pool_all", "независимо"),
       "бустинг пул":         ("oof_pool_all", "пул"),
       "GP":                  ("oof_gp", "GP"),
       "лес":                 ("oof_weak", "лес"),
       "гребневая":           ("oof_weak", "гребневая"),
       "kNN":                 ("oof_weak", "kNN")}


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    cache = {f: json.load(open(RES + f"preds/{f}.json"))["preds"]
             for f in {v[0] for v in SRC.values()}}

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]

    acc = {}
    for seed in SEEDS:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        # Сходство с ближайшим обучающим ЧЕРЕЗ фолд: та же величина, что в пункте 97.
        nn = np.zeros(len(rows))
        for f in range(5):
            te, trn = np.where(fold == f)[0], np.where(fold != f)[0]
            ref = [fps[j] for j in trn]
            for i in te:
                nn[i] = max(DataStructs.BulkTanimotoSimilarity(fps[i], ref))

        for c in CYPS:
            col = f"{c}_pIC50_direct_inhibition"
            m = tr[col].notna().to_numpy()
            y = tr.loc[m, col].to_numpy()
            v = nn[m]
            for name, (fn, arm) in SRC.items():
                p = np.asarray(cache[fn][f"{seed}|{arm}|{c}"], float)
                for k, (a, b) in enumerate(LAYERS):
                    sel = (v >= a) & (v < b)
                    if sel.sum() < 30:
                        continue
                    acc.setdefault((name, k), []).append(
                        (float(spearmanr(y[sel], p[sel]).statistic), int(sel.sum())))
        print(f"  сид {seed} готов", flush=True)

    print()
    print(f"{'модель':22s}" + "".join(f"{l:>12s}" for l in LAB))
    print(f"{'n на слой (сумма)':22s}"
          + "".join(f"{sum(x[1] for x in acc[(list(SRC)[0], k)]) // len(SEEDS):12d}"
                    if (list(SRC)[0], k) in acc else f"{'-':>12s}" for k in range(4)))
    print("-" * 70)
    tab = {}
    for name in SRC:
        cells = []
        for k in range(4):
            if (name, k) not in acc:
                cells.append(np.nan); continue
            cells.append(float(np.mean([x[0] for x in acc[(name, k)]])))
        tab[name] = cells
        print(f"{name:22s}" + "".join(f"{x:12.4f}" for x in cells))

    print()
    print("Разрыв с лучшим бустингом внутри слоя (отрицательное --- модель хуже):")
    best = np.array(tab["бустинг поферментно"])
    for name in SRC:
        if name == "бустинг поферментно":
            continue
        d = np.array(tab[name]) - best
        print(f"  {name:22s}" + "".join(f"{x:+12.4f}" for x in d))

    print("""
Как читать. Сравнивать столбцы между собой нельзя --- у слоёв разный разброс меток, и
корреляции в них несопоставимы. Сравнивать модели ВНУТРИ столбца можно, и только это здесь
и делается.

Что решает. Тест лежит на медиане 0.587, то есть в правом столбце. Если порядок моделей в
правом столбце тот же, что и в целом --- вопрос о режиме закрыт, и все прежние отказы
вынесены законно. Если kNN в правом столбце подтягивается к бустингу --- пункт 100 отверг
его в чужом режиме, и разностная модель на ближних соседях становится осмысленной.""")


if __name__ == "__main__":
    main()
