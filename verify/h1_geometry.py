"""Где тест лежит на оси «случайное разбиение — кластерное».

Кластерное даёт медианное сходство с ближайшим обучающим 0.435, тест 0.587. Вопрос:
сколько даёт случайное разбиение по молекулам и сколько даёт потолок — соединение против
всей остальной обучающей выборки. Это ограничивает сверху то, что вообще достижимо
перенарезкой, и отвечает на вопрос §12: можно ли построить «псевдотестовый» фолд.

Ответ: нельзя. Потолок 0.450 ниже, чем сам тест (0.587); тест ближе к обучению, чем
обучение к самому себе. Как страту он строится, как разбиение — нет."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D
import numpy as np, pandas as pd
from rdkit import DataStructs
from cypsplit import cluster_ids, fingerprints

rows = pd.read_csv(D + "rows.csv"); te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
btr = fingerprints(rows.SMILES); bte = fingerprints(te.SMILES)
n = len(btr)
S = np.zeros((n, n), dtype=np.float32)
for i in range(n):
    S[i] = DataStructs.BulkTanimotoSimilarity(btr[i], btr)
np.fill_diagonal(S, -1.0)


def med_nn(fold):
    v = np.empty(n)
    for f in np.unique(fold):
        he = np.where(fold == f)[0]; tr = np.where(fold != f)[0]
        v[he] = S[np.ix_(he, tr)].max(1)
    return v


print(f"{'режим':44s} {'медиана':>8s} {'75%':>7s} {'доля>0.7':>9s}")
ceil = S.max(1)
print(f"{'потолок: против всей обучающей выборки':44s} {np.median(ceil):8.3f} "
      f"{np.percentile(ceil,75):7.3f} {np.mean(ceil>0.7):9.3f}")
for seed in range(4):
    v = med_nn(np.random.default_rng(seed).integers(0, 5, n))
    print(f"{'случайное 5-кратное, зерно ' + str(seed):44s} {np.median(v):8.3f} "
          f"{np.percentile(v,75):7.3f} {np.mean(v>0.7):9.3f}")
cid, n_clusters = cluster_ids(list(rows.SMILES))
fold = np.random.default_rng(0).integers(0, 5, n_clusters)[cid]
v = med_nn(fold)
print(f"{'кластерное 0.35, зерно 0':44s} {np.median(v):8.3f} {np.percentile(v,75):7.3f} "
      f"{np.mean(v>0.7):9.3f}")
nte = np.array([max(DataStructs.BulkTanimotoSimilarity(f, btr)) for f in bte])
print(f"{'ТЕСТ против всего обучения':44s} {np.median(nte):8.3f} {np.percentile(nte,75):7.3f} "
      f"{np.mean(nte>0.7):9.3f}")
sizes = np.bincount(cid)
print(f"\nразмер кластеров: синглтонов {int((sizes==1).sum())} из {n_clusters}"
      f" | соединений в кластерах >1: {int(n - (sizes==1).sum())}")
