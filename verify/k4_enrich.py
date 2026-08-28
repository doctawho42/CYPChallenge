"""Обогащён ли тестовый набор активными? Первая, грубая попытка: метка ближайшего соседа.

Меток теста у нас нет, но у каждого тестового соединения есть ближайший обучающий сосед,
и его метку можно взять за прокси. Контроль обязателен: та же процедура для обучающих
соединений против остальной обучающей выборки — иначе неясно, что считать нормой.

Прокси смещён и это видно прямо в выдаче: медианное сходство с ближайшим у теста выше,
чем у обучения с обучением (0.47-0.54 против 0.35-0.46), так что соседи теста
информативнее, а их распределение — не распределение случайной выборки. Направление
при этом устойчиво: шум тянет метку соседа К СРЕДНЕМУ, а не от него, поэтому сдвиг
вверх нельзя объяснить шумом, и настоящее обогащение скорее сильнее измеренного.

Аккуратная версия — k5_shift.py, где та же модель применяется к обоим наборам и
сравниваются распределения её собственных предсказаний."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D
import numpy as np, pandas as pd
from rdkit import DataStructs
from scipy.stats import mannwhitneyu, ks_2samp
from cypsplit import fingerprints

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
rows = pd.read_csv(D + "rows.csv"); te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
        .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
btr = fingerprints(rows.SMILES); bte = fingerprints(te.SMILES)
n = len(btr)
Str = np.zeros((n, n), dtype=np.float32)
for i in range(n):
    Str[i] = DataStructs.BulkTanimotoSimilarity(btr[i], btr)
np.fill_diagonal(Str, -1.0)
Ste = np.zeros((len(bte), n), dtype=np.float32)
for i in range(len(bte)):
    Ste[i] = DataStructs.BulkTanimotoSimilarity(bte[i], btr)

print(f"{'фермент':8s} {'n_обуч':>7s} | {'обуч. ybar':>10s} {'сосед обуч.':>11s} "
      f"{'сосед тест':>10s} | {'сдвиг':>7s} {'MWU p':>9s} {'KS p':>8s} | "
      f"{'сх. обуч':>8s} {'сх. тест':>8s}")
SH = {}
for c in CYPS:
    col = f"{c}_pIC50_direct_inhibition"; m = tr[col].notna().to_numpy()
    y = tr[col].to_numpy(); idx = np.where(m)[0]
    A = Str[np.ix_(idx, idx)]
    lab_tr = y[idx][A.argmax(1)]; sim_tr = A.max(1)
    B = Ste[:, idx]
    lab_te = y[idx][B.argmax(1)]; sim_te = B.max(1)
    u = mannwhitneyu(lab_te, lab_tr, alternative="two-sided")
    ks = ks_2samp(lab_te, lab_tr)
    SH[c] = (lab_tr, lab_te)
    print(f"{c:8s} {m.sum():7d} | {y[idx].mean():10.3f} {lab_tr.mean():11.3f} "
          f"{lab_te.mean():10.3f} | {lab_te.mean()-lab_tr.mean():+7.3f} {u.pvalue:9.2e} "
          f"{ks.pvalue:8.2e} | {np.median(sim_tr):8.3f} {np.median(sim_te):8.3f}")

print("\nДоли по зонам (метка ближайшего соседа), обучение -> тест:")
print(f"{'фермент':8s} {'<4.5':>14s} {'4.5-5.5':>14s} {'>5.5':>14s}")
for c in CYPS:
    a, b = SH[c]
    f = lambda v: (np.mean(v < 4.5), np.mean((v >= 4.5) & (v <= 5.5)), np.mean(v > 5.5))
    fa, fb = f(a), f(b)
    print(f"{c:8s} " + " ".join(f"{fa[k]*100:5.1f}->{fb[k]*100:5.1f}%" for k in range(3)))

print("\nКонтроль: то же на парах с высоким сходством (>0.7), где прокси надёжнее")
for c in CYPS:
    col = f"{c}_pIC50_direct_inhibition"; m = tr[col].notna().to_numpy()
    y = tr[col].to_numpy(); idx = np.where(m)[0]
    A = Str[np.ix_(idx, idx)]; B = Ste[:, idx]
    ma, mb = A.max(1) > 0.7, B.max(1) > 0.7
    la, lb = y[idx][A.argmax(1)], y[idx][B.argmax(1)]
    if ma.sum() < 5 or mb.sum() < 5:
        print(f"{c:8s} n={ma.sum():4d}/{mb.sum():3d}  пар слишком мало, не считаем")
        continue
    print(f"{c:8s} n={ma.sum():4d}/{mb.sum():3d}  {la[ma].mean():6.3f} -> {lb[mb].mean():6.3f} "
          f"({lb[mb].mean()-la[ma].mean():+.3f})  p={mannwhitneyu(lb[mb],la[ma]).pvalue:.3f}")
