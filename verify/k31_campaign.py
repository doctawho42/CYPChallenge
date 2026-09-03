"""Две популяции в обучающей выборке, и почему CYP3A4 --- это два набора, а не один.

Найдено при проверке предложения оценивать tau по скрину. Арифметика: скрин покрывает 4375 из
4905 строк, значит у похожей пары обе молекулы должны быть в нём в 0.892^2 = 79.6 % случаев.
На сходстве 0.35-0.45 так и есть (79.9 %), а на 0.60-0.85 --- 14.0 %. Аналоговые пары приходят
НЕ из диверсити-скрина.

Этот файл печатает состав: 530 молекул вне скрининговой библиотеки, все размечены ровно по
одному ферменту и это CYP3A4, и у двух третей из них есть близкий сосед против 4.9 % у
остальных. Подробности и следствия --- пункт 129.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / 'src'))
import numpy as np, pandas as pd
from cyppaths import D
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
RDLogger.DisableLog("rdApp.*")
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
rows=pd.read_csv(D+"rows.csv")
tr=(pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
      .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
sc=pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
inscr=rows.Molecule_Name.isin(set(sc.Molecule_Name)).to_numpy()
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
fp=[gen.GetFingerprint(s) for s in map(Chem.MolFromSmiles, rows.SMILES)]
n=len(fp)
nb=np.zeros(n)
for i in range(n):
    v=np.array(DataStructs.BulkTanimotoSimilarity(fp[i],fp)); v[i]=0; nb[i]=v.max()
close=nb>=0.587
print(f"{'группа':>16s} {'молекул':>8s} {'с соседом >=0.587':>19s}")
print(f"{'в скрине':>16s} {int(inscr.sum()):8d} {100*close[inscr].mean():18.1f} %")
print(f"{'ВНЕ скрина':>16s} {int((~inscr).sum()):8d} {100*close[~inscr].mean():18.1f} %")
print(f"\n{'фермент':8s} {'меток всего':>12s} {'из них вне скрина':>19s} {'из них с соседом':>18s} "
      f"{'вне скрина И с соседом':>23s}")
for c in CYPS:
    m=tr[f"{c}_pIC50_direct_inhibition"].notna().to_numpy()
    print(f"{c:8s} {int(m.sum()):12d} {int((m&~inscr).sum()):19d} {int((m&close).sum()):18d} "
          f"{int((m&~inscr&close).sum()):23d}")
print("""
Если «вне скрина» и «с близким соседом» --- почти одно множество, то серийная часть обучающей
выборки локализована, и якорный сплит строится ИМЕННО НА НЕЙ, а не на всей выборке. Это не
отменяет вывод о невозможности общей серийной валидации, но меняет его форму: не «серий нет»,
а «серии есть, их мало, и они лежат отдельной группой, у которой нет скрининга».""")
