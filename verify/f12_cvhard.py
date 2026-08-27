"""Насколько наша кросс-валидация похожа по трудности на настоящий тест?
Сравниваем сходство «отложенное соединение -- ближайший сосед в обучающих фолдах»
со сходством «тестовое соединение -- ближайший сосед во всём обучении»."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
RDLogger.DisableLog('rdApp.*')

rows=pd.read_csv(D+"rows.csv"); te=pd.read_csv(D+"cyp-challenge-TEST-BLINDED.csv")
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
btr=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]
bte=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in te.SMILES]
d=[]
for i in range(1,len(btr)): d.extend([1-x for x in DataStructs.BulkTanimotoSimilarity(btr[i],btr[:i])])
cl=Butina.ClusterData(d,len(btr),0.35,isDistData=True)
cid=np.zeros(len(btr),int)
for k,c in enumerate(cl):
    for i in c: cid[i]=k
print("порог Бутины 0.35 ->",len(cl),"кластеров")
for seed in [0,1]:
    fold=np.random.default_rng(seed).integers(0,5,len(cl))[cid]
    nn=np.empty(len(btr))
    for i in range(len(btr)):
        other=np.where(fold!=fold[i])[0]
        s=DataStructs.BulkTanimotoSimilarity(btr[i],[btr[j] for j in other])
        nn[i]=max(s)
    q=np.percentile(nn,[50,75,90,95,99])
    print(f"  сид {seed}: сходство с ближайшим ИЗ ДРУГИХ фолдов -- медиана {q[0]:.3f}, 75% {q[1]:.3f}, "
          f"90% {q[2]:.3f}, 95% {q[3]:.3f}, 99% {q[4]:.3f}; доля >0.7: {np.mean(nn>0.7):.3f}")
    np.save(fRES + "nn_seed{seed}.npy",nn)
nte=np.array([max(DataStructs.BulkTanimotoSimilarity(f,btr)) for f in bte])
q=np.percentile(nte,[50,75,90,95,99])
print(f"\nТЕСТ против всего обучения -- медиана {q[0]:.3f}, 75% {q[1]:.3f}, 90% {q[2]:.3f}, "
      f"95% {q[3]:.3f}, 99% {q[4]:.3f}; доля >0.7: {np.mean(nte>0.7):.3f}")
nn=np.load(RES + "nn_seed0.npy")
print(f"\nразница медиан (кросс-валидация минус тест): {np.median(nn)-np.median(nte):+.3f}")
print("положительная разница = наша проверка ЛЕГЧЕ настоящего теста, числа оптимистичны")
