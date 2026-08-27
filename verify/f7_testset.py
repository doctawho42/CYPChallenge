"""Пересчёт утверждений про устройство тестового набора: 75 серий по 10, перцентили якорей."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
RDLogger.DisableLog('rdApp.*')
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
te=pd.read_csv(D+"cyp-challenge-TEST-BLINDED.csv"); tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
print("тест:",len(te),"строк,",te.SMILES.nunique(),"уникальных SMILES")
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
mte=[Chem.MolFromSmiles(s) for s in te.SMILES]; mtr=[Chem.MolFromSmiles(s) for s in tr.SMILES]
bte=[gen.GetFingerprint(m) for m in mte if m]; btr=[gen.GetFingerprint(m) for m in mtr if m]
oktr=[i for i,m in enumerate(mtr) if m]
d=[]
for i in range(1,len(bte)): d.extend([1-x for x in DataStructs.BulkTanimotoSimilarity(bte[i],bte[:i])])
for cut in [0.3,0.35,0.4,0.45]:
    cl=Butina.ClusterData(d,len(bte),cut,isDistData=True)
    sz=np.array([len(c) for c in cl])
    print(f"  порог {cut}: {len(cl)} групп | размер: медиана {np.median(sz):.0f}, "
          f"доля ровно по 10 = {np.mean(sz==10):.2f}, всего в группах >=5: {sz[sz>=5].sum()}")
# ближайший родственник каждой тестовой группы в обучении, и его перцентиль по активности
cl=Butina.ClusterData(d,len(bte),0.35,isDistData=True)
big=[c for c in cl if len(c)>=5]
print(f"\nгрупп размера >=5: {len(big)}; в них {sum(len(c) for c in big)} соединений")
anch=[]
for c in big:
    best=(-1,None)
    for i in c:
        s=DataStructs.BulkTanimotoSimilarity(bte[i],btr)
        j=int(np.argmax(s))
        if s[j]>best[0]: best=(s[j],oktr[j])
    anch.append(best)
print(f"медианное сходство якоря: {np.median([a[0] for a in anch]):.2f}")
idx=[a[1] for a in anch]
print("\nперцентиль якорей в обучающем распределении активности:")
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; full=tr[col].dropna().to_numpy()
    v=tr.loc[idx,col].dropna().to_numpy()
    pct=np.mean(full[None,:]<=v[:,None],axis=1)*100
    print(f"  {c}: якорей с меткой {len(v):3d}, медианный перцентиль {np.median(pct):5.1f}, средний {pct.mean():5.1f}")
