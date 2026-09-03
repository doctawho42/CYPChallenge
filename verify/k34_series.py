"""Расслепление серийного слоя: сравнить разброс НАШИХ предсказаний внутри тестовой серии
с тем, сколько разброса химия вообще допускает на таком сходстве.

Меток не нужно. Три исхода, все информативны:
  дисперсия модели >> tau  --- стягивание оправдано, коэффициент вычисляется, а не подбирается;
  дисперсия модели ~  tau  --- слой ничего не изменит, ставку можно не делать;
  дисперсия модели << tau  --- модель УЖЕ переглажена внутри серий, и слой сделает хуже.

Предсказания берутся из отправляемого файла, то есть ПОСЛЕ аффинной пары: она сжимает выход
с коэффициентом lambda, и сравнивать надо то, что реально уйдёт, а не сырое.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / 'src'))
import numpy as np, pandas as pd
from cyppaths import D, RES
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
RDLogger.DisableLog("rdApp.*")
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
rows=pd.read_csv(D+"rows.csv")
tr=(pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
      .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
sub=pd.read_csv(RES+"submission/activity_submission.csv")
te=pd.read_csv(D+"cyp-challenge-TEST-BLINDED.csv")
sub=sub.set_index("Molecule_Name").loc[te.Molecule_Name].reset_index()
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
ftr=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]
fte=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in te.SMILES]
par=np.array([int(np.argmax(DataStructs.BulkTanimotoSimilarity(f,ftr))) for f in fte])
u,inv,cnt=np.unique(par,return_inverse=True,return_counts=True)
print(f"серий {len(u)}, из них с >=2 членами {int((cnt>=2).sum())}, "
      f"покрыто тестовых {int(cnt[cnt>=2].sum())} из {len(te)}\n")

# tau по обучающим парам того же сходства, в единицах pIC50
I,J=[],[]
n=len(ftr)
for i in range(n):
    v=np.array(DataStructs.BulkTanimotoSimilarity(ftr[i],ftr[i+1:]))
    for kk in np.where((v>=0.60)&(v<=0.85))[0]: I.append(i); J.append(i+1+kk)
I,J=np.array(I),np.array(J)

col=f"{CYPS[0]}_pIC50_direct_inhibition"
print(f"{'фермент':8s} {'пар':>5s} {'tau (химия)':>12s} {'sd модели':>11s} {'отношение':>10s} {'вывод':>22s}")
for c in CYPS:
    y=tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float); lab=~np.isnan(y)
    both=lab[I]&lab[J]
    tau=float(np.std(y[I[both]]-y[J[both]])/np.sqrt(2)) if both.sum()>3 else float('nan')
    p=sub[f"{c}_pIC50_direct_inhibition"].to_numpy(float)
    sds=[np.std(p[inv==k]) for k in range(len(u)) if cnt[k]>=2]
    sdm=float(np.mean(sds))
    r=sdm/tau if tau==tau and tau>0 else float('nan')
    verdict=("стягивание оправдано" if r>1.3 else
             "слой сделает хуже" if r<0.7 else "слой ничего не изменит")
    print(f"{c:8s} {int(both.sum()):5d} {tau:12.3f} {sdm:11.3f} {r:10.2f} {verdict:>22s}")
print("""
Как читать. tau --- внутрисерийный разброс, который химия допускает: sd разности меток в паре
делится на корень из двух. sd модели --- средний разброс наших предсказаний внутри серии, уже
после аффинной пары. Отношение ниже единицы означает, что модель разводит членов серии слабее,
чем они на самом деле различаются, и дополнительное стягивание только усилит эту ошибку.""")
