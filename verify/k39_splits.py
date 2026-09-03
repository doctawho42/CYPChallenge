"""Куда уходят сплиты при подвыборке признаков и без неё.

Гипотеза: HistGB смотрит все колонки в каждом узле, и 217 плотных дескрипторов плюс 30
механистических конкурируют с 2048 почти всюду нулевыми битами. Подвыборка колонок на узел
физически убирает часть битов из рассмотрения, и явные колонки начинают выбираться чаще.

Если так, часть выигрыша своего бустинга --- не регуляризация вообще, а высвобождение
механистического блока из-под фингерпринта. Тогда max_features и мех-блок не складываются,
а один включает другой.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / 'src'))
import numpy as np, pandas as pd
from cyppaths import D
from sklearn.tree import DecisionTreeRegressor
from cypsplit import butina_folds
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
NTREE,LR,DEPTH=200,0.06,5
rows=pd.read_csv(D+"rows.csv")
tr=(pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
      .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
z=np.load(D+"feats.npz"); nFP,nDE,nME=z["FP"].shape[1],z["DESC"].shape[1],z["MECH"].shape[1]
X=np.hstack([z["FP"],z["DESC"],z["MECH"]]).astype(np.float32)
mn=[l.strip() for l in open(D+"mech_names.csv")]
IDX_PROT=nFP+nDE+mn.index("frac_prot_74")
fold,_=butina_folds(list(rows.SMILES))

def block(i):
    return "FP" if i<nFP else ("DESC" if i<nFP+nDE else "MECH")

print(f"доля колонок в матрице: FP {100*nFP/(nFP+nDE+nME):.1f} %, "
      f"DESC {100*nDE/(nFP+nDE+nME):.1f} %, MECH {100*nME/(nFP+nDE+nME):.1f} %\n")
print(f"{'фермент':8s} {'mf':>5s} {'сплитов':>8s} {'FP':>7s} {'DESC':>7s} {'MECH':>7s} {'frac_prot_74':>13s}")
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; m=tr[col].notna().to_numpy()
    y=tr.loc[m,col].to_numpy(); Xi=X[m]; fi=fold[m]
    trn=fi!=0
    for mf in (1.0,0.3):
        rng=np.random.default_rng(0); s=np.full(trn.sum(),y[trn].mean()); cnt={}
        prot=0; tot=0
        for t in range(NTREE):
            tree=DecisionTreeRegressor(max_depth=DEPTH,max_features=mf,
                                       random_state=int(rng.integers(1<<30)))
            tree.fit(Xi[trn],y[trn]-s)
            s=s+LR*tree.predict(Xi[trn])
            f_=tree.tree_.feature; f_=f_[f_>=0]
            tot+=len(f_); prot+=int((f_==IDX_PROT).sum())
            for i in f_: cnt[block(i)]=cnt.get(block(i),0)+1
        print(f"{c:8s} {mf:5.1f} {tot:8d} {100*cnt.get('FP',0)/tot:6.1f} % "
              f"{100*cnt.get('DESC',0)/tot:6.1f} % {100*cnt.get('MECH',0)/tot:6.1f} % "
              f"{100*prot/tot:12.3f} %")
