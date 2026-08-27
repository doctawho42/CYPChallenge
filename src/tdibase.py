"""Baseline for the TDI track + test of the 'regress both arms, then apply the rule' idea."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.metrics import matthews_corrcoef
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
RDLogger.DisableLog('rdApp.*')

z=np.load(D+"feats.npz"); X=np.hstack([z["FP"],z["DESC"],z["MECH"]])
rows=pd.read_csv(D+"rows.csv")
tdi=pd.read_csv(D+"cyp-challenge-TRAIN_TDI.csv").set_index("Molecule_Name")
idx=[n for n in rows.Molecule_Name if n in tdi.index]
keep=rows.Molecule_Name.isin(tdi.index).to_numpy()
T=tdi.loc[rows.Molecule_Name[keep]].reset_index(); Xt=X[keep]
print("соединений из TDI-файла, для которых есть признаки:", len(T), "из", len(tdi))
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
bits=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES[keep]]
dists=[]
for i in range(1,len(bits)): dists.extend([1-x for x in DataStructs.BulkTanimotoSimilarity(bits[i],bits[:i])])
cl=Butina.ClusterData(dists,len(bits),0.35,isDistData=True)
cid=np.zeros(len(bits),int)
for k,c in enumerate(cl):
    for i in c: cid[i]=k
fold=np.random.default_rng(0).integers(0,5,len(cl))[cid]
LOG2=np.log10(2)
for c in ["CYP3A4","CYP2D6"]:
    y=T[f"{c}_is_TDI"]; m=y.notna().to_numpy(); yv=y[m].astype(bool).to_numpy()
    Xi,fi=Xt[m],fold[m]
    print(f"\n=== {c}: n={m.sum()}, положительных {yv.mean()*100:.1f}% ===")
    # 1. прямая классификация
    p=np.zeros(len(yv))
    for f in range(5):
        a,b=fi!=f,fi==f
        if b.sum()==0: continue
        p[b]=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.06,random_state=0).fit(Xi[a],yv[a]).predict_proba(Xi[b])[:,1]
    best=max(((matthews_corrcoef(yv,p>t),t) for t in np.arange(0.05,0.95,0.01)))
    print(f"  прямая классификация: MCC при пороге 0.5 = {matthews_corrcoef(yv,p>0.5):.3f}; "
          f"лучший порог {best[1]:.2f} -> MCC {best[0]:.3f}")
    # 2. регрессия двух плеч -> правило
    dcol=f"{c}_pIC50_direct_inhibition"; tcol=f"{c}_pIC50_TDI_condition"
    md=T[dcol].notna().to_numpy()&m; mt=T[tcol].notna().to_numpy()&m
    pd_=np.full(len(T),np.nan); pt_=np.full(len(T),np.nan); psh=np.full(len(T),np.nan)
    both=md&mt
    for f in range(5):
        te=m&(fold==f)
        if te.sum()==0: continue
        a=md&(fold!=f); pd_[te]=HistGradientBoostingRegressor(max_iter=300,learning_rate=0.06,random_state=0).fit(X[keep][a],T.loc[a,dcol]).predict(X[keep][te])
        a=mt&(fold!=f); pt_[te]=HistGradientBoostingRegressor(max_iter=300,learning_rate=0.06,random_state=0).fit(X[keep][a],T.loc[a,tcol]).predict(X[keep][te])
        a=both&(fold!=f); psh[te]=HistGradientBoostingRegressor(max_iter=300,learning_rate=0.06,random_state=0).fit(X[keep][a],(T.loc[a,tcol]-T.loc[a,dcol])).predict(X[keep][te])
    di,td_,sh=pd_[m],pt_[m],psh[m]
    rule_sub=np.where(di>4,(td_-di)>LOG2,td_>4.301)
    rule_dif=np.where(di>4,sh>LOG2,td_>4.301)
    print(f"  правило по двум независимым регрессиям:  MCC {matthews_corrcoef(yv,rule_sub):.3f}")
    print(f"  правило, но сдвиг предсказан отдельной головой: MCC {matthews_corrcoef(yv,rule_dif):.3f}")
    bs=max(((matthews_corrcoef(yv,np.where(di>4,sh>t,td_>4.301)),t) for t in np.arange(0.0,1.0,0.02)))
    print(f"  то же с подобранным порогом сдвига {bs[1]:.2f}: MCC {bs[0]:.3f}")
    print(f"  всё-отрицательно:  MCC {matthews_corrcoef(yv,np.zeros(len(yv),bool)):.3f}")
