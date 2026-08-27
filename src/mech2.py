import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import pandas as pd, numpy as np, sys
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import rdFingerprintGenerator, Descriptors, Crippen
from rdkit.ML.Cluster import Butina
from sklearn.ensemble import HistGradientBoostingRegressor
from scipy.stats import spearmanr
RDLogger.DisableLog('rdApp.*')
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv"); sc=pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
mols=[Chem.MolFromSmiles(s) for s in tr.SMILES]
ok=[i for i,m in enumerate(mols) if m]; tr=tr.loc[ok].reset_index(drop=True); mols=[mols[i] for i in ok]
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=1024)
FP=np.array([gen.GetCountFingerprintAsNumPy(m) for m in mols],dtype=np.float32)
basicN=Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(N[#6]=[O,N,S]);!$(N[a]);!$(N[SX4](=O)=O);!$(N#*)]")
acid=Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
ph=np.array([[len(m.GetSubstructMatches(basicN)),len(m.GetSubstructMatches(acid)),Crippen.MolLogP(m),
              Descriptors.MolWt(m),Descriptors.NumAromaticRings(m),Descriptors.TPSA(m)] for m in mols],dtype=np.float32)
X=np.hstack([FP,ph])
bits=[gen.GetFingerprint(m) for m in mols]; dists=[]
for i in range(1,len(bits)):
    dists.extend([1-x for x in DataStructs.BulkTanimotoSimilarity(bits[i],bits[:i])])
cl=Butina.ClusterData(dists,len(bits),0.35,isDistData=True)
cid=np.zeros(len(bits),int)
for k,c in enumerate(cl):
    for i in c: cid[i]=k
rng=np.random.default_rng(0); fold=rng.integers(0,5,len(cl))[cid]
piv=sc.pivot_table(index="Molecule_Name",columns="enzyme",values="log2fc_estimate")
tr=tr.join(piv,on="Molecule_Name")
print("Matched Butina-cluster 5-fold split for BOTH models\n")
rows=[]
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"
    hy=tr[col].notna().to_numpy(); hl=tr[c].notna().to_numpy()
    pA=np.full(len(tr),np.nan); pB=np.full(len(tr),np.nan); pC=np.full(len(tr),np.nan)
    for k in range(5):
        te=hy&(fold==k)
        if te.sum()<30: continue
        a=hy&(fold!=k)
        pA[te]=HistGradientBoostingRegressor(max_iter=200,learning_rate=0.06,random_state=0).fit(X[a],tr.loc[a,col]).predict(X[te])
        b=hl&(fold!=k)
        mB=HistGradientBoostingRegressor(max_iter=200,learning_rate=0.06,random_state=0).fit(X[b],tr.loc[b,c])
        pB[te]=-mB.predict(X[te])
        aux=mB.predict(X).astype(np.float32).reshape(-1,1)
        X2=np.hstack([X,aux])
        pC[te]=HistGradientBoostingRegressor(max_iter=200,learning_rate=0.06,random_state=0).fit(X2[a],tr.loc[a,col]).predict(X2[te])
    y=tr.loc[hy,col].to_numpy(); lo=tr.loc[hy,col+"_conf_low"].to_numpy(); hi=tr.loc[hy,col+"_conf_high"].to_numpy()
    def sc_(p):
        p=p[hy]; return round(spearmanr(p,y).statistic,3)
    # rescale log2fc-model output to pIC50 units by isotonic-free linear fit on train folds (rank metric only for B)
    rows.append(dict(cyp=c, rho_DRConly=sc_(pA), rho_log2fc_only=sc_(pB), rho_stacked=sc_(pC),
                     STRAE_DRConly=round(strae(y,pA[hy],y_true_upper=hi,y_true_lower=lo),4),
                     STRAE_stacked=round(strae(y,pC[hy],y_true_upper=hi,y_true_lower=lo),4)))
r=pd.DataFrame(rows); r.loc[len(r)]=dict(cyp="MACRO",**{k:round(r[k].mean(),4) for k in r.columns if k!="cyp"})
print(r.to_string(index=False))
