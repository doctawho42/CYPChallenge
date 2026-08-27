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
from sklearn.isotonic import IsotonicRegression
from scipy.stats import spearmanr
RDLogger.DisableLog('rdApp.*')
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv"); sc=pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
mols=[Chem.MolFromSmiles(s) for s in tr.SMILES]
ok=[i for i,m in enumerate(mols) if m]; tr=tr.loc[ok].reset_index(drop=True); mols=[mols[i] for i in ok]
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=1024)
FP=np.array([gen.GetCountFingerprintAsNumPy(m) for m in mols],dtype=np.float32)
bN=Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(N[#6]=[O,N,S]);!$(N[a]);!$(N[SX4](=O)=O);!$(N#*)]")
ac=Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
ph=np.array([[len(m.GetSubstructMatches(bN)),len(m.GetSubstructMatches(ac)),Crippen.MolLogP(m),
              Descriptors.MolWt(m),Descriptors.NumAromaticRings(m),Descriptors.TPSA(m),
              Descriptors.FractionCSP3(m),Descriptors.NumHDonors(m)] for m in mols],dtype=np.float32)
X=np.hstack([FP,ph])
bits=[gen.GetFingerprint(m) for m in mols]; dists=[]
for i in range(1,len(bits)): dists.extend([1-x for x in DataStructs.BulkTanimotoSimilarity(bits[i],bits[:i])])
cl=Butina.ClusterData(dists,len(bits),0.35,isDistData=True)
cid=np.zeros(len(bits),int)
for k,c in enumerate(cl):
    for i in c: cid[i]=k
rng=np.random.default_rng(0); fold=rng.integers(0,5,len(cl))[cid]
piv=sc.pivot_table(index="Molecule_Name",columns="enzyme",values="log2fc_estimate")
tr=tr.join(piv,on="Molecule_Name")
def gbm(): return HistGradientBoostingRegressor(max_iter=250,learning_rate=0.06,random_state=0)
rows=[]
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"
    hy=tr[col].notna().to_numpy(); hl=tr[c].notna().to_numpy()
    pA=np.full(len(tr),np.nan); pB=np.full(len(tr),np.nan); pE=np.full(len(tr),np.nan)
    for k in range(5):
        te=hy&(fold==k); a=hy&(fold!=k); b=hl&(fold!=k)
        if te.sum()<30: continue
        pA[te]=gbm().fit(X[a],tr.loc[a,col]).predict(X[te])                     # A: direct pIC50 regression
        mB=gbm().fit(X[b],tr.loc[b,c])                                          # B: predict log2fc from structure
        iso=IsotonicRegression(increasing=False,out_of_bounds='clip').fit(tr.loc[a&hl,c],tr.loc[a&hl,col])
        pB[te]=iso.predict(mB.predict(X[te]))                                   #    then invert via isotonic calibration
        pE[te]=0.5*pA[te]+0.5*pB[te]                                            # E: 50/50 blend
    y=tr.loc[hy,col].to_numpy(); lo=tr.loc[hy,col+"_conf_low"].to_numpy(); hi=tr.loc[hy,col+"_conf_high"].to_numpy()
    f=lambda p: round(strae(y,p[hy],y_true_upper=hi,y_true_lower=lo),4)
    g=lambda p: round(spearmanr(p[hy],y).statistic,3)
    rows.append(dict(cyp=c,n=int(hy.sum()),n_log2fc=int(hl.sum()),
        STRAE_A_pIC50direct=f(pA), STRAE_B_via_log2fc=f(pB), STRAE_blend=f(pE),
        rho_A=g(pA), rho_B=g(pB), rho_blend=g(pE)))
r=pd.DataFrame(rows); r.loc[len(r)]=dict(cyp="MACRO",n=r.n.sum(),n_log2fc=r.n_log2fc.sum(),
    **{k:round(r[k].mean(),4) for k in r.columns if k not in("cyp","n","n_log2fc")})
print(r.to_string(index=False))
