"""Control for the sample-size confound: predict log2fc using only as many compounds
as that isoform has pIC50 labels. If rho stays high, the gap is about label RANGE,
not about how many labels there are."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor
from scipy.stats import spearmanr
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina
RDLogger.DisableLog('rdApp.*')
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
z=np.load(D+"feats.npz"); X=np.hstack([z["FP"],z["DESC"],z["MECH"]])
rows=pd.read_csv(D+"rows.csv")
tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv").set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index()
sc=pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
tr=tr.join(sc.pivot_table(index="Molecule_Name",columns="enzyme",values="log2fc_estimate"),on="Molecule_Name")
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
bits=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]
dists=[]
for i in range(1,len(bits)): dists.extend([1-x for x in DataStructs.BulkTanimotoSimilarity(bits[i],bits[:i])])
cl=Butina.ClusterData(dists,len(bits),0.35,isDistData=True)
cid=np.zeros(len(bits),int)
for k,c in enumerate(cl):
    for i in c: cid[i]=k
fold=np.random.default_rng(0).integers(0,5,len(cl))[cid]
rng=np.random.default_rng(7); res=[]
for c in CYPS:
    npic=int(tr[f"{c}_pIC50_direct_inhibition"].notna().sum())
    full=np.where(tr[c].notna().to_numpy())[0]
    sub=np.sort(rng.choice(full,npic,replace=False))
    msk=np.zeros(len(tr),bool); msk[sub]=True
    y=tr.loc[msk,c].to_numpy(); Xi,fi=X[msk],fold[msk]; p=np.zeros_like(y)
    for f in range(5):
        a,b=fi!=f,fi==f
        if b.sum()==0: continue
        p[b]=HistGradientBoostingRegressor(max_iter=300,learning_rate=0.06,random_state=0).fit(Xi[a],y[a]).predict(Xi[b])
    res.append(dict(cyp=c,n=npic,rho_log2fc_nmatched=round(spearmanr(p,y).statistic,3),
                    R2=round(1-((y-p)**2).sum()/((y-y.mean())**2).sum(),3),
                    iqr=round(float(np.percentile(y,75)-np.percentile(y,25)),2)))
print(pd.DataFrame(res).to_string(index=False))
