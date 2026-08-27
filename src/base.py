import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import pandas as pd, numpy as np, sys
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, Descriptors
from cypsplit import butina_folds
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.metrics import mean_absolute_error, r2_score
from scipy.stats import spearmanr
RDLogger.DisableLog('rdApp.*')
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
mols=[Chem.MolFromSmiles(s) for s in tr.SMILES]
ok=[i for i,m in enumerate(mols) if m is not None]; tr=tr.loc[ok].reset_index(drop=True); mols=[mols[i] for i in ok]
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048,countSimulation=False)
FP=np.array([gen.GetCountFingerprintAsNumPy(m) for m in mols],dtype=np.float32)
dl=Descriptors.CalcMolDescriptors
desc=pd.DataFrame([dl(m) for m in mols]).replace([np.inf,-np.inf],np.nan)
desc=desc.loc[:,desc.isna().mean()<0.05].fillna(0).to_numpy(np.float32)
X=np.hstack([FP,desc]); print("X",X.shape)
folds,n_clusters=butina_folds(list(tr.SMILES))
print("clusters:",n_clusters)
res=[]
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; msk=tr[col].notna().to_numpy()
    y=tr.loc[msk,col].to_numpy(); lo=tr.loc[msk,col+"_conf_low"].to_numpy(); hi=tr.loc[msk,col+"_conf_high"].to_numpy()
    Xi=X[msk]; fi=folds[msk]; pred=np.zeros_like(y)
    for f in range(5):
        a,b=fi!=f, fi==f
        if b.sum()==0: continue
        mdl=HistGradientBoostingRegressor(max_iter=400,learning_rate=0.06,max_leaf_nodes=31,l2_regularization=1.0,random_state=0)
        mdl.fit(Xi[a],y[a]); pred[b]=mdl.predict(Xi[b])
    mean_pred=np.full_like(y,y.mean())
    res.append(dict(cyp=c,n=len(y),
        ST_RAE=round(strae(y,pred,y_true_upper=hi,y_true_lower=lo),4),
        ST_RAE_meanbaseline=round(strae(y,mean_pred,y_true_upper=hi,y_true_lower=lo),4),
        MAE=round(mean_absolute_error(y,pred),3), R2=round(r2_score(y,pred),3),
        rho=round(spearmanr(y,pred).statistic,3),
        frac_in_band=round(float(np.mean((pred>=lo)&(pred<=hi))),3),
        frac_true_band_gt1=round(float(np.mean((hi-lo)>1)),3)))
r=pd.DataFrame(res); r.loc[len(r)]=dict(cyp="MACRO",n=r.n.sum(),**{k:round(r[k].mean(),4) for k in r.columns if k not in("cyp","n")})
print(r.to_string(index=False))
