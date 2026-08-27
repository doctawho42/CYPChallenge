import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import pandas as pd, numpy as np, sys
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, Descriptors, Crippen
from sklearn.ensemble import HistGradientBoostingRegressor
from scipy.stats import spearmanr
RDLogger.DisableLog('rdApp.*')
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv"); sc=pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
mols=[Chem.MolFromSmiles(s) for s in tr.SMILES]
ok=[i for i,m in enumerate(mols) if m]; tr=tr.loc[ok].reset_index(drop=True); mols=[mols[i] for i in ok]
basicN=Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(N[#6]=[O,N,S]);!$(N[a]);!$(N[SX4](=O)=O);!$(N#*)]")
amidine=Chem.MolFromSmarts("[NX3][CX3]=[NX2]"); acid=Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
tetraz=Chem.MolFromSmarts("c1nnn[nH]1")
f=pd.DataFrame({
 "nBasicN":[len(m.GetSubstructMatches(basicN))+len(m.GetSubstructMatches(amidine)) for m in mols],
 "nAcid":[len(m.GetSubstructMatches(acid))+len(m.GetSubstructMatches(tetraz)) for m in mols],
 "logP":[Crippen.MolLogP(m) for m in mols],
 "MW":[Descriptors.MolWt(m) for m in mols],
 "nAr":[Descriptors.NumAromaticRings(m) for m in mols],
 "TPSA":[Descriptors.TPSA(m) for m in mols]})
tr=pd.concat([tr,f],axis=1)
print("=== pIC50 by basic-nitrogen count (median / n) ===")
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; d=tr[tr[col].notna()]
    g=d.groupby(d.nBasicN.clip(0,2))[col].agg(['median','count']).round(2)
    print(f"{c}: "+" | ".join(f"nB={i}: med={r['median']} (n={int(r['count'])})" for i,r in g.iterrows()))
print("\n=== pIC50 by acid group (median / n) ===")
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; d=tr[tr[col].notna()]
    g=d.groupby(d.nAcid.clip(0,1))[col].agg(['median','count']).round(2)
    print(f"{c}: "+" | ".join(f"nA={i}: med={r['median']} (n={int(r['count'])})" for i,r in g.iterrows()))
print("\n=== Spearman of simple physchem descriptors with pIC50 ===")
print(pd.DataFrame({c: {k: round(spearmanr(tr.loc[tr[f'{c}_pIC50_direct_inhibition'].notna(),k],
     tr.loc[tr[f'{c}_pIC50_direct_inhibition'].notna(),f'{c}_pIC50_direct_inhibition']).statistic,3)
     for k in ["nBasicN","nAcid","logP","MW","nAr","TPSA"]} for c in CYPS}).to_string())

# --- transfer test: model trained on log2fc only, evaluated on DRC pIC50 ranking ---
print("\n=== Transfer: GBM trained on single-conc log2fc -> ranks DRC pIC50? ===")
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
FP=np.array([gen.GetCountFingerprintAsNumPy(m) for m in mols],dtype=np.float32)
X=np.hstack([FP,f.to_numpy(np.float32)])
piv=sc.pivot_table(index="Molecule_Name",columns="enzyme",values="log2fc_estimate")
tr=tr.join(piv,on="Molecule_Name",rsuffix="_l2")
rng=np.random.default_rng(0); fold=rng.integers(0,5,len(tr))
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"
    hasy=tr[col].notna().to_numpy(); hasl=tr[c].notna().to_numpy()
    rho=[]
    for k in range(5):
        te=hasy&(fold==k); trn=hasl&(fold!=k)          # train on log2fc only, never on pIC50
        if te.sum()<30: continue
        m1=HistGradientBoostingRegressor(max_iter=300,learning_rate=0.06,random_state=0).fit(X[trn],tr.loc[trn,c])
        p=-m1.predict(X[te])
        rho.append(spearmanr(p,tr.loc[te,col]).statistic)
    print(f"  {c}: spearman(pred_from_log2fc_model, true pIC50) = {np.mean(rho):.3f}  (DRC-supervised baseline was in eda: see base.py)")
