import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import pandas as pd, numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, Descriptors, Crippen
from rdkit import DataStructs
RDLogger.DisableLog('rdApp.*')

tr = pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
te = pd.read_csv(D+"cyp-challenge-TEST-BLINDED.csv")
sc = pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
print("test names sample:", te.Molecule_Name.head(3).tolist(), "| train names:", tr.Molecule_Name.head(3).tolist())
print("name overlap train/test:", len(set(tr.Molecule_Name)&set(te.Molecule_Name)))
print("single-conc names in train_inhib:", len(set(sc.Molecule_Name)&set(tr.Molecule_Name)), "of", sc.Molecule_Name.nunique())
print("single-conc names in test:", len(set(sc.Molecule_Name)&set(te.Molecule_Name)))

def canon(s):
    m=Chem.MolFromSmiles(s); return Chem.MolToSmiles(m) if m else None
tr["can"]=tr.SMILES.map(canon); te["can"]=te.SMILES.map(canon)
print("SMILES overlap train/test:", len(set(tr.can.dropna())&set(te.can.dropna())))

gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
def fps(smis):
    out=[]
    for s in smis:
        m=Chem.MolFromSmiles(s); out.append(gen.GetFingerprint(m) if m else None)
    return out
trf=fps(tr.SMILES); tef=fps(te.SMILES)
tr_ok=[i for i,f in enumerate(trf) if f is not None]
trf2=[trf[i] for i in tr_ok]

nn_idx=[]; nn_sim=[]
for f in tef:
    sims=np.array(DataStructs.BulkTanimotoSimilarity(f, trf2))
    j=int(sims.argmax()); nn_idx.append(tr_ok[j]); nn_sim.append(float(sims[j]))
te["nn_train_row"]=nn_idx; te["nn_sim"]=nn_sim
print("\n=== Test -> nearest train neighbour (ECFP4 Tanimoto) ===")
print(pd.Series(nn_sim).describe().round(3).to_string())
vc=pd.Series(nn_idx).value_counts()
print("distinct train anchors hit:", vc.nunique(), "| n unique anchors:", len(vc))
print("top anchor counts:", vc.head(15).tolist())
print("compounds whose NN is one of the top-%d anchors: %d" % (len(vc.head(80)), vc.head(80).sum()))

anchors=tr.loc[vc.index[:120]]
print("\n=== pIC50 profile of the most-hit train anchors (top 80 by count) ===")
a=tr.loc[vc.index[:80], [f"{c}_pIC50_direct_inhibition" for c in CYPS]]
a.columns=CYPS
print(a.describe().round(2).to_string())
print("\nHow do those anchors rank within train per CYP (percentile of their pIC50)?")
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"
    allv=tr[col].dropna()
    sub=tr.loc[vc.index[:80], col].dropna()
    if len(sub)==0: continue
    pct=[ (allv<v).mean() for v in sub ]
    print(f"  {c}: n={len(sub)}, median percentile={np.median(pct):.3f}, frac in top-25 of train={np.mean([v>=allv.nlargest(25).min() for v in sub]):.2f}")
te.to_csv(D+"test_with_nn.csv",index=False)
