import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import pandas as pd, numpy as np
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from scipy.cluster.hierarchy import linkage, fcluster
from scipy.spatial.distance import squareform
RDLogger.DisableLog('rdApp.*')

tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv"); te=pd.read_csv(D+"cyp-challenge-TEST-BLINDED.csv")
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
def fps(sm):
    o=[];
    for s in sm:
        m=Chem.MolFromSmiles(s); o.append(gen.GetFingerprint(m) if m else None)
    return o
tef=fps(te.SMILES); trf=fps(tr.SMILES)
n=len(tef); S=np.zeros((n,n))
for i in range(n):
    S[i,:]=DataStructs.BulkTanimotoSimilarity(tef[i],tef)
Dm=1-S; np.fill_diagonal(Dm,0); Dm=(Dm+Dm.T)/2
Z=linkage(squareform(Dm,checks=False),method='average')
for t in [0.45,0.5,0.55,0.6]:
    cl=fcluster(Z,t=t,criterion='distance')
    sizes=pd.Series(cl).value_counts()
    print(f"cut={t}: n_clusters={len(sizes)}, size dist: {sizes.value_counts().sort_index().to_dict()}")
cl=fcluster(Z,t=0.5,criterion='distance'); te["series"]=cl
sizes=pd.Series(cl).value_counts()
print(f"\nchosen cut 0.5 -> {len(sizes)} series; median size {sizes.median()}")

# anchor per series = train compound with max mean similarity to series members
tr_ok=[i for i,f in enumerate(trf) if f is not None]
trf2=[trf[i] for i in tr_ok]
rows=[]
for s in sorted(set(cl)):
    idx=np.where(cl==s)[0]
    M=np.zeros((len(idx),len(trf2)))
    for k,i in enumerate(idx): M[k,:]=DataStructs.BulkTanimotoSimilarity(tef[i],trf2)
    mean_sim=M.mean(0); j=int(mean_sim.argmax()); anchor=tr_ok[j]
    rows.append(dict(series=s,size=len(idx),anchor=anchor,anchor_meansim=round(float(mean_sim[j]),3),
        anchor_maxsim=round(float(M[:,j].max()),3),
        **{c: tr.loc[anchor,f"{c}_pIC50_direct_inhibition"] for c in CYPS}))
A=pd.DataFrame(rows)
print("\n=== Series anchors ===")
print("unique anchors:", A.anchor.nunique(), "of", len(A))
print("anchor mean-sim distribution:"); print(A.anchor_meansim.describe().round(3).to_string())
print("\nanchor has measured pIC50 for:"); print(A[CYPS].notna().sum().to_string())
print("\nanchor pIC50 medians:"); print(A[CYPS].median().round(2).to_string())

# within-series SAR spread: from TRAIN, |dpIC50| for pairs with Tanimoto in [0.5,0.7] and >=0.7
print("\n=== |ΔpIC50| vs similarity, measured on TRAIN pairs ===")
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"
    sub=tr[tr[col].notna()].reset_index(drop=True)
    f2=fps(sub.SMILES); ok=[i for i,x in enumerate(f2) if x is not None]
    sub=sub.loc[ok].reset_index(drop=True); f2=[f2[i] for i in ok]
    m=min(len(sub),900); rng=np.random.default_rng(0); pick=rng.choice(len(sub),m,replace=False)
    bins={"0.4-0.5":[], "0.5-0.6":[],"0.6-0.7":[],"0.7-0.85":[],"0.85-1":[]}
    for a in range(m):
        i=pick[a]; sims=np.array(DataStructs.BulkTanimotoSimilarity(f2[i],[f2[j] for j in pick]))
        d=np.abs(sub.loc[pick,col].to_numpy()-sub.loc[i,col])
        for k,(lo,hi) in zip(bins,[(0.4,0.5),(0.5,0.6),(0.6,0.7),(0.7,0.85),(0.85,1.0)]):
            msk=(sims>=lo)&(sims<hi); bins[k].extend(d[msk].tolist())
    print(f"  {c}: " + ", ".join(f"{k}: med={np.median(v):.2f} n={len(v)}" for k,v in bins.items() if len(v)>5))
te.to_csv(D+"test_series.csv",index=False); A.to_csv(D+"series_anchors.csv",index=False)
