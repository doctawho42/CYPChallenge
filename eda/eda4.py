import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import pandas as pd, numpy as np

tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv"); sc=pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
print(sc.head(3).to_string())
print("\nconcentration:", sc.concentration_M.unique(), " plates:", sc.plate_id.nunique())
piv=sc.pivot_table(index="Molecule_Name",columns="enzyme",values="log2fc_estimate")
m=tr.set_index("Molecule_Name").join(piv,how="left",rsuffix="_l2fc")
print("\n=== log2fc @49.5uM  vs  DRC availability / pIC50 ===")
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; l=c
    has=m[col].notna(); l2=m[l]
    print(f"\n{c}: compounds with l2fc={l2.notna().sum()}")
    print(f"   l2fc | DRC present : mean={l2[has].mean():.2f} med={l2[has].median():.2f}  n={l2[has].notna().sum()}")
    print(f"   l2fc | DRC absent  : mean={l2[~has].mean():.2f} med={l2[~has].median():.2f}  n={l2[~has].notna().sum()}")
    both=m[has & l2.notna()]
    if len(both)>10:
        r=np.corrcoef(both[l],both[col])[0,1]; rs=both[[l,col]].corr(method='spearman').iloc[0,1]
        print(f"   corr(l2fc, pIC50): pearson={r:.3f} spearman={rs:.3f}  n={len(both)}")
        # how well does l2fc alone separate active/inactive
        for thr in [-1.0,-1.5,-2.0]:
            sel=both[l]<thr
            print(f"     l2fc<{thr}: n={sel.sum():4d}  median pIC50={both.loc[sel,col].median():.2f} | l2fc>={thr}: n={(~sel).sum():4d} median pIC50={both.loc[~sel,col].median():.2f}")
# what fraction of DRC-measured compounds were "hits"
print("\n=== Selection function: P(DRC | log2fc) for CYP3A4 ===")
c="CYP3A4"; col=f"{c}_pIC50_direct_inhibition"
d=m[m[c].notna()].copy(); d["has"]=d[col].notna()
d["bin"]=pd.cut(d[c],bins=[-10,-3,-2,-1.5,-1,-0.5,-0.25,0,10])
print(d.groupby("bin",observed=True).agg(n=("has","size"),frac_DRC=("has","mean"),med_pIC50=(col,"median")).round(3).to_string())
