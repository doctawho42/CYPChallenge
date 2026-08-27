import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import pandas as pd, numpy as np
pd.set_option('display.width', 200)

tr = pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
te = pd.read_csv(D+"cyp-challenge-TEST-BLINDED.csv")
tdi= pd.read_csv(D+"cyp-challenge-TRAIN_TDI.csv")
sc = pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
em = pd.read_csv(D+"cyp-challenge-TRAIN_Emax.csv")
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]

print("=== TRAIN inhibition: n =", len(tr))
rows=[]
for c in CYPS:
    v=tr[f"{c}_pIC50_direct_inhibition"]; lo=tr[f"{c}_pIC50_direct_inhibition_conf_low"]; hi=tr[f"{c}_pIC50_direct_inhibition_conf_high"]
    w=hi-lo
    rows.append(dict(cyp=c, n=v.notna().sum(), nan=v.isna().sum(),
        min=round(v.min(),2), q25=round(v.quantile(.25),2), med=round(v.median(),2),
        q75=round(v.quantile(.75),2), max=round(v.max(),2), mean=round(v.mean(),2),
        frac_lt4=round((v<4).mean(),3), frac_lt45=round((v<4.5).mean(),3), frac_gt6=round((v>6).mean(),3),
        CIw_med=round(w.median(),3), CIw_med_act=round(w[v>=5].median(),3), CIw_med_inact=round(w[v<4].median(),3),
        CIw_q90=round(w.quantile(.9),2)))
print(pd.DataFrame(rows).to_string(index=False))

print("\n=== Correlations between isoforms (pIC50, pairwise complete) ===")
P = tr[[f"{c}_pIC50_direct_inhibition" for c in CYPS]]; P.columns=CYPS
print(P.corr(method='pearson').round(3).to_string())
print("\nSpearman:"); print(P.corr(method='spearman').round(3).to_string())
print("\nCo-occurrence (both measured):"); print(P.notna().astype(int).T.dot(P.notna().astype(int)).to_string())

print("\n=== TDI train: n =", len(tdi))
for c in ["CYP2D6","CYP3A4"]:
    y=tdi[f"{c}_is_TDI"]
    print(f"{c}: n={y.notna().sum()}, pos={int((y==True).sum())} ({(y==True).mean():.3f}), neg={int((y==False).sum())}, nan={y.isna().sum()}")
print("\nTDI vs direct pIC50 relation (CYP3A4):")
d=tdi[["CYP3A4_is_TDI","CYP3A4_pIC50_direct_inhibition","CYP3A4_pIC50_TDI_condition"]].dropna()
d["shift"]=d.CYP3A4_pIC50_TDI_condition-d.CYP3A4_pIC50_direct_inhibition
print(d.groupby("CYP3A4_is_TDI")[["CYP3A4_pIC50_direct_inhibition","CYP3A4_pIC50_TDI_condition","shift"]].describe().round(2).to_string())

print("\n=== single-conc: n =", len(sc))
print(sc.enzyme.value_counts().to_string())
print("compounds:", sc.Molecule_Name.nunique(), "| concentrations:", sc.concentration_M.unique()[:5])
print(sc[["log2fc_estimate","log2fc_std_error","cohens_d"]].describe().round(3).to_string())

print("\n=== Emax train: n =", len(em))
for c in CYPS:
    y=em[f"{c}_is_TDI"]
    print(f"  {c}_is_TDI: pos={int((y==True).sum())} ({(y==True).mean():.3f}) nan={y.isna().sum()}")
