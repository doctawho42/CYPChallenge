import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import numpy as np, pandas as pd, json
from sklearn.ensemble import HistGradientBoostingClassifier
from cypsplit import butina_folds

z=np.load(D+"feats.npz"); X=np.hstack([z["FP"],z["DESC"],z["MECH"]])
rows=pd.read_csv(D+"rows.csv")
tdi=pd.read_csv(D+"cyp-challenge-TRAIN_TDI.csv").set_index("Molecule_Name")
keep=rows.Molecule_Name.isin(tdi.index).to_numpy()
T=tdi.loc[rows.Molecule_Name[keep]].reset_index(); Xt=X[keep]
fold,_=butina_folds(list(rows.SMILES[keep]))
out={}
for c in ["CYP3A4","CYP2D6"]:
    y=T[f"{c}_is_TDI"]; m=y.notna().to_numpy(); yv=y[m].astype(bool).to_numpy()
    Xi,fi=Xt[m],fold[m]; p=np.zeros(len(yv))
    for f in range(5):
        a,b=fi!=f,fi==f
        if b.sum()==0: continue
        p[b]=HistGradientBoostingClassifier(max_iter=300,learning_rate=0.06,random_state=0).fit(Xi[a],yv[a]).predict_proba(Xi[b])[:,1]
    out[c]={"y":yv.astype(int).tolist(),"p":p.tolist()}
    print(c,"done",len(yv),flush=True)
json.dump(out,open(RES+"preds/tdi_probs.json","w")); print("saved")
