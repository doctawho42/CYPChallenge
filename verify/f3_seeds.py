"""Переживает ли вывод смену случайного разбиения на фолды?
Меняем ТОЛЬКО сид, которым кластеры раскидываются по фолдам. Кластеризация и признаки те же."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd, sys, time, json
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from sklearn.ensemble import HistGradientBoostingRegressor
from scipy.stats import spearmanr
from cypsplit import cluster_ids
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
z=np.load(D+"feats.npz"); FP,DESC,MECH=z["FP"],z["DESC"],z["MECH"]
rows=pd.read_csv(D+"rows.csv")
tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv").set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index()
cid,n_clusters=cluster_ids(list(rows.SMILES))
print("кластеров:",n_clusters,flush=True)
SETS={"FP+DESC":np.hstack([FP,DESC]),"FP+DESC+MECH":np.hstack([FP,DESC,MECH])}
res=[]; oofall={}
for seed in [1,2,3]:
    fold=np.random.default_rng(seed).integers(0,5,n_clusters)[cid]
    print(f"сид {seed}: размеры фолдов {np.bincount(fold)}",flush=True)
    for name,X in SETS.items():
        t0=time.time(); r={"seed":seed,"features":name}
        for c in CYPS:
            col=f"{c}_pIC50_direct_inhibition"; m=tr[col].notna().to_numpy()
            y=tr.loc[m,col].to_numpy(); lo=tr.loc[m,col+"_conf_low"].to_numpy(); hi=tr.loc[m,col+"_conf_high"].to_numpy()
            Xi,fi=X[m],fold[m]; p=np.zeros_like(y)
            for f in range(5):
                a,b=fi!=f,fi==f
                if b.sum()==0: continue
                mdl=HistGradientBoostingRegressor(max_iter=300,learning_rate=0.06,max_leaf_nodes=31,
                                                  l2_regularization=1.0,random_state=0)
                mdl.fit(Xi[a],y[a]); p[b]=mdl.predict(Xi[b])
            r[c]=round(strae(y,p,y_true_upper=hi,y_true_lower=lo),4)
            r["rho_"+c]=round(spearmanr(p,y).statistic,3)
            oofall[f"{seed}|{name}|{c}"]=p.tolist()
        r["MACRO"]=round(np.mean([r[c] for c in CYPS]),4)
        r["MACRO_rho"]=round(np.mean([r["rho_"+c] for c in CYPS]),3)
        res.append(r); print(f"  {name:14s} {time.time()-t0:.0f}s  макро {r['MACRO']:.4f}",flush=True)
        json.dump(oofall,open(RES+"preds/oof_seeds.json","w"))
        pd.DataFrame(res).to_csv(RES + "seeds.csv",index=False)
print(pd.DataFrame(res).to_string(index=False))
