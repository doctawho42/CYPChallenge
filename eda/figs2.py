import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import numpy as np, pandas as pd, json, sys
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.metrics import matthews_corrcoef
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import rdFingerprintGenerator
RDLogger.DisableLog('rdApp.*')
F=RES + "fig/"
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]; SHORT=["1A2","2C9","2D6","3A4"]
INK="#1A1F3A"; MUT="#5C6386"; LINE="#C8CEE0"
ACC="#2A32C8"; HEME="#9C3A25"; TEAL="#0E6E5C"; GREY="#8A90AC"
plt.rcParams.update({
 "font.family":"PT Sans","font.size":9,"axes.labelsize":9,"axes.titlesize":9.5,
 "xtick.labelsize":8.2,"ytick.labelsize":8.2,"legend.fontsize":8.2,
 "axes.edgecolor":LINE,"axes.labelcolor":INK,"text.color":INK,
 "xtick.color":MUT,"ytick.color":MUT,"axes.linewidth":0.8,
 "figure.facecolor":"white","axes.facecolor":"white",
 "axes.spines.top":False,"axes.spines.right":False,
 "axes.grid":True,"grid.color":"#E8EBF3","grid.linewidth":0.7,
 "legend.frameon":False,"pdf.fonttype":42,"axes.axisbelow":True})
def save(fig,n): fig.savefig(F+n+".pdf",bbox_inches="tight",pad_inches=0.02); plt.close(fig); print("  ",n)

# --- 10. абляция признаков ---------------------------------------------------
sets=["фингер-\nпринт","дескрип-\nторы","механисти-\nческий блок","ФП + деск-\nрипторы","всё\nвместе"]
ST={"1A2":[0.923,0.893,1.067,0.880,0.879],"2C9":[0.794,0.742,1.050,0.702,0.691],
    "2D6":[1.024,1.006,1.095,0.993,0.980],"3A4":[0.593,0.588,0.928,0.518,0.519]}
RHO={"1A2":[0.470,0.463,0.263,0.497,0.496],"2C9":[0.526,0.552,0.257,0.583,0.597],
     "2D6":[0.391,0.331,0.324,0.351,0.403],"3A4":[0.712,0.720,0.377,0.764,0.765]}
fig,axes=plt.subplots(1,2,figsize=(7.2,2.6))
x=np.arange(5); w=.2
for k,(s,col) in enumerate(zip(SHORT,[ACC,TEAL,HEME,GREY])):
    axes[0].bar(x+(k-1.5)*w,ST[s],w,color=col,label=s)
    axes[1].bar(x+(k-1.5)*w,RHO[s],w,color=col,label=s)
axes[0].axhline(1,color=HEME,ls="--",lw=1)
axes[0].set_ylabel("ST-RAE  (меньше — лучше)"); axes[0].set_ylim(0,1.18)
axes[1].set_ylabel(r"Спирмен $\rho$  (больше — лучше)"); axes[1].set_ylim(0,.85)
for ax in axes:
    ax.set_xticks(x,sets,fontsize=7.4); ax.grid(axis="x",visible=False)
axes[1].legend(ncol=4,loc="upper left")
save(fig,"ablation")

# --- 11. порог решения в TDI -------------------------------------------------
tp=json.load(open(RES+"preds/tdi_probs.json"))
fig,ax=plt.subplots(figsize=(4.8,2.5))
for c,s,col in [("CYP3A4","3A4",ACC),("CYP2D6","2D6",HEME)]:
    y=np.array(tp[c]["y"],bool); p=np.array(tp[c]["p"])
    th=np.arange(.05,.95,.01)
    m=[matthews_corrcoef(y,p>t) for t in th]
    ax.plot(th,m,lw=1.7,color=col,label=s)
    j=int(np.argmax(m)); ax.scatter([th[j]],[m[j]],s=48,color=col,zorder=4)
    ax.annotate(f"{m[j]:.3f} при {th[j]:.2f}",(th[j],m[j]),textcoords="offset points",
                xytext=(6,6),fontsize=8,color=col)
    k=int(np.argmin(np.abs(th-.5))); ax.scatter([th[k]],[m[k]],s=34,facecolor="white",
                edgecolor=col,zorder=4,linewidth=1.4)
ax.axvline(.5,color=GREY,ls=":",lw=1)
ax.text(.505,.02,"порог по умолчанию",fontsize=7.6,color=MUT)
ax.set_xlabel("порог принятия решения"); ax.set_ylabel("MCC")
ax.legend(loc="upper right"); ax.set_xlim(.05,.9)
save(fig,"threshold")

# --- 12. обрывы активности ---------------------------------------------------
tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
edges=[(0.3,0.4),(0.4,0.5),(0.5,0.6),(0.6,0.7),(0.7,1.0)]
res={}
rng=np.random.default_rng(0)
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"
    sub=tr[tr[col].notna()].reset_index(drop=True)
    fps=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in sub.SMILES]
    m=min(len(sub),800); pick=rng.choice(len(sub),m,replace=False)
    buckets=[[] for _ in edges]
    for a in range(m):
        i=pick[a]
        sims=np.array(DataStructs.BulkTanimotoSimilarity(fps[i],[fps[j] for j in pick]))
        d=np.abs(sub.loc[pick,col].to_numpy()-sub.loc[i,col])
        for k,(lo,hi) in enumerate(edges):
            msk=(sims>=lo)&(sims<hi); buckets[k].extend(d[msk].tolist())
    res[c]=[(np.median(b) if len(b)>8 else np.nan) for b in buckets]
    print("  cliffs",c,[f"{v:.2f}" if v==v else "-" for v in res[c]])
fig,ax=plt.subplots(figsize=(4.8,2.4))
xs=[np.mean(e) for e in edges]
for c,s,col in zip(CYPS,SHORT,[ACC,TEAL,HEME,GREY]):
    ax.plot(xs,res[c],marker="o",ms=4,lw=1.5,color=col,label=s)
ax.set_xlabel("сходство пары соединений (Танимото)")
ax.set_ylabel("медианная разница\npIC$_{50}$ в паре")
ax.legend(ncol=4,loc="upper right"); ax.set_ylim(0,1.35)
save(fig,"cliffs")
print("done")
