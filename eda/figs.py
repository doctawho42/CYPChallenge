"""Все графики документа. Считаются из настоящих данных челленджа."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import numpy as np, pandas as pd, json, sys
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
F=RES + "fig/"
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
SHORT=["1A2","2C9","2D6","3A4"]
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
def save(fig,name):
    fig.savefig(F+name+".pdf",bbox_inches="tight",pad_inches=0.02); plt.close(fig); print("  ",name)

tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
sc=pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
piv=sc.pivot_table(index="Molecule_Name",columns="enzyme",values="log2fc_estimate")
trj=tr.join(piv,on="Molecule_Name")

# --- 1. заполненность обучающей матрицы -------------------------------------
fig,ax=plt.subplots(figsize=(5.4,2.0))
n=[int(tr[f"{c}_pIC50_direct_inhibition"].notna().sum()) for c in CYPS]
tot=len(tr)
ax.barh(SHORT,[tot]*4,color="#EDEFF6",height=.62)
ax.barh(SHORT,n,color=ACC,height=.62)
for i,(v,s) in enumerate(zip(n,SHORT)):
    ax.text(v+70,i,f"{v}  ({v/tot*100:.0f}%)",va="center",fontsize=8.4,color=INK)
ax.set_xlim(0,tot*1.28); ax.set_xlabel("соединений из 4905 в обучающей выборке")
ax.grid(axis="y",visible=False); ax.invert_yaxis()
save(fig,"fill")

# --- 2. корреляции изоформ ---------------------------------------------------
P=tr[[f"{c}_pIC50_direct_inhibition" for c in CYPS]]; P.columns=SHORT
C=P.corr(method="spearman").values
fig,ax=plt.subplots(figsize=(3.5,3.1))
im=ax.imshow(C,cmap="RdBu_r",vmin=-.8,vmax=.8)
ax.set_xticks(range(4),SHORT); ax.set_yticks(range(4),SHORT); ax.grid(False)
for i in range(4):
    for j in range(4):
        ax.text(j,i,f"{C[i,j]:.2f}",ha="center",va="center",fontsize=8.6,
                color="white" if abs(C[i,j])>.5 else INK)
cb=fig.colorbar(im,ax=ax,fraction=.046,pad=.04); cb.outline.set_linewidth(0)
cb.ax.tick_params(labelsize=7.5)
ax.set_title("ранговая корреляция активностей",pad=8)
save(fig,"corr")

# --- 3. распределения активностей -------------------------------------------
fig,axes=plt.subplots(1,4,figsize=(7.2,1.9),sharey=True)
for ax,c,s in zip(axes,CYPS,SHORT):
    v=tr[f"{c}_pIC50_direct_inhibition"].dropna()
    ax.hist(v,bins=26,range=(2,7.5),color=ACC,alpha=.82,edgecolor="none")
    q1,q3=v.quantile(.25),v.quantile(.75)
    ax.axvspan(q1,q3,color=HEME,alpha=.13,zorder=0)
    ax.set_title(f"{s}   размах {q3-q1:.2f}",fontsize=8.8)
    ax.set_xlabel("pIC$_{50}$"); ax.set_xlim(2,7.5); ax.grid(axis="x",visible=False)
axes[0].set_ylabel("соединений")
save(fig,"dists")

# --- 4. ширина интервала достоверности vs активность -------------------------
fig,ax=plt.subplots(figsize=(5.4,2.4))
for c,s,col in zip(CYPS,SHORT,[ACC,TEAL,HEME,GREY]):
    v=tr[f"{c}_pIC50_direct_inhibition"]; w=tr[f"{c}_pIC50_direct_inhibition_conf_high"]-tr[f"{c}_pIC50_direct_inhibition_conf_low"]
    d=pd.DataFrame({"v":v,"w":w}).dropna()
    b=pd.cut(d.v,bins=np.arange(2,7.6,.5))
    g=d.groupby(b,observed=True).w.median()
    x=[iv.mid for iv in g.index]
    ax.plot(x,g.values,marker="o",ms=3.4,lw=1.5,color=col,label=s)
ax.axvline(4,color=HEME,ls="--",lw=1,alpha=.7)
ax.text(4.05,2.55,"граница надёжности прибора",fontsize=7.6,color=HEME)
ax.set_xlabel("pIC$_{50}$"); ax.set_ylabel("медианная ширина\nинтервала, лог. ед.")
ax.legend(ncol=4,loc="upper right",fontsize=8)
save(fig,"ciwidth")

# --- 5. связь одноточечного скрининга и pIC50 --------------------------------
fig,axes=plt.subplots(1,4,figsize=(7.2,2.0),sharex=True)
bins=np.array([2,3,3.5,4,4.31,4.6,5,5.5,6,6.5,8])
for ax,c,s in zip(axes,CYPS,SHORT):
    col=f"{c}_pIC50_direct_inhibition"
    d=trj[[col,c]].dropna(); b=pd.cut(d[col],bins)
    g=d.groupby(b,observed=True)[c].median()
    x=[iv.mid for iv in g.index]
    ax.plot(x,g.values,marker="o",ms=3.6,lw=1.6,color=ACC)
    ax.axvline(4.31,color=TEAL,ls="--",lw=1)
    ax.set_title(s,fontsize=8.8); ax.set_xlabel("pIC$_{50}$"); ax.set_xlim(2,7.2)
    ax.grid(axis="x",visible=False)
axes[0].set_ylabel("медианный log2fc\nпри 49.5 мкМ")
axes[3].text(4.4,-3.6,"C$_0$",fontsize=8,color=TEAL)
save(fig,"hill")

# --- 6. базовая линия --------------------------------------------------------
st=[0.881,0.702,0.994,0.518]
fig,ax=plt.subplots(figsize=(4.4,2.2))
bars=ax.bar(SHORT,st,color=[HEME if v>.9 else ACC for v in st],width=.6)
ax.axhline(1.0,color=HEME,ls="--",lw=1.1)
ax.text(3.46,1.10,"уровень предсказания среднего",fontsize=7.6,color=HEME,ha="right")
ax.axhline(np.mean(st),color=TEAL,ls=":",lw=1.1)
ax.text(3.46,np.mean(st)+.025,f"среднее по треку {np.mean(st):.3f}",fontsize=7.6,color=TEAL,ha="right")
for b,v in zip(bars,st): ax.text(b.get_x()+b.get_width()/2,v/2,f"{v:.3f}",ha="center",va="center",fontsize=9,color="white",fontweight="bold")
ax.set_ylim(0,1.22); ax.set_ylabel("ST-RAE  (меньше — лучше)"); ax.grid(axis="x",visible=False)
save(fig,"baseline")

# --- 7. потолок надёжности ---------------------------------------------------
rel=[0.948,0.905,0.934,0.898]; r2=[0.216,0.350,0.122,0.589]
x=np.arange(4); w=.36
fig,ax=plt.subplots(figsize=(4.8,2.2))
ax.bar(x-w/2,rel,w,color="#D6DAF0",label="потолок (надёжность меток)")
ax.bar(x+w/2,r2,w,color=ACC,label="достигнутый R²")
for i,(a,b) in enumerate(zip(rel,r2)):
    ax.text(i+w/2,b+.02,f"{b/a*100:.0f}%",ha="center",fontsize=8,color=ACC)
ax.set_xticks(x,SHORT); ax.set_ylim(0,1.08); ax.set_ylabel("доля дисперсии")
ax.legend(loc="upper left",ncol=1); ax.grid(axis="x",visible=False)
save(fig,"ceiling")

# --- 8. ограничение диапазона ------------------------------------------------
iqr=[1.00,0.92,0.83,1.50]; gain=[0.027,0.075,0.128,-0.028]
fig,ax=plt.subplots(figsize=(4.6,2.6))
ax.axhline(0,color=LINE,lw=1)
ax.plot(iqr,gain,ls="--",lw=1,color=GREY,zorder=1)
for i,s in enumerate(SHORT):
    col=HEME if gain[i]<0 else ACC
    ax.scatter(iqr[i],gain[i],s=64,color=col,zorder=3)
    ax.annotate(s,(iqr[i],gain[i]),textcoords="offset points",xytext=(9,4),fontsize=8.6,color=col)
ax.set_xlabel("межквартильный размах меток pIC$_{50}$")
ax.set_ylabel("прирост $\\rho$ при переходе\nк широкому показателю")
ax.set_xlim(.75,1.62); ax.set_ylim(-.055,.16)
save(fig,"restriction")

# --- 9. сходство тест-трейн и обрывы активности ------------------------------
te=pd.read_csv(D+"test_with_nn.csv")
fig,ax=plt.subplots(figsize=(4.6,2.1))
ax.hist(te.nn_sim,bins=34,color=ACC,alpha=.85,edgecolor="none")
ax.axvline(te.nn_sim.median(),color=HEME,lw=1.2)
ax.text(te.nn_sim.median()+.008,ax.get_ylim()[1]*.86,f"медиана {te.nn_sim.median():.3f}",fontsize=8,color=HEME)
ax.set_xlabel("сходство тестового соединения с ближайшим обучающим (Танимото)")
ax.set_ylabel("соединений"); ax.grid(axis="x",visible=False)
save(fig,"nnsim")
print("figures done")
