# -*- coding: utf-8 -*-
"""Графики для общего документа. Один стиль на все, никаких украшений."""
import numpy as np, pandas as pd, json, sys, warnings
warnings.filterwarnings("ignore")
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from scipy.stats import spearmanr
from scipy.optimize import least_squares
from sklearn.isotonic import IsotonicRegression

D="/tmp/cyp/"; F="/tmp/doc/fig/"
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
# фиксированный порядок цветов: цвет закреплён за ферментом во всех графиках
COL={"CYP1A2":"#3260C4","CYP2C9":"#C0432A","CYP2D6":"#118A6C","CYP3A4":"#B98A10"}
INK="#1a1a1a"; MUT="#6b6b6b"; GRID="#d8d8d8"; SURF="#ffffff"

plt.rcParams.update({
    "font.family":"PT Sans", "font.size":9,
    "axes.edgecolor":"#9a9a9a", "axes.linewidth":0.7, "axes.labelcolor":INK,
    "axes.titlesize":9.5, "axes.titleweight":"bold", "axes.titlelocation":"left",
    "axes.titlepad":6, "axes.labelsize":9, "axes.facecolor":SURF, "axes.axisbelow":True,
    "axes.spines.top":False, "axes.spines.right":False,
    "xtick.color":MUT,"ytick.color":MUT,"xtick.labelsize":8.2,"ytick.labelsize":8.2,
    "xtick.major.width":0.7,"ytick.major.width":0.7,"xtick.major.size":3,"ytick.major.size":3,
    "grid.color":GRID,"grid.linewidth":0.55,"legend.frameon":False,"legend.fontsize":8.4,
    "figure.facecolor":SURF,"savefig.facecolor":SURF,"savefig.bbox":"tight","savefig.pad_inches":0.02,
    "lines.linewidth":1.7, "pdf.fonttype":42,
})
def save(fig,name):
    fig.savefig(F+name+".pdf"); fig.savefig(F+name+".png",dpi=185); plt.close(fig)
    print("  ->",name)

inh=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
sc=pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
piv=sc.pivot_table(index="Molecule_Name",columns="enzyme",values="log2fc_estimate")
mj=inh.set_index("Molecule_Name").join(piv)
oof=json.load(open(D+"oof.json")); rows=pd.read_csv(D+"rows.csv")
tr=inh.set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index()
pC0=-np.log10(4.95049505e-05)

# ------------------------------------------------------------------ 1. распределения
print("1 распределения")
fig,axes=plt.subplots(2,2,figsize=(6.10,3.75),sharey=True,sharex=True)
for ax,c in zip(axes.ravel(),CYPS):
    y=inh[f"{c}_pIC50_direct_inhibition"].dropna().to_numpy()
    q1,q3=np.percentile(y,[25,75])
    ax.hist(y,bins=34,color="#e3e3e3",edgecolor="#c2c2c2",linewidth=.4)
    ax.axvspan(q1,q3,color=COL[c],alpha=.22,lw=0)
    ax.axvline(np.median(y),color=COL[c],lw=1.8)
    ax.set_title(c,color=COL[c])
    ax.annotate(f"IQR {q3-q1:.2f}\nn = {len(y)}",xy=(.97,.93),xycoords="axes fraction",ha="right",va="top",
                fontsize=8.4,color=INK,linespacing=1.45)
    ax.set_xlim(2.4,7.6); ax.xaxis.set_major_locator(MultipleLocator(1))
    ax.grid(axis="y")
axes[0,0].set_ylim(0,300)
for ax in axes[:,0]: ax.set_ylabel("соединений")
for ax in axes[1,:]: ax.set_xlabel("pIC$_{50}$")
fig.subplots_adjust(hspace=.42,wspace=.10)
save(fig,"dists")

# ------------------------------------------------------------------ 2. калибровка
print("2 калибровка")
fig,axes=plt.subplots(2,2,figsize=(6.10,4.35),sharey=True,sharex=True)
CAL={}
for ax,c in zip(axes.ravel(),CYPS):
    col=f"{c}_pIC50_direct_inhibition"; k=mj[col].notna()&mj[c].notna()
    pi=mj.loc[k,col].to_numpy(); y=mj.loc[k,c].to_numpy()
    def resid(t_):
        E,h=t_; I=E/(1+10**(h*(pC0-pi)))
        return np.log2(np.clip(1-I,1e-3,None))-y
    E,h=least_squares(resid,[1.,1.],bounds=([.2,.2],[1.2,4.])).x
    CAL[c]=(E,h)
    ax.scatter(pi,y,s=3.0,color=COL[c],alpha=.28,lw=0,rasterized=True)
    g=np.linspace(2.0,7.9,300)
    ax.plot(g,np.log2(np.clip(1-E/(1+10**(h*(pC0-g))),1e-3,None)),color=INK,lw=1.7)
    iso=IsotonicRegression(increasing=False,out_of_bounds="clip").fit(pi,y)
    ax.plot(g,iso.predict(g),color=INK,lw=1.1,ls=(0,(3.5,2.5)))
    ax.axvline(pC0,color=MUT,lw=.8,ls=":")
    ax.set_title(c,color=COL[c])
    ax.annotate(f"$E$ = {E:.2f}   $h$ = {h:.2f}\n$\\rho$ = {spearmanr(pi,y).statistic:.2f}",
                xy=(.97,.93),xycoords="axes fraction",ha="right",va="top",fontsize=8.2,color=INK,linespacing=1.45)
    ax.grid(axis="y"); ax.set_xlim(1.9,8.0); ax.set_ylim(-5.4,3.2)
for ax in axes[:,0]: ax.set_ylabel("log$_2$fc скрининга")
for ax in axes[1,:]: ax.set_xlabel("pIC$_{50}$")
fig.subplots_adjust(hspace=.34,wspace=.09)
from matplotlib.lines import Line2D as _L
fig.legend(handles=[_L([],[],color=INK,lw=1.7,label="калибровка кривой Хилла"),
                    _L([],[],color=INK,lw=1.1,ls=(0,(3.5,2.5)),label="свободная монотонная (изотоника)"),
                    _L([],[],color=MUT,lw=.9,ls=":",label="концентрация скрининга p$C_0$")],
           loc="lower center",ncols=2,fontsize=8.2,bbox_to_anchor=(.5,-.115))
save(fig,"calib")

# ------------------------------------------------------------------ 3. информация
print("3 информация")
fig,ax=plt.subplots(figsize=(4.75,2.55))
pi=np.linspace(2,8,900)
for hh,ls,lab in [(1.0,"-","$h$ = 1.0"),(2.0,(0,(4,2)),"$h$ = 2.0"),(0.6,(0,(1.4,1.6)),"$h$ = 0.6")]:
    x=10**(hh*(pC0-pi)); ax.plot(pi,x/(1+x)**2,color=INK if hh==1 else MUT,ls=ls,lw=1.7 if hh==1 else 1.2,label=lab)
ax.axvline(pC0,color="#C0432A",lw=1.0)
ax.annotate(f"p$C_0$ = {pC0:.2f}\n49.5 мкМ",xy=(pC0+1.45,.175),fontsize=8.2,color="#C0432A",linespacing=1.3)
ax.axhline(.25,color=GRID,lw=.7,zorder=0)
ax.annotate("1/4",xy=(2.05,.253),fontsize=8,color=MUT)
ax.set_xlabel("истинная pIC$_{50}$ соединения"); ax.set_ylabel("$x/(1+x)^2$")
ax.set_xlim(2,8); ax.set_ylim(0,.285); ax.legend(loc="upper right"); ax.grid(axis="y")
save(fig,"fisher")

# ------------------------------------------------------------------ 4. Emax
print("4 Emax")
em=pd.read_csv(D+"cyp-challenge-TRAIN_Emax.csv")
fig,axes=plt.subplots(1,2,figsize=(6.10,2.75),gridspec_kw={"width_ratios":[1.45,1],"wspace":.10})
ax=axes[0]
for i,c in enumerate(CYPS):
    v=em[f"{c}_EmaxVsPosCtrl_direct_inhibition"].dropna().to_numpy()
    lo=em[f"{c}_EmaxVsPosCtrl_direct_inhibition_conf_low"]; hi=em[f"{c}_EmaxVsPosCtrl_direct_inhibition_conf_high"]
    w=np.nanmean((hi-lo).to_numpy())
    q5,q95=np.percentile(v,[5,95]); med=np.median(v)
    ax.plot([q5,q95],[i-.13,i-.13],color=COL[c],lw=7,solid_capstyle="butt",alpha=.9)
    ax.plot([med-w/2,med+w/2],[i+.22,i+.22],color="#8a8a8a",lw=7,solid_capstyle="butt")
    ax.plot(med,i-.13,"o",color=SURF,ms=4,mec=COL[c],mew=1.6,zorder=5)
    ax.annotate(f"{q95-q5:.2f}",xy=(q95+.015,i-.13),va="center",fontsize=8.2,color=COL[c],fontweight="bold")
    ax.annotate(f"{w:.2f}",xy=(med+w/2+.015,i+.22),va="center",fontsize=8.2,color="#6b6b6b")
ax.set_yticks(range(4)); ax.set_yticklabels(CYPS); ax.invert_yaxis()
ax.set_ylim(3.75,-1.35)
ax.set_xlabel("$E$ относительно положительного контроля")
ax.annotate("цветом — где лежат 90 % соединений\nсерым — ширина полосы одного измерения",
            xy=(.015,.985),xycoords="axes fraction",va="top",fontsize=8.2,color=INK,linespacing=1.4)
ax.set_xlim(-1.30,-0.68); ax.grid(axis="x")
ax=axes[1]
rel=[]
for c in CYPS:
    v=em[f"{c}_EmaxVsPosCtrl_direct_inhibition"]; lo=em[f"{c}_EmaxVsPosCtrl_direct_inhibition_conf_low"]; hi=em[f"{c}_EmaxVsPosCtrl_direct_inhibition_conf_high"]
    k=v.notna()&lo.notna()&hi.notna(); s=((hi[k]-lo[k])/(2*1.96)).to_numpy()
    rel.append(1-np.mean(s**2)/v[k].to_numpy().var(ddof=1))
ax.barh(range(4),rel,color=[COL[c] for c in CYPS],height=.5)
ax.axvline(0,color=INK,lw=1.0)
for i,r in enumerate(rel):
    ax.annotate(f"{r:+.2f}",xy=(r+.06,i),va="center",ha="left",fontsize=8.4,color=SURF,fontweight="bold")
ax.set_yticks(range(4)); ax.set_yticklabels([]); ax.invert_yaxis(); ax.set_ylim(3.75,-1.35)
ax.tick_params(axis="y",length=0)
ax.set_xlabel("надёжность $E$"); ax.set_xlim(-1.72,.30); ax.grid(axis="x")
ax.annotate("ноль означает, что весь наблюдаемый\nразброс объясняется шумом; ниже нуля —\nчто шума заявлено даже больше",
            xy=(.02,.985),xycoords="axes fraction",va="top",fontsize=8.2,color=INK,linespacing=1.4)
save(fig,"emax")

# ------------------------------------------------------------------ 5. ограничение диапазона
print("5 ограничение диапазона")
RG={"CYP2D6":(0.399,0.527,0.83),"CYP2C9":(0.599,0.674,0.92),"CYP1A2":(0.507,0.534,1.00),"CYP3A4":(0.763,0.735,1.50)}
fig,ax=plt.subplots(figsize=(4.95,3.05))
for c,(a,b,iqr) in RG.items():
    ax.annotate("",xy=(iqr,b),xytext=(iqr,a),
                arrowprops=dict(arrowstyle="-|>,head_width=.16,head_length=.34",color=COL[c],lw=1.7))
    ax.plot(iqr,a,"o",color=SURF,ms=5.5,mec=COL[c],mew=1.5,zorder=4)
    ax.plot(iqr,b,"o",color=COL[c],ms=5.5,zorder=4)
    dy=.018 if b>a else -.030
    ax.annotate(f"{c}\n{b-a:+.3f}",xy=(iqr,max(a,b)+.022),ha="center",fontsize=8.2,
                color=COL[c],fontweight="bold",linespacing=1.3)
ax.set_xlabel("межквартильный размах меток pIC$_{50}$")
ax.set_ylabel("Спирмен $\\rho$ на кросс-валидации")
ax.set_xlim(.72,1.62); ax.set_ylim(.36,.86); ax.grid(axis="y")
ax.plot([],[],"o",color=SURF,mec=MUT,mew=1.5,label="метки кривых (узкий диапазон)")
ax.plot([],[],"o",color=MUT,label="скрининг, столько же примеров")
ax.legend(loc="lower right")
save(fig,"range")

# ------------------------------------------------------------------ 6. сиды
print("6 сиды")
from matplotlib.patches import Patch
S=pd.read_csv("/tmp/verify/seeds_all.csv")
fig,axes=plt.subplots(1,2,figsize=(6.10,2.65),gridspec_kw={"width_ratios":[1,1.15],"wspace":.30})
ax=axes[0]
P=S.pivot(index="seed",columns="features",values="MACRO")
xs=np.arange(4)
for i in xs:
    ax.plot([i,i],[P["FP+DESC"].iloc[i],P["FP+DESC+MECH"].iloc[i]],color="#dcdcdc",lw=1.2,zorder=0)
ax.plot(xs,P["FP+DESC"],"o-",color="#8a8a8a",ms=6,zorder=3)
ax.plot(xs,P["FP+DESC+MECH"],"s-",color="#118A6C",ms=5.5,zorder=3)
ax.annotate("без блока",xy=(0.06,P["FP+DESC"].iloc[0]+.0011),fontsize=8.4,color="#6b6b6b",fontweight="bold")
ax.annotate("с блоком",xy=(0.06,P["FP+DESC+MECH"].iloc[0]-.0022),fontsize=8.4,color="#118A6C",fontweight="bold")
ax.annotate(f"ско {P['FP+DESC'].std(ddof=1):.4f}",xy=(3.14,P["FP+DESC"].iloc[3]),fontsize=8.2,color="#6b6b6b",va="center")
ax.annotate(f"ско {P['FP+DESC+MECH'].std(ddof=1):.4f}",xy=(3.14,P["FP+DESC+MECH"].iloc[3]),fontsize=8.2,color="#118A6C",va="center")
ax.set_xticks(xs); ax.set_xticklabels([str(i) for i in range(4)])
ax.set_xlabel("зерно разбиения")
ax.set_ylabel("макро ST-RAE"); ax.grid(axis="y")
ax.set_xlim(-.25,4.75); ax.set_ylim(.7630,.7860)
ax=axes[1]
w=.36
for j,c in enumerate(CYPS):
    a=S[S.features=="FP+DESC"][c].to_numpy(); b=S[S.features=="FP+DESC+MECH"][c].to_numpy()
    ax.bar(j-w/2-.012,a.std(ddof=1),w,color="#c9c9c9")
    ax.bar(j+w/2+.012,b.std(ddof=1),w,color=COL[c])
    ax.annotate(f"{a.std(ddof=1)/b.std(ddof=1):.1f}×",xy=(j,max(a.std(ddof=1),b.std(ddof=1))+.0018),
                ha="center",fontsize=8.6,color=INK,fontweight="bold")
ax.set_xticks(range(4)); ax.set_xticklabels(CYPS)
ax.set_ylabel("ско ST-RAE по зёрнам"); ax.grid(axis="y"); ax.set_ylim(0,.0205)
ax.legend(handles=[Patch(facecolor="#c9c9c9",label="без блока"),
                   Patch(facecolor="#8a8a8a",label="с блоком, цвет фермента")],
          loc="upper right",fontsize=8.1,bbox_to_anchor=(1.02,1.03))
save(fig,"seeds")

# ------------------------------------------------------------------ 7. трудность проверки
print("7 трудность")
nn=np.load("/tmp/verify/nn_seed0.npy"); te=pd.read_csv(D+"test_with_nn.csv") if False else None
from rdkit import Chem, RDLogger, DataStructs
from rdkit.Chem import rdFingerprintGenerator
RDLogger.DisableLog("rdApp.*")
gen=rdFingerprintGenerator.GetMorganGenerator(radius=2,fpSize=2048)
btr=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]
tef=pd.read_csv(D+"cyp-challenge-TEST-BLINDED.csv")
bte=[gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in tef.SMILES]
nte=np.array([max(DataStructs.BulkTanimotoSimilarity(f,btr)) for f in bte])
fig,ax=plt.subplots(figsize=(5.15,2.75))
bins=np.linspace(.15,.95,45)
ax.hist(nn,bins=bins,density=True,color="#cfcfcf",edgecolor="#a8a8a8",lw=.4,label="отложенный кластерный фолд")
ax.hist(nte,bins=bins,density=True,histtype="step",color="#C0432A",lw=1.9,label="настоящий тест")
mn,mt=np.median(nn),np.median(nte)
ax.axvline(mn,color="#6b6b6b",lw=1.1,ls=(0,(3,2)))
ax.axvline(mt,color="#C0432A",lw=1.1,ls=(0,(3,2)))
ax.annotate("",xy=(mt,6.35),xytext=(mn,6.35),
            arrowprops=dict(arrowstyle="<|-|>,head_width=.13,head_length=.32",color=INK,lw=1))
ax.annotate(f"{mn:.3f}",xy=(mn-.010,.28),ha="right",va="bottom",fontsize=8.2,color="#6b6b6b",fontweight="bold")
ax.annotate(f"{mt:.3f}",xy=(mt+.010,.28),va="bottom",fontsize=8.2,color="#C0432A",fontweight="bold")
ax.annotate("0.153",xy=((mn+mt)/2,6.62),ha="center",fontsize=8.2,color=INK,fontweight="bold")
ax.set_xlabel("сходство Танимото с ближайшим обучающим соединением")
ax.set_ylabel("плотность"); ax.legend(loc="upper left",fontsize=8.2,bbox_to_anchor=(-.015,1.02)); ax.grid(axis="y"); ax.set_ylim(0,7.6)
save(fig,"nnsim")

# ------------------------------------------------------------------ 8. корреляции и факторы
print("8 корреляции")
M=inh[[f"{c}_pIC50_direct_inhibition" for c in CYPS]]; M.columns=CYPS
R=np.eye(4); N=np.zeros((4,4),int)
for i in range(4):
    for j in range(i+1,4):
        k=M[CYPS[i]].notna()&M[CYPS[j]].notna(); N[i,j]=N[j,i]=k.sum()
        R[i,j]=R[j,i]=spearmanr(M.loc[k,CYPS[i]],M.loc[k,CYPS[j]]).statistic
def fa(R,r=2,it=300):
    h=np.ones(4)*.5
    for _ in range(it):
        Rr=R.copy(); np.fill_diagonal(Rr,h)
        v,V=np.linalg.eigh(Rr); idx=np.argsort(v)[::-1]; v,V=v[idx],V[:,idx]
        L=V[:,:r]*np.sqrt(np.clip(v[:r],1e-9,None)); h=np.clip((L**2).sum(1),1e-3,.999)
    return L,1-h
L,kap=fa(R)
from matplotlib.colors import LinearSegmentedColormap
cmap=LinearSegmentedColormap.from_list("dv",["#C0432A","#efefec","#3260C4"])
fig,axes=plt.subplots(1,2,figsize=(6.10,2.85),gridspec_kw={"width_ratios":[1,.95],"wspace":.50})
ax=axes[0]
im=ax.imshow(R,cmap=cmap,vmin=-1,vmax=1)
for i in range(4):
    for j in range(4):
        if i==j: continue
        ax.text(j,i,f"{R[i,j]:+.2f}",ha="center",va="center",fontsize=8.6,
                color=SURF if abs(R[i,j])>.55 else INK,fontweight="bold")
    ax.text(i,i,"1",ha="center",va="center",fontsize=8.6,color="#e8eef8")
ax.set_xticks(range(4)); ax.set_xticklabels([c[3:] for c in CYPS])
ax.set_yticks(range(4)); ax.set_yticklabels([c[3:] for c in CYPS])
ax.tick_params(length=0)
cb=fig.colorbar(im,ax=ax,fraction=.042,pad=.035); cb.outline.set_linewidth(0); cb.ax.tick_params(length=2,labelsize=7.4)
ax=axes[1]
ax.axhline(0,color=MUT,lw=.8); ax.axvline(0,color=MUT,lw=.8)
for i,c in enumerate(CYPS):
    ax.annotate("",xy=(L[i,0],L[i,1]),xytext=(0,0),
                arrowprops=dict(arrowstyle="-|>,head_width=.15,head_length=.3",color=COL[c],lw=1.8))
    off={"CYP1A2":(-.10,.075),"CYP2C9":(-.09,.055),"CYP2D6":(.055,-.075),"CYP3A4":(-.055,-.085)}[c]
    ax.annotate(c[3:],xy=(L[i,0]+off[0],L[i,1]+off[1]),ha="center",va="center",fontsize=8.6,
                color=COL[c],fontweight="bold")
ax.set_xlabel("ось 1 — липофильность"); ax.set_ylabel("ось 2 — заряд")
ax.set_xlim(-1.18,.34); ax.set_ylim(-1.00,.40); ax.grid(True,lw=.5)
save(fig,"corr")

# ------------------------------------------------------------------ 9. порог MCC
print("9 порог")
from sklearn.metrics import matthews_corrcoef
from sklearn.linear_model import LogisticRegression
from matplotlib.lines import Line2D
Pp=json.load(open(D+"tdi_probs.json"))
def mccc(tp,tn,fp,fn):
    d=np.sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)); return 0. if d<=0 else (tp*tn-fp*fn)/d
ts=np.linspace(.03,.9,180)
fig,axes=plt.subplots(1,2,figsize=(6.10,2.75),gridspec_kw={"wspace":.26})
rng=np.random.default_rng(3)
for ax,key,col in zip(axes,["CYP3A4","CYP2D6"],["#B98A10","#118A6C"]):
    p_=np.asarray(Pp[key]["p"],float); y=np.asarray(Pp[key]["y"],int); n=len(y)
    g=rng.integers(0,5,n); pc=np.zeros(n)
    z=lambda q: np.log(np.clip(q,1e-6,1-1e-6)/(1-np.clip(q,1e-6,1-1e-6)))
    for f in range(5):
        a,b=g!=f,g==f
        lr=LogisticRegression(C=1e6).fit(z(p_[a]).reshape(-1,1),y[a])
        pc[b]=lr.predict_proba(z(p_[b]).reshape(-1,1))[:,1]
    emp=np.array([matthews_corrcoef(y,(p_>t).astype(int)) for t in ts])
    empc=np.array([matthews_corrcoef(y,(pc>t).astype(int)) for t in ts])
    ax.plot(ts,emp,color="#8a8a8a",lw=1.6)
    ax.plot(ts,empc,color=col,lw=1.9)
    pl=lambda q: ts[int(np.argmax([mccc(q[q>t].sum(),(1-q[q<=t]).sum(),(1-q[q>t]).sum(),q[q<=t].sum()) for t in ts]))]
    t0,t1=pl(p_),pl(pc)
    ax.plot(t0,emp[np.argmin(abs(ts-t0))],"o",color=SURF,mec="#8a8a8a",mew=1.8,ms=8,zorder=5)
    ax.plot(t1,empc[np.argmin(abs(ts-t1))],"o",color=col,ms=8,zorder=5)
    ax.plot(ts[emp.argmax()],emp.max(),"X",color=INK,ms=9,zorder=5)
    lo_,hi_=ax.get_ylim(); ax.set_ylim(lo_,hi_+.10*(hi_-lo_))
    ax.set_title(key,color=col)
    ax.set_xlabel("порог"); ax.grid(axis="y")
axes[0].set_ylabel("MCC")
fig.legend(handles=[
    Line2D([],[],color="#8a8a8a",lw=1.6,label="MCC от порога, сырые вероятности"),
    Line2D([],[],color="#6b6b6b",lw=1.9,label="то же после вложенной калибровки"),
    Line2D([],[],ls="",marker="o",mfc=SURF,mec="#8a8a8a",mew=1.8,ms=8,label="порог, выбранный plug-in по сырым"),
    Line2D([],[],ls="",marker="o",color="#6b6b6b",ms=8,label="порог plug-in после калибровки"),
    Line2D([],[],ls="",marker="X",color=INK,ms=9,label="достижимый максимум"),
], loc="lower center", ncols=3, fontsize=7.7, bbox_to_anchor=(.5,-.32))
save(fig,"mcc")

# ------------------------------------------------------------------ 10. базовая линия
print("10 базовая линия")
from matplotlib.patches import Patch
fig,ax=plt.subplots(figsize=(4.95,2.55))
base=[0.8797,0.7015,0.9927,0.5178]; mech=[0.8786,0.6908,0.9803,0.5194]
xs=np.arange(4); w=.34
ax.bar(xs-w/2-.012,base,w,color="#c9c9c9")
ax.bar(xs+w/2+.012,mech,w,color=[COL[c] for c in CYPS])
ax.axhline(1.0,color="#C0432A",lw=1.2)
ax.annotate("предсказание одного числа для всех соединений",xy=(-.42,1.135),fontsize=8,color="#C0432A")
for i,(a,b) in enumerate(zip(base,mech)):
    ax.annotate(f"{a:.3f}",xy=(i-w/2-.012,a+.014),ha="center",fontsize=7.7,color=MUT)
    ax.annotate(f"{b:.3f}",xy=(i+w/2+.012,b+.014),ha="center",fontsize=7.7,color=INK,fontweight="bold")
ax.set_xticks(xs); ax.set_xticklabels(CYPS); ax.set_ylabel("ST-RAE"); ax.set_ylim(0,1.30)
fig.legend(handles=[Patch(facecolor="#c9c9c9",label="фингерпринт и дескрипторы"),
                    Patch(facecolor=MUT,label="плюс механистический блок, цвет фермента")],
           loc="lower center",ncols=2,fontsize=8.1,bbox_to_anchor=(.5,-.11))
ax.grid(axis="y")
save(fig,"baseline")
print("готово")
