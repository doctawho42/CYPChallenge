"""Пересчёт КАЖДОГО числа, которое стоит в документах, прямо из исходных таблиц."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd, json, sys
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from scipy.stats import spearmanr
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
res=[]
def say(claim, got, want, tol=None, note=""):
    if tol is None: ok = (str(got)==str(want))
    else: ok = abs(float(got)-float(want))<=tol
    res.append(ok)
    print(f"{'OK  ' if ok else 'РАСХОЖДЕНИЕ'}  {claim:52s} в тексте {want!s:>8s}  посчитано {got!s:>8s} {note}")

inh=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
tst=pd.read_csv(D+"cyp-challenge-TEST-BLINDED.csv")
tdi=pd.read_csv(D+"cyp-challenge-TRAIN_TDI.csv")
ema=pd.read_csv(D+"cyp-challenge-TRAIN_Emax.csv")
sc =pd.read_csv(D+"cyp-challenge-single-concentration-TRAIN.csv")
print("=== размеры таблиц ==="); 
say("строк в TRAIN_inhibition", len(inh), 4905)
say("строк в TEST-BLINDED",     len(tst), 750)
say("строк в TRAIN_TDI",        len(tdi), 6145)
say("строк в TRAIN_Emax",       len(ema), 6146)
say("строк в single-concentration", len(sc), 17504)
say("уникальных соединений в скрине", sc.Molecule_Name.nunique(), 4376)
_col = next((c for c in ("concentration_M", "concentration") if c in sc), None)
say("концентрация скрина, М", f"{sc[_col].iloc[0]:.2e}" if _col else "н/д", "4.95e-05",
    note=f"(колонка {_col})" if _col else "(колонки концентрации нет)")

print("\n=== число меток по ферментам (таблица базовой линии) ===")
for c,want in zip(CYPS,[1412,1285,1493,2335]):
    say(f"меток pIC50 {c}", int(inh[f"{c}_pIC50_direct_inhibition"].notna().sum()), want)
say("сумма меток", int(sum(inh[f'{c}_pIC50_direct_inhibition'].notna().sum() for c in CYPS)), 6525)

print("\n=== межквартильный размах pIC50 ===")
iqr={}
# Порядок CYPS: 1A2, 2C9, 2D6, 3A4. Раньше здесь стояло [0.83, 0.92, 1.00, 1.50], где
# первое и третье переставлены относительно таблицы в docs/tex/s14.tex; документ и
# данные согласны, отставал список. Проверка сообщала о двух расхождениях, которых нет.
for c,want in zip(CYPS,[1.00,0.92,0.83,1.50]):
    y=inh[f"{c}_pIC50_direct_inhibition"].dropna().to_numpy()
    iqr[c]=np.percentile(y,75)-np.percentile(y,25)
    say(f"IQR {c}", round(iqr[c],2), want, tol=0.005)

print("\n=== баланс классов TDI ===")
for col,want_pct,want_n in [("CYP3A4_is_TDI",21.3,3584),("CYP2D6_is_TDI",21.6,1497)]:
    if col in tdi:
        d=tdi[col].dropna()
        say(f"размечено {col}", len(d), want_n)
        say(f"доля положительных {col}, %", round(100*d.mean(),1), want_pct, tol=0.06)
    else: print("  нет колонки", col, "| есть:", [x for x in tdi.columns if "TDI" in x])

print("\n=== асимметрия доверительной полосы (b-a), медиана ===")
mx=0
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; d=inh[inh[col].notna()]
    y=d[col].to_numpy(); a=y-d[col+"_conf_low"].to_numpy(); b=d[col+"_conf_high"].to_numpy()-y
    for lab,m in [("y>=5",y>=5),("4<=y<5",(y>=4)&(y<5)),("y<4",y<4)]:
        if m.sum()<10: continue
        v=float(np.median(b[m]-a[m])); mx=max(mx,abs(v))
        print(f"    {c} {lab:8s} n={m.sum():4d}  медиана b-a = {v:+.3f}")
say("максимум |медиана b-a| по всем ячейкам", round(mx,2), 0.13, tol=0.005)

print("\n=== модель только на механистическом блоке: доля ранжирующей способности ===")
oof=json.load(open(RES+"preds/oof.json"))
rows=pd.read_csv(D+"rows.csv")
tr=inh.set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index()
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; m=tr[col].notna().to_numpy(); y=tr.loc[m,col].to_numpy()
    rM=spearmanr(np.asarray(oof[f"MECH|{c}"]),y).statistic
    rA=spearmanr(np.asarray(oof[f"FP+DESC|{c}"]),y).statistic
    rF=spearmanr(np.asarray(oof[f"FP+DESC+MECH|{c}"]),y).statistic
    print(f"    {c}: MECH {rM:.3f} | FP+DESC {rA:.3f} ({100*rM/rA:.0f}%) | FP+DESC+MECH {rF:.3f} ({100*rM/rF:.0f}%)")
print(json.dumps({"итог":"проверено выше"},ensure_ascii=False)[:0])
print(f"\nСОШЛОСЬ {sum(res)} из {len(res)}")
