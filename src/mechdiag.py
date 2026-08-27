"""Where exactly does the mechanistic block act? Per-isoform diagnostics."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import numpy as np, pandas as pd, json
CYPS=["CYP1A2","CYP2C9","CYP2D6","CYP3A4"]
oof=json.load(open(RES+"preds/oof.json")); rows=pd.read_csv(D+"rows.csv")
tr=pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv").set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index()
from scipy.stats import spearmanr
print("=== Сколько ранговой способности полной модели воспроизводят 28 механистических чисел? ===")
for c in CYPS:
    col=f"{c}_pIC50_direct_inhibition"; m=tr[col].notna().to_numpy(); y=tr.loc[m,col].to_numpy()
    r_mech=spearmanr(np.array(oof[f"MECH|{c}"]),y).statistic
    r_full=spearmanr(np.array(oof[f"FP+DESC|{c}"]),y).statistic
    print(f"  {c}: MECH rho={r_mech:.3f} / FP+DESC rho={r_full:.3f} = {r_mech/r_full*100:.0f}%   "
          f"(28 признаков против 2265)")
print("\n=== Выигрыш от блока отдельно по группам соединений (CYP2D6) ===")
z=np.load(D+"feats.npz"); MECH=z["MECH"]; names=[l.strip() for l in open(D+"mech_names.csv")]
ib=names.index("is_base_74")
c="CYP2D6"; col=f"{c}_pIC50_direct_inhibition"; m=tr[col].notna().to_numpy()
y=tr.loc[m,col].to_numpy(); base=MECH[m][:,ib]>0.5
pa=np.array(oof[f"FP+DESC|{c}"]); pb=np.array(oof[f"FP+DESC+MECH|{c}"])
for lab,sel in [("основание при pH 7.4",base),("не основание",~base)]:
    print(f"  {lab}: n={sel.sum():4d}  медиана pIC50={np.median(y[sel]):.2f}  "
          f"MAE {np.abs(y[sel]-pa[sel]).mean():.3f} -> {np.abs(y[sel]-pb[sel]).mean():.3f}  "
          f"rho {spearmanr(pa[sel],y[sel]).statistic:.3f} -> {spearmanr(pb[sel],y[sel]).statistic:.3f}")
print("\n=== Важность механистических признаков (permutation, MECH-модель, CYP2D6) ===")
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.inspection import permutation_importance
X=MECH[m]
mdl=HistGradientBoostingRegressor(max_iter=300,learning_rate=0.06,random_state=0).fit(X,y)
pi=permutation_importance(mdl,X,y,n_repeats=8,random_state=0,scoring="r2")
imp=pd.DataFrame({"признак":names,"важность":pi.importances_mean.round(4)}).sort_values("важность",ascending=False)
print(imp.head(12).to_string(index=False))
