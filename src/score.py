"""Score the ablation, with a paired bootstrap over compounds for every delta vs FP+DESC."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd, json, sys
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from scipy.stats import spearmanr
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
oof = json.load(open(RES + "preds/oof.json"))
rows = pd.read_csv(D + "rows.csv")
tr = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv").set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index()
SETS = ["FP", "DESC", "MECH", "DESC+MECH", "FP+DESC", "FP+DESC+MECH"]
truth = {}
for c in CYPS:
    col = f"{c}_pIC50_direct_inhibition"; m = tr[col].notna().to_numpy()
    truth[c] = (tr.loc[m, col].to_numpy(), tr.loc[m, col + "_conf_low"].to_numpy(), tr.loc[m, col + "_conf_high"].to_numpy())

print("=== ST-RAE (ниже — лучше; 1.000 = предсказание среднего) ===")
tab = []
for s in SETS:
    r = {"features": s}
    for c in CYPS:
        y, lo, hi = truth[c]; p = np.array(oof[f"{s}|{c}"])
        r[c] = round(strae(y, p, y_true_upper=hi, y_true_lower=lo), 4)
    r["MACRO"] = round(np.mean([r[c] for c in CYPS]), 4)
    tab.append(r)
T = pd.DataFrame(tab); print(T.to_string(index=False))

print("\n=== Spearman rho ===")
tab = []
for s in SETS:
    r = {"features": s}
    for c in CYPS:
        y, _, _ = truth[c]; r[c] = round(spearmanr(np.array(oof[f"{s}|{c}"]), y).statistic, 3)
    r["MACRO"] = round(np.mean([r[c] for c in CYPS]), 3)
    tab.append(r)
print(pd.DataFrame(tab).to_string(index=False))

print("\n=== Парный бутстрап по соединениям: Δ ST-RAE относительно FP+DESC ===")
print("    (отрицательное = добавление блока помогает; 2000 ресэмплов)")
rng = np.random.default_rng(1); B = 2000
out = []
for s in SETS:
    if s == "FP+DESC": continue
    row = {"features": s}
    for c in CYPS:
        y, lo, hi = truth[c]
        pa = np.array(oof[f"FP+DESC|{c}"]); pb = np.array(oof[f"{s}|{c}"])
        n = len(y); d = np.empty(B)
        for b in range(B):
            i = rng.integers(0, n, n)
            d[b] = (strae(y[i], pb[i], y_true_upper=hi[i], y_true_lower=lo[i])
                    - strae(y[i], pa[i], y_true_upper=hi[i], y_true_lower=lo[i]))
        row[c] = f"{d.mean():+.3f} [{np.percentile(d,2.5):+.3f},{np.percentile(d,97.5):+.3f}]"
        row[c + "_p"] = round(float(np.mean(d > 0)), 3)
    out.append(row)
O = pd.DataFrame(out)
for c in CYPS:
    print(f"\n-- {c} --")
    print(O[["features", c, c + "_p"]].to_string(index=False, header=["набор", "ΔST-RAE [95% ДИ]", "P(хуже)"]))
