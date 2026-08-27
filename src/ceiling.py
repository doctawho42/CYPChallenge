"""How much signal is there to learn at all? Reliability from the per-compound curve-fit std."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd, sys
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
tr = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
rows = []
for c in CYPS:
    col = f"{c}_pIC50_direct_inhibition"
    d = tr[tr[col].notna()]
    y = d[col].to_numpy(); s = d[col + "_std"].to_numpy()
    lo = d[col + "_conf_low"].to_numpy(); hi = d[col + "_conf_high"].to_numpy()
    var_obs = y.var(ddof=1); var_noise = np.nanmean(s ** 2)
    rel = max(0.0, 1 - var_noise / var_obs)
    # ST-RAE of a perfect-but-noisy oracle: predict the true value plus one draw of assay noise
    rng = np.random.default_rng(0)
    orc = [strae(y, y + rng.normal(0, np.nan_to_num(s, nan=np.nanmedian(s))),
                 y_true_upper=hi, y_true_lower=lo) for _ in range(20)]
    # ST-RAE of the constant predictor = 1.0 by construction; what does a "shrunk" constant give?
    rows.append(dict(cyp=c, n=len(y), sd_y=round(float(np.sqrt(var_obs)), 3),
                     iqr=round(float(np.percentile(y, 75) - np.percentile(y, 25)), 2),
                     median_std=round(float(np.nanmedian(s)), 3),
                     rms_noise=round(float(np.sqrt(var_noise)), 3),
                     reliability=round(rel, 3), max_R2=round(rel, 3),
                     STRAE_noisy_oracle=round(float(np.mean(orc)), 3),
                     mean_CI_width=round(float(np.nanmean(hi - lo)), 2),
                     denom_share_lt4=round(float(
                         np.sum(np.clip(y.mean() - hi, 0, None)[y < 4] + np.clip(lo - y.mean(), 0, None)[y < 4]) /
                         np.sum(np.clip(y.mean() - hi, 0, None) + np.clip(lo - y.mean(), 0, None))), 3)))
r = pd.DataFrame(rows)
print(r.to_string(index=False))
print("""
reliability = 1 - E[sigma^2]/Var(y): доля дисперсии меток, не объяснимая шумом фита кривой.
Это потолок R^2 для любой модели на этих метках.
STRAE_noisy_oracle: ST-RAE «идеального» предсказателя, знающего истину с точностью одного
повторного измерения. Ниже этого не спуститься даже теоретически.
denom_share_lt4: какая доля ЗНАМЕНАТЕЛЯ ST-RAE приходится на соединения с pIC50 < 4.""")
