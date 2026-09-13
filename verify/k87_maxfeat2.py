"""max_features arm 2 over the ensemble, fresh seeds 4-7. Item 286 (pre-registered, blind).

The per-enzyme member with max_features=0.1 helps CYP2D6 by +0.0149 at member level (item 273/k82),
but member-level gains die over the ensemble (items 269, 273, 281). Since the SOLO changes (items
282-285) the per-enzyme member ships only on CYP1A2 (пофе+GP+ствол) and CYP2D6 (all five), so this
tests those two cells.

Rebuilds ONLY the per-enzyme member with max_features=0.1 and its dead-zone pass, under scikit-learn
1.8.0, and substitutes it into the shipped composition using the four other members cached in
members_seed{4-7}.json (built under 1.3.2 = mf=1.0 under 1.8.0, bit-identical, item 273). Baseline is
the cached plain per-enzyme member. Claim is per-enzyme on CYP2D6 (floor 0.0049) and CYP1A2 (0.0061).

MUST run under the overlay:  uv run --with 'scikit-learn==1.8.0' python verify/k87_maxfeat2.py
Writes results/logs/k87_maxfeat2.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json, time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor

from cypsplit import butina_folds
import submit as S
from submit import CYPS, _keep, DZ_KW

MEM = ["поферментно", "пул", "GP", "гребневая", "ствол"]
MF = 0.1
# gbm_reg() pins plus max_features, and DZ_KW plus max_features for the dead-zone refit -- the
# whole member is "the mf=0.1 booster", base fit and refit alike.
BASE_KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31, l2_regularization=1.0,
               random_state=0, max_features=MF)
DZ_MF = dict(DZ_KW, max_features=MF)


def perenz_mf(X, Y, LO, HI, mask, fold):
    """Per-enzyme member with max_features=0.1, dead-zone-passed, per enzyme."""
    out = []
    for e in range(4):
        m = mask[:, e]
        Xi, yy, fi = X[m], Y[m, e], fold[m]
        P = np.zeros(len(yy))
        for f in range(5):                       # base OOF with mf
            a, b = fi != f, fi == f
            if b.sum() == 0:
                continue
            P[b] = HistGradientBoostingRegressor(**BASE_KW).fit(Xi[a], yy[a]).predict(Xi[b])
        t = np.clip(P, LO[m, e], HI[m, e])        # dead-zone target
        Q = np.zeros(len(yy))
        for f in range(5):                       # dz refit with mf under absolute_error
            a, b = fi != f, fi == f
            if b.sum() == 0:
                continue
            Q[b] = HistGradientBoostingRegressor(**DZ_MF).fit(Xi[a], t[a]).predict(Xi[b])
        out.append(Q)
    return out


def compose(members, e):
    return np.mean([np.asarray(members[k][e], float) for k in MEM if _keep(e, k)], axis=0)


def main():
    import sklearn
    if tuple(int(x) for x in sklearn.__version__.split(".")[:2]) < (1, 4):
        raise SystemExit(f"scikit-learn {sklearn.__version__} без max_features; запускайте под "
                         f"overlay --with 'scikit-learn==1.8.0'")
    print(f"scikit-learn {sklearn.__version__}\n", flush=True)
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(Y)
    S.LO, S.HI = LO, HI
    FLOOR = {"CYP1A2": 0.0061, "CYP2D6": 0.0049}
    LIVE = ["CYP1A2", "CYP2D6"]                    # cells that still use the per-enzyme member

    out = {"seeds": {}}
    print(f"  {'сид':>3s} {'фермент':>8s} {'подаётся':>9s} {'mf=0.1':>9s} {'Δранг':>8s}", flush=True)
    for s in (4, 5, 6, 7):
        mem = json.load(open(RES + f"preds/members_seed{s}.json"))
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        t0 = time.time()
        pm = perenz_mf(X, Y, LO, HI, mask, fold)
        memM = dict(mem); memM["поферментно"] = [pm[e].tolist() for e in range(4)]
        res = {}
        for c in LIVE:
            e = CYPS.index(c); m = mask[:, e]; yy = Y[m, e]
            r0 = float(spearmanr(yy, compose(mem, e)).statistic)
            r1 = float(spearmanr(yy, compose(memM, e)).statistic)
            res[c] = [r0, r1]
            print(f"  {s:3d} {c:>8s} {r0:9.4f} {r1:9.4f} {r1-r0:+8.4f}", flush=True)
        out["seeds"][str(s)] = res
        print(f"      (член пересобран за {time.time()-t0:.0f} c)", flush=True)
        json.dump(out, open(RES + "logs/k87_maxfeat2.json", "w"), ensure_ascii=False, indent=1)

    ks = list(out["seeds"])
    print(f"\n  {'фермент':>8s} {'Δранг ср.':>10s} {'sd':>8s} {'знак':>6s} {'пол':>8s} {'вердикт':>10s}",
          flush=True)
    adopt = False
    for c in LIVE:
        d = np.array([out["seeds"][k][c][1] - out["seeds"][k][c][0] for k in ks])
        p = d.mean() > FLOOR[c] and (d > 0).sum() >= 3
        adopt = adopt or p
        print(f"  {c:>8s} {d.mean():+10.4f} {d.std(ddof=1):8.4f} {int((d>0).sum())}/{len(ks):<2d}"
              f" {FLOOR[c]:8.4f} {'ПРОХОДИТ' if p else 'нет':>10s}", flush=True)
    print(f"\n  условие пункта 286 (2D6>0.0049 или 1A2>0.0061, знак>=3/4): "
          f"{'ДА, ПРИНЯТЬ' if adopt else 'НЕТ'}", flush=True)
    out["adopt"] = bool(adopt)
    json.dump(out, open(RES + "logs/k87_maxfeat2.json", "w"), ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
