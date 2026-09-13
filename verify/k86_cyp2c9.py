"""CYP2C9 composition: fixed subset from seeds 0-3, tested on FRESH seeds 4-7. Item 282.

Item 280 found CYP2C9's shipped all-five composition ninth of 31 subsets; the best on seeds 0-3 is
GP+ствол (+0.0124 in-sample). That is a hypothesis from seeds 0-3, so the honest test fixes the
subset and scores it on seeds it never saw. Builds the five dead-zone-passed members for seeds
4,5,6,7 (reusing k84_subsets.build), saves them, and reports:
  * primary: fixed GP+ствол vs shipped all-five on CYP2C9, four fresh seeds, against floor 0.0071;
  * secondary: full 31-subset enumeration on the fresh seeds -- is GP+ствол still the winner;
  * by-product: does GP+ствол beat GP-alone on CYP3A4 (item 218 vs item 280 contradiction).
Acceptance in item 282. Writes results/logs/k86_cyp2c9.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json, time
from itertools import combinations
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from cypsplit import butina_folds
from k84_subsets import build, MEMBERS
from submit import CYPS

SUBSETS = [s for k in range(1, 6) for s in combinations(range(5), k)]
GP_ST = (MEMBERS.index("GP"), MEMBERS.index("ствол"))
GP = (MEMBERS.index("GP"),)
ALL5 = tuple(range(5))
FRESH = [4, 5, 6, 7]


def rank_subset(members, e, idx, yy):
    p = np.mean([np.asarray(members[MEMBERS[i]][e], float) for i in idx], axis=0)
    return float(spearmanr(yy, p).statistic)


def main():
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)
    import submit as S
    S.LO, S.HI = LO, HI
    e2c9 = CYPS.index("CYP2C9"); e3a4 = CYPS.index("CYP3A4")
    y2 = y[mask[:, e2c9], e2c9]; y3 = y[mask[:, e3a4], e3a4]

    out = {"seeds": {}}
    print(f"  {'сид':>3s} {'2C9 all5':>9s} {'2C9 GP+ст':>10s} {'Δ2C9':>8s} | "
          f"{'3A4 GP':>8s} {'3A4 GP+ст':>10s} {'Δ3A4':>8s}", flush=True)
    for s in FRESH:
        mempath = RES + f"preds/members_seed{s}.json"
        t0 = time.time()
        if _pl.Path(mempath).exists():
            members = json.load(open(mempath)); src = "кэш"
        else:
            members, _ = build(s, X, y, LO, HI, mask)
            json.dump(members, open(mempath, "w")); src = f"{time.time()-t0:.0f} c"
        r_all = rank_subset(members, e2c9, ALL5, y2)
        r_gpst = rank_subset(members, e2c9, GP_ST, y2)
        r3_gp = rank_subset(members, e3a4, GP, y3)
        r3_gpst = rank_subset(members, e3a4, GP_ST, y3)
        # full fresh enumeration on 2C9
        best = max(SUBSETS, key=lambda idx: rank_subset(members, e2c9, idx, y2))
        out["seeds"][str(s)] = {"2c9_all5": r_all, "2c9_gpst": r_gpst,
                                "3a4_gp": r3_gp, "3a4_gpst": r3_gpst,
                                "2c9_fresh_best": list(best)}
        print(f"  {s:3d} {r_all:9.4f} {r_gpst:10.4f} {r_gpst-r_all:+8.4f} | "
              f"{r3_gp:8.4f} {r3_gpst:10.4f} {r3_gpst-r3_gp:+8.4f}   [{src}] "
              f"свеж.лучший 2C9: {'+'.join(MEMBERS[i][:4] for i in best)}", flush=True)
        json.dump(out, open(RES + "logs/k86_cyp2c9.json", "w"), ensure_ascii=False, indent=1)

    ks = list(out["seeds"])
    d2 = np.array([out["seeds"][k]["2c9_gpst"] - out["seeds"][k]["2c9_all5"] for k in ks])
    d3 = np.array([out["seeds"][k]["3a4_gpst"] - out["seeds"][k]["3a4_gp"] for k in ks])
    print(f"\n  CYP2C9 GP+ствол минус all5:  {d2.mean():+.4f}  sd {d2.std(ddof=1):.4f}  "
          f"знак {int((d2>0).sum())}/{len(ks)}  пол 0.0071", flush=True)
    print(f"  условие пункта 282 (>0.0071 при знаке >=3/4): "
          f"{'ДА, ПРИНЯТЬ' if d2.mean()>0.0071 and (d2>0).sum()>=3 else 'НЕТ'}", flush=True)
    freshbest = [tuple(out["seeds"][k]["2c9_fresh_best"]) for k in ks]
    print(f"  свежая энумерация 2C9 выбирает GP+ствол в {sum(b==GP_ST for b in freshbest)}/{len(ks)} сидах",
          flush=True)
    print(f"\n  побочно, CYP3A4 GP+ствол минус GP-один: {d3.mean():+.4f}  sd {d3.std(ddof=1):.4f}  "
          f"знак {int((d3>0).sum())}/{len(ks)}  (пункт 280 против 218; пол 0.0033)", flush=True)
    out["verdict"] = {"c9_gain": d2.tolist(), "a4_gain": d3.tolist(),
                      "adopt_c9": bool(d2.mean() > 0.0071 and (d2 > 0).sum() >= 3)}
    json.dump(out, open(RES + "logs/k86_cyp2c9.json", "w"), ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
