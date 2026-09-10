"""Honest nested per-enzyme subset enumeration over the five dead-zone-passed members. Item 278.

Item 274 knocked each member out of the finished model and found the pooled member (+0.0008) and
the ridge (+0.0035) net-negative, sign 0/4. This asks the full question: for each enzyme, which of
the 31 non-empty subsets of the five members is best -- chosen HONESTLY, on outer-training folds and
scored on the held-out fold, so the choice never sees the rows it is graded on (item 245).

Members (dead-zone-passed, order fixed by submit.oof_members): поферментно, пул, GP, гребневая,
ствол. A subset's prediction is the unweighted mean of its members, matching the ensemble.

Per seed the five members are rebuilt once (the expensive part) and their per-compound predictions
saved to results/preds/members_seed{s}.json, so the enumeration, the nested selection and the
item-275 error-correlation matrix are all cheap arithmetic on the saved arrays.

Reports, per seed and averaged over four:
  * подаётся      --- shipped composition scored honestly (SOLO GP on CYP3A4, five-mean else);
  * нест. выбор   --- honest nested per-enzyme subset selection;
  * потолок       --- in-sample per-enzyme best subset (item 245 contamination, ceiling only);
  * все пять      --- reference.
Acceptance conditions are in item 278 and not repeated. Writes results/logs/k84_subsets.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse, json, time
from itertools import combinations
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from cypsplit import butina_folds
import submit as S
from submit import CYPS, oof_members, dz_pass, _keep

MEMBERS = ["поферментно", "пул", "GP", "гребневая", "ствол"]
SUBSETS = [s for k in range(1, 6) for s in combinations(range(5), k)]   # 31 non-empty


def build(seed, X, y, LO, HI, mask):
    """Five dead-zone-passed members, per-compound, as {name: [per-enzyme arrays]}."""
    fold, _ = butina_folds(list(pd.read_csv(D + "rows.csv").SMILES), seed=seed)
    parts = oof_members(X, y, mask, fold, "ансамбль5", seed=seed, dead=True)
    parts = dz_pass(parts, X, mask, fold, (LO, HI))
    d = {name: [np.asarray(P[e], float).tolist() for e in range(4)] for name, P in parts}
    return d, fold


def mean_subset(members, e, idx):
    return np.mean([np.asarray(members[MEMBERS[i]][e], float) for i in idx], axis=0)


def rho(members, e, idx, ytrue, sel):
    p = mean_subset(members, e, idx)[sel]
    return float(spearmanr(ytrue[sel], p).statistic)


def score_config(members, y, mask, fold, chooser):
    """chooser(e, selrows) -> subset idx, applied nested. Returns (macro rank, per-enzyme rank)."""
    rk = []
    for e in range(4):
        m = mask[:, e]
        yy, fi = y[m, e], fold[m]
        out = np.zeros(len(yy))
        for f in range(5):
            tr = fi != f
            te = fi == f
            if te.sum() == 0:
                continue
            idx = chooser(e, tr)
            out[te] = mean_subset(members, e, idx)[te]
        rk.append(float(spearmanr(yy, out).statistic))
    return float(np.mean(rk)), rk


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    a = ap.parse_args()
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)
    S.LO, S.HI = LO, HI
    SHIP = {e: [i for i, nm in enumerate(MEMBERS) if _keep(e, nm)] for e in range(4)}

    dst = RES + "logs/k84_subsets.json"
    out = json.load(open(dst)) if _pl.Path(dst).exists() else {"seeds": {}}
    for s in [int(x) for x in a.seeds.split(",") if x.strip()]:
        t0 = time.time()
        mempath = RES + f"preds/members_seed{s}.json"
        if _pl.Path(mempath).exists():
            members = json.load(open(mempath)); fold, _ = butina_folds(list(rows.SMILES), seed=s)
            print(f"сид {s}: члены из кэша", flush=True)
        else:
            members, fold = build(s, X, y, LO, HI, mask)
            json.dump(members, open(mempath, "w"))
            print(f"сид {s}: члены собраны за {time.time()-t0:.0f} c", flush=True)

        # baseline (shipped, no selection), nested selection, in-sample ceiling
        ship_r, ship_pe = score_config(members, y, mask, fold, lambda e, sel: SHIP[e])
        allp_r, _ = score_config(members, y, mask, fold, lambda e, sel: list(range(5)))

        picks = {e: [] for e in range(4)}
        def nested(e, sel):
            yy = y[mask[:, e], e]
            best = max(SUBSETS, key=lambda idx: rho(members, e, idx, yy, sel))
            picks[e].append(best); return list(best)
        nest_r, nest_pe = score_config(members, y, mask, fold, nested)

        def insample(e, sel):
            yy = y[mask[:, e], e]
            allrows = np.ones(len(yy), bool)
            return list(max(SUBSETS, key=lambda idx: rho(members, e, idx, yy, allrows)))
        ceil_r, _ = score_config(members, y, mask, fold, insample)

        res = {"подаётся": [ship_r, ship_pe], "нест": [nest_r, nest_pe],
               "потолок": ceil_r, "все5": allp_r,
               "picks": {str(e): [list(p) for p in picks[e]] for e in range(4)},
               "ship_subset": {str(e): SHIP[e] for e in range(4)}}
        out["seeds"][str(s)] = res
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(dst, "w"), ensure_ascii=False, indent=1)
        print(f"  подаётся {ship_r:.4f}  нест {nest_r:.4f}  Δ {nest_r-ship_r:+.4f}   "
              f"потолок {ceil_r:.4f}  все5 {allp_r:.4f}", flush=True)
        for e, c in enumerate(CYPS):
            names = ["+".join(MEMBERS[i][:4] for i in p) for p in picks[e]]
            print(f"    {c}: подаётся {'+'.join(MEMBERS[i][:4] for i in SHIP[e])}  "
                  f"нест-выбор по фолдам [{', '.join(names)}]", flush=True)

    ks = [k for k in out["seeds"]]
    if ks:
        d = np.array([out["seeds"][k]["нест"][0] - out["seeds"][k]["подаётся"][0] for k in ks])
        print(f"\nнест. минус подаётся, макро-ранг: {d.mean():+.4f}  sd "
              f"{d.std(ddof=1) if len(d)>1 else 0:.4f}  знак {int((d>0).sum())}/{len(d)}  "
              f"пол 0.0052", flush=True)
        print("условия пункта 278: (1) >0.0052 --- {}; (2) знак>=3/4 --- {}".format(
            "ДА" if d.mean() > 0.0052 else "НЕТ",
            "ДА" if (d > 0).sum() >= 0.75 * len(d) else "НЕТ"), flush=True)
        # стабильность выбора
        print("\nустойчивость поферментного выбора (субсет: сколько из 5x{} фолдо-сидов)".format(len(ks)),
              flush=True)
        from collections import Counter
        for e, c in enumerate(CYPS):
            cnt = Counter(tuple(p) for k in ks for p in out["seeds"][k]["picks"][str(e)])
            top, n = cnt.most_common(1)[0]
            tot = sum(cnt.values())
            print(f"  {c}: чаще всего {'+'.join(MEMBERS[i][:4] for i in top)}  {n}/{tot}", flush=True)

    # item-275 correlation matrix from the last seed's saved members
    print("\nматрица корреляций ОШИБОК членов (пункт 275), усреднённая по ферментам:", flush=True)
    corr = np.zeros((5, 5)); w = 0
    for e in range(4):
        m = mask[:, e]; yy = y[m, e]
        E = np.stack([np.asarray(members[MEMBERS[i]][e], float) - yy for i in range(5)])
        corr += np.corrcoef(E); w += 1
    corr /= w
    print("        " + "  ".join(f"{n[:4]:>6s}" for n in MEMBERS), flush=True)
    for i, n in enumerate(MEMBERS):
        print(f"  {n[:6]:>6s} " + "  ".join(f"{corr[i, j]:6.3f}" for j in range(5)), flush=True)
    out["err_corr"] = corr.tolist()
    json.dump(out, open(dst, "w"), ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
