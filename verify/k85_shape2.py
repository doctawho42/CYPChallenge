"""Shape block, arm 2 over the ensemble. Item 279 (pre-registered, blind).

The +0.0087 on CYP2D6 (item 165) is on the BARE per-enzyme member. Items 269/273/274 show
member-level gains dying over the ensemble, so the question is whether a shape-augmented per-enzyme
member helps 2D6's rank inside the shipped five-member composition.

Reuses the members cached by k84 (results/preds/members_seed{s}.json, dead-zone-passed). Only the
per-enzyme member is rebuilt on FP+DESC+MECH+shape3d and dead-zone-passed; the other four members
and the folds are identical, so this is the paired shared-structure comparison item 280 argued for
(the sd is small because everything but the one member cancels).

Claim is PER-ENZYME on 2D6 against floor 0.0049, not macro. 3A4 ships GP alone (SOLO), so the
per-enzyme member -- shape or plain -- does not enter it; 2C9 is near zero at member level. Live
cells: 2D6, 1A2. Acceptance in item 279. Writes results/logs/k85_shape2.json.
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

from cypsplit import butina_folds
import submit as S
from submit import CYPS, _oof_one, _dz_oof, _keep

MEM = ["поферментно", "пул", "GP", "гребневая", "ствол"]


def shape_member(Xs, y, LO, HI, mask, fold):
    """Dead-zone-passed per-enzyme member on the shape-augmented matrix."""
    P = _oof_one(Xs, y, mask, fold, False)            # plain OOF per enzyme
    out = []
    for e in range(4):
        m = mask[:, e]
        t = np.clip(P[e], LO[m, e], HI[m, e])
        out.append(_dz_oof("поферментно", Xs[m], t, fold[m], e))
    return out


def compose(members, e):
    """Shipped rule: mean of members kept on enzyme e (SOLO applies)."""
    keep = [np.asarray(members[k][e], float) for k in MEM if _keep(e, k)]
    return np.mean(keep, axis=0)


def main():
    z = np.load(D + "feats.npz")
    Xs = np.hstack([z["FP"], z["DESC"], z["MECH"], np.load(D + "shape3d.npz")["train"]])
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)
    S.LO, S.HI = LO, HI
    FLOOR = [0.0061, 0.0071, 0.0049, 0.0033]

    out = {"seeds": {}}
    print(f"  {'сид':>3s} {'фермент':>8s} {'подаётся':>9s} {'+форма':>9s} {'Δранг':>8s}", flush=True)
    for s in range(4):
        mem = json.load(open(RES + f"preds/members_seed{s}.json"))
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        t0 = time.time()
        shp = shape_member(Xs, y, LO, HI, mask, fold)
        memS = dict(mem); memS["поферментно"] = [shp[e].tolist() for e in range(4)]
        res = {}
        for e, c in enumerate(CYPS):
            m = mask[:, e]; yy = y[m, e]
            r0 = float(spearmanr(yy, compose(mem, e)).statistic)
            r1 = float(spearmanr(yy, compose(memS, e)).statistic)
            res[c] = [r0, r1]
            print(f"  {s:3d} {c:>8s} {r0:9.4f} {r1:9.4f} {r1-r0:+8.4f}", flush=True)
        out["seeds"][str(s)] = res
        print(f"      (форменный член пересобран за {time.time()-t0:.0f} c)", flush=True)
        json.dump(out, open(RES + "logs/k85_shape2.json", "w"), ensure_ascii=False, indent=1)

    print(f"\n  {'фермент':>8s} {'Δранг ср.':>10s} {'sd':>8s} {'знак':>6s} {'пол':>8s} {'вердикт':>10s}",
          flush=True)
    ks = list(out["seeds"])
    verdict = {}
    for e, c in enumerate(CYPS):
        d = np.array([out["seeds"][k][c][1] - out["seeds"][k][c][0] for k in ks])
        passes = d.mean() > FLOOR[e] and (d > 0).sum() >= 0.75 * len(d)
        verdict[c] = bool(passes)
        print(f"  {c:>8s} {d.mean():+10.4f} {d.std(ddof=1):8.4f} {int((d>0).sum())}/{len(d):<2d}"
              f" {FLOOR[e]:8.4f} {'ПРОХОДИТ' if passes else 'нет':>10s}", flush=True)
    d2 = np.array([out["seeds"][k]["CYP2D6"][1] - out["seeds"][k]["CYP2D6"][0] for k in ks])
    print(f"\n  условие пункта 279 (2D6 > 0.0049 при знаке >=3/4): "
          f"{'ДА' if verdict['CYP2D6'] else 'НЕТ'}   (2D6 Δ {d2.mean():+.4f}, член давал +0.0087)",
          flush=True)
    out["verdict"] = verdict
    json.dump(out, open(RES + "logs/k85_shape2.json", "w"), ensure_ascii=False, indent=1)


if __name__ == "__main__":
    main()
