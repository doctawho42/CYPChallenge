"""Co-crystal overlay contrast over the ensemble. Item 295 (pre-registered, blind).

Item 199 built `src/overlay.py`, passed its design check on three enzymes of four, and said in so
many words that the rank ablation was worth running. It never was -- `data/overlay.npz` has sat on
disk since 1 September with nothing reading it. This is that ablation, and it is the cheap first
rung of the structural ladder: it decides whether a POSE carries cavity information that a guessed
cavity (`ablsite`, item 179, CYP2D6 -0.0094) did not, before anyone spends ~20000 docking runs.

Why the quantity is not in the class that keeps returning zero: an overlay onto a cavity's own
BOUND co-crystal ligand is a function of the (ligand, cavity) PAIR. Item 189's six nulls are all
functions of the ligand alone, and item 168 says the interaction form -- ligand columns conditioned
on a known site -- is the one that works.

THE FEATURE IS THE CONTRAST, not the score. The common part of the four overlap scores is how large
and how greasy the molecule is, and the ligand block already carries that in 2295 columns. Centring
each molecule across the four cavities cancels it and leaves complementarity to a particular pocket.
Measured: centring lifts CYP1A2 from +0.149 to +0.173 and collapses CYP2C9/CYP3A4 from +0.077/+0.082
to +0.002 -- on those two the raw score was mostly bulk, which is item 199's size check from the
other side.

Structure of the experiment, mirroring `verify/k85_shape2.py` exactly so the comparison is paired:
only the PER-ENZYME member is rebuilt on FP+DESC+MECH+overlay and dead-zone-passed; the other four
members and the folds come from the cached `members_seed{s}.json`. Everything but the one member
cancels, so the sd is small (item 280).

LIVE CELLS ARE TWO. After items 282-285 the per-enzyme member is kept only on CYP1A2
(поферментно+GP+ствол) and CYP2D6 (all five). CYP2C9 and CYP3A4 ship GP+ствол, so this arm cannot
move them by construction -- they are a harness control that must read exactly 0.0000.

Arms (all append to the shared block; acceptance in item 295):
    A  ЦЕЛЕВАЯ             контраст СВОЕЙ полости, одна колонка
    B  НЕВЕРНАЯ ИЗОФОРМА   контраст ЧУЖОЙ полости (сдвиг на 1) -- the control item 199 demanded,
                           and NOT a permutation: only a wrong-CAVITY column separates real pocket
                           complementarity from a volume descriptor any cavity would supply
    C  НЕНАПРАВЛЕННАЯ      все четыре колонки контраста каждому ферменту -- because item 179 found
                           the undirected arm beat the targeted one on the shape block

Writes results/logs/k93_overlay.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse, json, time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from cypsplit import butina_folds
import submit as S
from submit import CYPS, _oof_one, _dz_oof, _keep

MEM = ["поферментно", "пул", "GP", "гребневая", "ствол"]
FLOOR = [0.0061, 0.0071, 0.0049, 0.0033]


def contrast(S_raw):
    """Per-molecule centring across the four cavities: removes size/greasiness, keeps pocket fit."""
    return S_raw - S_raw.mean(axis=1, keepdims=True)


def per_enzyme_member(Xs, y, LO, HI, mask, fold):
    """Dead-zone-passed per-enzyme member on an augmented matrix. Mirrors k85_shape2.shape_member."""
    P = _oof_one(Xs, y, mask, fold, False)
    out = []
    for e in range(len(CYPS)):
        m = mask[:, e]
        t = np.clip(P[e], LO[m, e], HI[m, e])
        out.append(_dz_oof("поферментно", Xs[m], t, fold[m], e))
    return out


def compose(members, e):
    keep = [np.asarray(members[k][e], float) for k in MEM if _keep(e, k)]
    return np.mean(keep, axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="4,5,6,7",
                    help="сиды с готовым кэшем членов; свежие для ЭТОЙ гипотезы (пункт 295)")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    z = np.load(D + "feats.npz")
    base = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    C = contrast(np.load(D + "overlay.npz", allow_pickle=True)["train"].astype(np.float64))
    print(f"контраст перекрытия: {C.shape}, NaN {int(np.isnan(C).sum())}", flush=True)

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)
    S.LO, S.HI = LO, HI

    # Arm A gives enzyme e its OWN cavity column; arm B gives it the NEXT cavity's (wrong isoform);
    # arm C gives every enzyme all four. A and B differ in exactly one column, which is the point.
    ARMS = {"A целевая": [C[:, [e]] for e in range(4)],
            "B неверная изоформа": [C[:, [(e + 1) % 4]] for e in range(4)],
            "C ненаправленная": [C for _ in range(4)]}

    dst = RES + "logs/k93_overlay.json"
    out = json.load(open(dst)) if _pl.Path(dst).exists() else {"seeds": {}}
    print(f"\n  {'сид':>3s} {'арм':>20s} {'фермент':>8s} {'подаётся':>9s} {'+перекр':>9s} "
          f"{'Δранг':>8s}", flush=True)
    for s in seeds:
        cache = RES + f"preds/members_seed{s}.json"
        if not _pl.Path(cache).exists():
            print(f"  сид {s}: нет {cache}, пропуск", flush=True)
            continue
        mem = json.load(open(cache))
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        res = out["seeds"].get(str(s), {})

        ref = {c: float(spearmanr(y[mask[:, e], e], compose(mem, e)).statistic)
               for e, c in enumerate(CYPS)}
        res["эталон"] = ref

        for arm, blocks in ARMS.items():
            t0 = time.time()
            # Cost control, and it is not cosmetic. k85 could call _oof_one ONCE because its matrix
            # was shared by all four enzymes. Arms A and B give each enzyme a DIFFERENT column, and
            # _oof_one loops over all four internally, so the naive form computes 4x what it needs
            # and throws three quarters away. Two economies, both exact:
            #   1. mask out the other three enzymes -- _oof_one's own `b.sum() == 0` guard then
            #      skips them at no cost, and the enzyme we want is computed identically;
            #   2. skip any enzyme where the per-enzyme member is NOT kept by the shipped SOLO.
            #      There compose() ignores the rebuilt member entirely, so the delta is 0.0000 by
            #      CONSTRUCTION rather than by measurement -- which is what item 295's condition 4
            #      actually asserts. CYP2C9 and CYP3A4 ship GP+ствол, so this halves arms A and B.
            aug = [np.asarray(mem["поферментно"][e], float) for e in range(len(CYPS))]
            shared = all(b is blocks[0] for b in blocks)
            if shared:
                Xs = np.hstack([base, blocks[0]])
                P = _oof_one(Xs, y, mask, fold, False)
                for e in range(len(CYPS)):
                    if not _keep(e, "поферментно"):
                        continue
                    m = mask[:, e]
                    t = np.clip(P[e], LO[m, e], HI[m, e])
                    aug[e] = _dz_oof("поферментно", Xs[m], t, fold[m], e)
            else:
                for e in range(len(CYPS)):
                    if not _keep(e, "поферментно"):
                        continue
                    Xs = np.hstack([base, blocks[e]])
                    one = np.zeros_like(mask); one[:, e] = mask[:, e]
                    P = _oof_one(Xs, y, one, fold, False)
                    m = mask[:, e]
                    t = np.clip(P[e], LO[m, e], HI[m, e])
                    aug[e] = _dz_oof("поферментно", Xs[m], t, fold[m], e)
            memA = dict(mem); memA["поферментно"] = [np.asarray(v).tolist() for v in aug]
            res[arm] = {}
            for e, c in enumerate(CYPS):
                r1 = float(spearmanr(y[mask[:, e], e], compose(memA, e)).statistic)
                res[arm][c] = r1
                print(f"  {s:3d} {arm:>20s} {c:>8s} {ref[c]:9.4f} {r1:9.4f} {r1-ref[c]:+8.4f}",
                      flush=True)
            print(f"      ({time.time()-t0:.0f} c)", flush=True)

        out["seeds"][str(s)] = res
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(dst, "w"), ensure_ascii=False, indent=1)

    # ---- verdict against item 295 -----------------------------------------------------------
    ks = [k for k in out["seeds"] if "эталон" in out["seeds"][k]]
    if not ks:
        return
    print(f"\n  {'арм':>20s} {'фермент':>8s} {'Δранг ср.':>10s} {'sd':>8s} {'знак':>6s} "
          f"{'пол':>8s} {'вердикт':>10s}", flush=True)
    summ = {}
    for arm in ARMS:
        summ[arm] = {}
        for e, c in enumerate(CYPS):
            d = np.array([out["seeds"][k][arm][c] - out["seeds"][k]["эталон"][c] for k in ks
                          if arm in out["seeds"][k]])
            if not len(d):
                continue
            ok = d.mean() > FLOOR[e] and (d > 0).sum() >= 0.75 * len(d)
            summ[arm][c] = {"mean": float(d.mean()),
                            "sd": float(d.std(ddof=1)) if len(d) > 1 else 0.0,
                            "sign": int((d > 0).sum()), "n": int(len(d)), "passes": bool(ok)}
            print(f"  {arm:>20s} {c:>8s} {d.mean():+10.4f} "
                  f"{d.std(ddof=1) if len(d) > 1 else 0:8.4f} {int((d>0).sum())}/{len(d):<2d} "
                  f"{FLOOR[e]:8.4f} {'ПРОХОДИТ' if ok else 'нет':>10s}", flush=True)

    print("\n  === условия пункта 295 ===", flush=True)
    for c in ("CYP1A2", "CYP2D6"):
        A = summ.get("A целевая", {}).get(c)
        B = summ.get("B неверная изоформа", {}).get(c)
        if not A or not B:
            continue
        beats = A["mean"] > B["mean"]
        print(f"  {c}: A {A['mean']:+.4f} (знак {A['sign']}/{A['n']}) против "
              f"B {B['mean']:+.4f} -> целевая {'бьёт' if beats else 'НЕ бьёт'} неверную изоформу; "
              f"порог {'взят' if A['passes'] else 'не взят'}   "
              f"=> {'ПРИНЯТО' if (A['passes'] and beats) else 'ОТКЛОНЕНО'}", flush=True)
    for c in ("CYP2C9", "CYP3A4"):
        A = summ.get("A целевая", {}).get(c)
        if A:
            print(f"  контроль стенда {c}: Δ {A['mean']:+.6f} "
                  f"(обязан быть 0.000000 --- поферментный член туда не входит)", flush=True)

    out["итог"] = summ
    json.dump(out, open(dst, "w"), ensure_ascii=False, indent=1)
    print(f"\n  сохранено: {dst}", flush=True)


if __name__ == "__main__":
    main()
