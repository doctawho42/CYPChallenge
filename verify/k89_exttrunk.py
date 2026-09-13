"""Stage 0: re-score the external-CYP auxiliary-head trunk by RANK over the ensemble. Item 66 gap.

The one artefact that most directly instantiates the "one shared trunk, two tasks, the second
task EXTERNAL" proposal is src/trunkext.py: the twohead trunk with 8004 external ChEMBL CYP pIC50
rows in the auxiliary head instead of the screen. Item 66 measured it only in RAW ST-RAE
(channel -0.0176, 4/4) and never converted it to rank after the affine pair, and never inserted
it into the ensemble as a member. CLAUDE.md / item 77 are explicit that a raw ST-RAE number is
not evidence -- the affine pair rewrites it, reversing three of the five item-77 arms. So the
decisive measurement for this proposal was never taken, exactly the item-267 pattern (shelved by
assumption, the deciding number never computed).

This reads the COMMITTED predictions (results/preds/trunk_ext.json, seeds 0-3, lambda 0 and 3.0)
and puts them through the SAME rank/affine-pair path the knockout ledger uses (verify/k83_loo.py:
score() = mean Spearman as rank, mean ST-RAE after fit_apply as pair), against the SHIPPED
composition (submit.oof_members ансамбль5 + dz_pass + _keep/SOLO). No training; costs seconds.

PRE-REGISTERED before running (stated in chat first):
  criterion for adoption : macro ensemble RANK gain > 0.007, sign 4/4.
  prediction             : does NOT pass -- it lands like the sibling screen trunk (item 79):
                           survives the affine pair standalone, near-parity as a member, and
                           sub-floor over the ensemble, because contrast is already vested in the
                           pooled member (пул knockout +0.0008, item 274) and a second neural trunk
                           lands in the most-tied niche (поферментно-пул rho 0.989, item 280).
  caveat                 : seeds 0-3 are the SAME seeds item 66 used; a PASS would only warrant a
                           fresh-seed confirmation (Stage 1), which needs the external CSVs that are
                           no longer in the tree. A FAIL closes the external-trunk proposal honestly.

Arms measured per seed (all against the shipped 5-member reference):
  standalone       : the external trunk alone (lam 3.0 and lam 0.0) and, for comparison, the
                     shipped screen trunk (ствол) alone.
  INSERT lam3      : add the external trunk as an EXTRA member on every enzyme (the "would a new
                     member survive" question) -> Δrank vs reference.
  INSERT lam0      : control -- external rows present, auxiliary term OFF (bookkeeping only). The
                     channel is INSERT(lam3) - INSERT(lam0).
  SWAP lam3        : replace the screen trunk with the external trunk where ствол is kept.

Writes results/logs/k89_exttrunk.json after each seed.
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
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply
import submit as S
from submit import CYPS, oof_members, dz_pass, _keep, _trunk_clip


def score(P, y, LO, HI, mask, fold):
    """Identical to verify/k83_loo.py: rank = mean Spearman (affine-pair-invariant),
    pair = mean ST-RAE after the per-fold affine pair."""
    rk, pr = [], []
    for e in range(len(CYPS)):
        m = mask[:, e]
        yy, lo, hi = y[m, e], LO[m, e], HI[m, e]
        u = np.ones(len(yy)) / len(yy)
        rk.append(float(spearmanr(yy, P[e]).statistic))
        pr.append(float(strae(yy, fit_apply(P[e], lo, hi, fold[m], u),
                              y_true_upper=hi, y_true_lower=lo)))
    return float(np.mean(rk)), float(np.mean(pr)), rk, pr


def combine_keep(parts, mask):
    """Shipped mean: kept members per _keep/SOLO."""
    return [np.mean([P[e] for k, P in parts if _keep(e, k)], axis=0) for e in range(len(CYPS))]


def combine_insert(parts, extra, mask):
    """Shipped mean PLUS the external trunk as an extra member on EVERY enzyme."""
    return [np.mean([P[e] for k, P in parts if _keep(e, k)] + [extra[e]], axis=0)
            for e in range(len(CYPS))]


def combine_swap(parts, extra, mask):
    """Shipped mean with the external trunk substituted for ствол where ствол is kept;
    on enzymes where ствол is not kept it is inserted (so the trunk-type slot is the external one)."""
    out = []
    for e in range(len(CYPS)):
        keep = [P[e] for k, P in parts if _keep(e, k) and k != "ствол"]
        out.append(np.mean(keep + [extra[e]], axis=0))
    return out


def ext_member(preds, key, y, mask):
    """External-trunk OOF predictions for one (seed|lambda) key, per enzyme, clipped like _oof_trunk."""
    A = np.asarray(preds[key], float)
    return [_trunk_clip(A[mask[:, e], e], y[mask[:, e], e]) for e in range(len(CYPS))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    a = ap.parse_args()

    z = np.load(D + "feats.npz")
    XFULL = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)
    S.LO, S.HI = LO, HI

    ext = json.load(open(RES + "preds/trunk_ext.json"))["preds"]

    dst = RES + "logs/k89_exttrunk.json"
    out = json.load(open(dst)) if _pl.Path(dst).exists() else {}
    print(f"  {'сид':>3s} {'арм':>22s} {'ранг':>8s} {'пара':>8s} {'Δранг':>9s} {'Δпара':>9s}",
          flush=True)
    for s in [int(x) for x in a.seeds.split(",") if x.strip()]:
        if f"{s}|3.0" not in ext:
            print(f"  сид {s}: нет ключа {s}|3.0 в trunk_ext.json, пропуск", flush=True)
            continue
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        res = out.get(str(s), {})
        t0 = time.time()
        parts = oof_members(XFULL, y, mask, fold, "ансамбль5", seed=s, dead=True)
        pdz = dz_pass(parts, XFULL, mask, fold, (LO, HI))
        print(f"      сборка+мёртвая зона: {time.time() - t0:.0f} c", flush=True)

        ref = combine_keep(pdz, mask)
        res["эталон"] = score(ref, y, LO, HI, mask, fold)
        base, basep = res["эталон"][0], res["эталон"][1]
        print(f"  {s:3d} {'эталон (подаётся)':>22s} {base:8.4f} {basep:8.4f} "
              f"{'':>9s} {'':>9s}", flush=True)

        ex3 = ext_member(ext, f"{s}|3.0", y, mask)
        ex0 = ext_member(ext, f"{s}|0.0", y, mask)
        stv = [P for k, P in pdz if k == "ствол"][0]  # shipped screen trunk, per enzyme

        arms = [
            ("standalone ext lam3", [ex3[e] for e in range(4)]),
            ("standalone ext lam0", [ex0[e] for e in range(4)]),
            ("standalone ствол(scr)", [stv[e] for e in range(4)]),
            ("INSERT lam3", combine_insert(pdz, ex3, mask)),
            ("INSERT lam0", combine_insert(pdz, ex0, mask)),
            ("SWAP  lam3", combine_swap(pdz, ex3, mask)),
        ]
        for nm, P in arms:
            sc = score(P, y, LO, HI, mask, fold)
            res[nm] = sc
            standalone = nm.startswith("standalone")
            dr = "" if standalone else f"{sc[0]-base:+9.4f}"
            dp = "" if standalone else f"{sc[1]-basep:+9.4f}"
            print(f"  {s:3d} {nm:>22s} {sc[0]:8.4f} {sc[1]:8.4f} {dr:>9s} {dp:>9s}", flush=True)

        out[str(s)] = res
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(dst, "w"), ensure_ascii=False, indent=1)

    # ---- aggregate over seeds ---------------------------------------------------------------
    ins_arms = ["INSERT lam3", "INSERT lam0", "SWAP  lam3"]
    print(f"\n  {'арм (над ансамблем)':>22s} {'Δранг ср.':>10s} {'sd':>8s} {'знак>0':>8s} "
          f"{'Δпара ср.':>10s}   поле 0.007", flush=True)
    ks = [k for k in out if "эталон" in out[k]]
    for nm in ins_arms:
        kk = [k for k in ks if nm in out[k]]
        if not kk:
            continue
        d = np.array([out[k][nm][0] - out[k]["эталон"][0] for k in kk])
        q = np.array([out[k][nm][1] - out[k]["эталон"][1] for k in kk])
        verdict = "ПРОХОДИТ" if (d.mean() > 0.007 and (d > 0).sum() == len(d)) else "не проходит"
        print(f"  {nm:>22s} {d.mean():+10.4f} {d.std(ddof=1) if len(d) > 1 else 0:8.4f} "
              f"{int((d > 0).sum()):5d}/{len(d):<2d} {q.mean():+10.4f}   {verdict}", flush=True)

    # channel = INSERT(lam3) - INSERT(lam0), the auxiliary head's own contribution over the ensemble
    kk = [k for k in ks if "INSERT lam3" in out[k] and "INSERT lam0" in out[k]]
    if kk:
        ch = np.array([out[k]["INSERT lam3"][0] - out[k]["INSERT lam0"][0] for k in kk])
        print(f"\n  канал вспом. головы над ансамблем (INSERT lam3 - lam0): "
              f"{ch.mean():+.4f} ранга, знак {int((ch > 0).sum())}/{len(ch)}", flush=True)

    # standalone: does the external channel survive the pair standalone? (first half of the prediction)
    kk = [k for k in ks if "standalone ext lam3" in out[k] and "standalone ext lam0" in out[k]]
    if kk:
        sr = np.array([out[k]["standalone ext lam3"][0] - out[k]["standalone ext lam0"][0] for k in kk])
        sv = np.array([out[k]["standalone ствол(scr)"][0] for k in kk])
        ex = np.array([out[k]["standalone ext lam3"][0] for k in kk])
        print(f"  standalone: внешний канал {sr.mean():+.4f} ранга (знак {int((sr > 0).sum())}/{len(sr)}); "
              f"внешний ствол {ex.mean():.4f} против screen-ствола {sv.mean():.4f}", flush=True)


if __name__ == "__main__":
    main()
