"""Stage 0, second half: re-score the MULTI-TASK trees (item 188) by RANK over the ENSEMBLE.

Item 188 measured multi-task trees (one split structure, four values per leaf) at +0.0069 macro
rank over the per-enzyme reference AT MEMBER LEVEL, four seeds -- a genuine positive. But it was
never inserted into the CURRENT shipped ensemble and scored by rank over it. That is the symmetric
gap to k89 (the external trunk): "does a multi-task representation earn a place in the ensemble?"
measured in the decisive currency, not by a member-level number plus a ceiling argument.

The multi-task predictions are COMMITTED at results/preds/oof_multi.json (seeds 0-3, arms
независимо / пул / многозадачно / многозадачно+масштаб). This reads them and puts the honest-form
member through the SAME rank/affine-pair path k89/k83_loo use, against the SHIPPED composition.

HONEST FORM (item 188): multi-task on CYP1A2/CYP2C9/CYP2D6, but INDEPENDENT (per-enzyme) on
CYP3A4, where the data-rich task shows ordinary negative transfer (независимо beats многозадачно
on 3A4 in every seed of oof_multi.json's table).

WIRING NOTE. ablmulti.py's learner is DecisionTree-based boosting (NTREE=200, DEPTH=5, MAXFEAT=0.3),
NOT the shipped HistGradientBoosting "поферментно" member, and its predictions carry no dead-zone
pass. So this measures INSERTION with a control rather than a clean swap:
  INSERT multi   : add the honest-form multi member as an EXTRA member on every enzyme -> Δrank.
  INSERT незав   : same learner, per-enzyme -> the control. The multi-task-specific channel over
                   the ensemble is INSERT(multi) - INSERT(незав), exactly parallel to k89's λ3-λ0.
  SWAP           : replace the shipped "поферментно" member with the multi member where поферментно
                   is kept (1A2, 2D6 only; 2C9/3A4 do not keep it, so SWAP=INSERT there).

PRE-REGISTERED before running (stated in chat first):
  criterion for adoption : macro ensemble RANK gain > 0.007, sign 4/4.
  prediction             : INSERT multi does NOT pass -- the multi member sits in the boosters' niche
                           (ρ≈0.95, item 289) next to the per-enzyme booster and the pool, so the
                           +0.0069 member-level gain is eaten by redundancy over the ensemble (item
                           269: escaping dilution is not sufficient; the contrast is already vested,
                           пул knockout +0.0008, item 274). The control INSERT(незав) should be <= 0
                           and the channel INSERT(multi)-INSERT(незав) small and positive.
  caveat                 : seeds 0-3 are the same seeds item 188 used; a PASS warrants a fresh-seed
                           confirmation. A FAIL closes "multi-task earns an ensemble place" honestly.

Writes results/logs/k90_multitask.json after each seed.
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
from submit import CYPS, oof_members, dz_pass, _keep

# item 188 honest form: multi on these, independent on the rest
MULTI_ENZ = {"CYP1A2", "CYP2C9", "CYP2D6"}


def score(P, y, LO, HI, mask, fold):
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
    return [np.mean([P[e] for k, P in parts if _keep(e, k)], axis=0) for e in range(len(CYPS))]


def combine_insert(parts, extra, mask):
    return [np.mean([P[e] for k, P in parts if _keep(e, k)] + [extra[e]], axis=0)
            for e in range(len(CYPS))]


def combine_swap(parts, extra, mask):
    out = []
    for e in range(len(CYPS)):
        keep = [P[e] for k, P in parts if _keep(e, k) and k != "поферментно"]
        out.append(np.mean(keep + [extra[e]], axis=0))
    return out


def multi_member(preds, seed, honest, mask):
    """Per-enzyme member from oof_multi.json. honest=True -> multi on MULTI_ENZ, независимо on 3A4;
    honest=False -> многозадачно on all four; 'независимо' arm -> the control member."""
    out = []
    for e, c in enumerate(CYPS):
        if honest == "control":
            arm = "независимо"
        elif honest and c not in MULTI_ENZ:
            arm = "независимо"
        else:
            arm = "многозадачно"
        A = np.asarray(preds[f"{seed}|{arm}|{c}"], float)  # already masked, length = mask[:,e].sum()
        out.append(A)
    return out


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

    mp = json.load(open(RES + "preds/oof_multi.json"))["preds"]

    dst = RES + "logs/k90_multitask.json"
    out = json.load(open(dst)) if _pl.Path(dst).exists() else {}
    print(f"  {'сид':>3s} {'арм':>24s} {'ранг':>8s} {'пара':>8s} {'Δранг':>9s} {'Δпара':>9s}",
          flush=True)
    for s in [int(x) for x in a.seeds.split(",") if x.strip()]:
        if f"{s}|многозадачно|CYP1A2" not in mp:
            print(f"  сид {s}: нет ключей в oof_multi.json, пропуск", flush=True)
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
        print(f"  {s:3d} {'эталон (подаётся)':>24s} {base:8.4f} {basep:8.4f} {'':>9s} {'':>9s}",
              flush=True)

        m_honest = multi_member(mp, s, True, mask)
        m_pure = multi_member(mp, s, False, mask)
        m_ctrl = multi_member(mp, s, "control", mask)

        arms = [
            ("standalone multi(honest)", m_honest),
            ("standalone незав(control)", m_ctrl),
            ("INSERT multi(honest)", combine_insert(pdz, m_honest, mask)),
            ("INSERT multi(pure)", combine_insert(pdz, m_pure, mask)),
            ("INSERT незав(control)", combine_insert(pdz, m_ctrl, mask)),
            ("SWAP поф->multi(honest)", combine_swap(pdz, m_honest, mask)),
        ]
        for nm, P in arms:
            sc = score(P, y, LO, HI, mask, fold)
            res[nm] = sc
            standalone = nm.startswith("standalone")
            dr = "" if standalone else f"{sc[0]-base:+9.4f}"
            dp = "" if standalone else f"{sc[1]-basep:+9.4f}"
            print(f"  {s:3d} {nm:>24s} {sc[0]:8.4f} {sc[1]:8.4f} {dr:>9s} {dp:>9s}", flush=True)

        out[str(s)] = res
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(dst, "w"), ensure_ascii=False, indent=1)

    ins_arms = ["INSERT multi(honest)", "INSERT multi(pure)", "INSERT незав(control)",
                "SWAP поф->multi(honest)"]
    print(f"\n  {'арм (над ансамблем)':>24s} {'Δранг ср.':>10s} {'sd':>8s} {'знак>0':>8s} "
          f"{'Δпара ср.':>10s}   поле 0.007", flush=True)
    ks = [k for k in out if "эталон" in out[k]]
    for nm in ins_arms:
        kk = [k for k in ks if nm in out[k]]
        if not kk:
            continue
        d = np.array([out[k][nm][0] - out[k]["эталон"][0] for k in kk])
        q = np.array([out[k][nm][1] - out[k]["эталон"][1] for k in kk])
        verdict = "ПРОХОДИТ" if (d.mean() > 0.007 and (d > 0).sum() == len(d)) else "не проходит"
        print(f"  {nm:>24s} {d.mean():+10.4f} {d.std(ddof=1) if len(d) > 1 else 0:8.4f} "
              f"{int((d > 0).sum()):5d}/{len(d):<2d} {q.mean():+10.4f}   {verdict}", flush=True)

    kk = [k for k in ks if "INSERT multi(honest)" in out[k] and "INSERT незав(control)" in out[k]]
    if kk:
        ch = np.array([out[k]["INSERT multi(honest)"][0] - out[k]["INSERT незав(control)"][0]
                       for k in kk])
        print(f"\n  канал многозадачности над ансамблем (INSERT multi - INSERT незав): "
              f"{ch.mean():+.4f} ранга, знак {int((ch > 0).sum())}/{len(ch)}", flush=True)


if __name__ == "__main__":
    main()
