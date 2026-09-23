"""Item 318's arm, run: the three-head trunk REPLACING the trunk member in place, by rank
over the shipped ensemble, against the pre-registered gate of +0.0105 macro Spearman.

Three arms, scored over the ENSEMBLE, not over the member:

  A  REFERENCE : the shipped five-member ensemble, ствол = results/preds/trunk_twohead_dead.json.
  B  CONTROL   : the same ensemble with ствол replaced by a THREE-head trunk trained with the
                 8004 external rows PRESENT and lam_ext = 0. This separates the cost of
                 inserting the rows (standardisation and step count -- trunkext.py's docstring
                 measured that at 6.49 pIC50 maximum, it is not nothing) from the value of the
                 channel. Item 290 paid -0.0189 for insertion alone; without this arm the
                 result is uninterpretable.
  C  ARM       : the same, with lam_ext > 0.

Everything else about the replaced member is the shipped member's: mode twohead, lam 3.0,
BLOCKS DESC+MECH, the dead-zone pass with the target projected from trunk_twohead.json at the
SAME seed, l1 on the pIC50 head, and _trunk_clip. Only the external rows and the third head
differ, which is what makes B - A the vehicle's cost and C - B the channel.

lam_ext is NOT searched -- see the LAM_EXT comment below.

Resumable by construction, because the previous attempt at this run lost 90 minutes to a
network drop with nothing written down: fold-level predictions land in
results/preds/trunk_thirdhead.json the moment each fold finishes, and each seed's scores land
in results/logs/k103_thirdhead.json the moment that seed finishes. Re-running skips whatever
is already there.
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

from cypsplit import butina_folds, fold_digest
from shrinkchoice import fit_apply
import submit as S
from submit import (CYPS, oof_members, dz_pass, _keep, _trunk_clip, _oof_trunk,
                    TRUNK_LAM, TRUNK_MODE)
import trunk as T

# lam_ext = 3.0, carried in, NOT chosen here. The discipline is src/trunk.py's: everything the
# two arms SHARE was fixed on the channel-off arm so the choice cannot favour the experiment,
# and the lambda itself was never tuned against the outcome at all -- submit.TRUNK_LAM's own
# comment reads "значение, на котором пункт 79 мерил канал", i.e. it is the value the channel
# was already measured at, not a value picked to make it look good.
#
# Three independent reasons 3.0 is that value here:
#   1. trunk.py states that lam_ext weights its term through the SAME masked_mse, normalised by
#      observed cells, so lam and lam_ext are the same kind of number. The shipped screen head
#      runs at 3.0.
#   2. src/trunkext.py's default is --lams 0,3.0; items 66 and 290 measured this very external
#      channel at 3.0, in the other vehicle.
#   3. Nothing else in the tree prices it.
#
# NO GRID IS RUN. A grid over lam_ext could only be scored on the channel-ON arm, which is the
# definition of selection contamination (item 245); one arm at a carried-in value is the honest
# form, and the cost of that honesty is that this measures lam_ext = 3.0 and not the best
# lam_ext.
LAM_EXT = 3.0

PRED = RES + "preds/trunk_thirdhead.json"
LOG = RES + "logs/k103_thirdhead.json"
GATE = 0.0105


def score(P, y, LO, HI, mask, fold):
    """verify/k83_loo.py's and k89's score(), unchanged: rank = mean Spearman
    (affine-pair-invariant), pair = mean ST-RAE AFTER the per-fold affine pair."""
    rk, pr = [], []
    for e in range(len(CYPS)):
        m = mask[:, e]
        yy, lo, hi = y[m, e], LO[m, e], HI[m, e]
        u = np.ones(len(yy)) / len(yy)
        rk.append(float(spearmanr(yy, P[e]).statistic))
        pr.append(float(strae(yy, fit_apply(P[e], lo, hi, fold[m], u),
                              y_true_upper=hi, y_true_lower=lo)))
    return float(np.mean(rk)), float(np.mean(pr)), rk, pr


def combine_keep(parts):
    return [np.mean([P[e] for k, P in parts if _keep(e, k)], axis=0) for e in range(len(CYPS))]


def combine_replace(parts, new):
    """The shipped mean with ствол REPLACED by `new`. Not a sixth member: item 316 measured a
    sixth member at any fitted weight as exactly +0.0000, so adding is dead and replacing is the
    untested cell."""
    out = []
    for e in range(len(CYPS)):
        keep = [P[e] for k, P in parts if _keep(e, k) and k != "ствол"]
        if _keep(e, "ствол"):
            keep = keep + [new[e]]
        out.append(np.mean(keep, axis=0))
    return out


MEMBERS = ("поферментно", "пул", "GP", "гребневая", "ствол")


def members(s, XFULL, y, mask, fold, LO, HI):
    """The five dead-zone-passed members, from the committed cache when it exists.

    results/preds/members_seed{s}.json is written by verify/k84_subsets.py as exactly
    dz_pass(oof_members(..., "ансамбль5", seed=s, dead=True)), which is what this needs, and
    rebuilding it costs about fifteen minutes a seed.

    The control that makes the cache usable has to exercise the clause that could be wrong --
    "the cache was built on a different split or a different trunk file" -- rather than merely
    prove the file parses. The cached ствол is the one member that is READ rather than trained,
    so it can be re-derived here in milliseconds from the digest-guarded trunk file for THIS
    seed, and it must match to the bit. A cache from another split fails this.
    """
    p = RES + f"preds/members_seed{s}.json"
    if not _pl.Path(p).exists():
        parts = oof_members(XFULL, y, mask, fold, "ансамбль5", seed=s, dead=True)
        return dz_pass(parts, XFULL, mask, fold, (LO, HI))
    mem = json.load(open(p))
    if set(mem) != set(MEMBERS):
        raise SystemExit(f"кэш членов сида {s} содержит {sorted(mem)}, а не {sorted(MEMBERS)}")
    pdz = [(k, [np.asarray(mem[k][e], float) for e in range(len(CYPS))]) for k in MEMBERS]
    fresh = _oof_trunk(y, mask, fold, s, dead=True)
    cached = dict(pdz)["ствол"]
    d = max(float(np.max(np.abs(cached[e] - fresh[e]))) for e in range(len(CYPS)))
    if d != 0.0:
        raise SystemExit(f"сид {s}: кэш членов расходится со стволом на {d:.3e}; пересоберите")
    print(f"      члены из кэша {p} (контроль: ствол совпал побитово, {d:.1e})", flush=True)
    return pdz


def load_pred_cache():
    return json.load(open(PRED)) if _pl.Path(PRED).exists() else {}


def train_arm(key, XX, yy, ss, ee, ff, tgt, lam_ext, seed, device, n):
    """Five folds of the three-head trunk, checkpointed after EVERY fold."""
    cache = load_pred_cache()
    rec = cache.get(key) or {"pred": [[None] * 4 for _ in range(n)], "done": []}
    pred = np.array([[np.nan if v is None else v for v in r] for r in rec["pred"]], float)
    done = set(int(f) for f in rec["done"])
    for f in range(5):
        if f in done:
            continue
        t0 = time.time()
        te = ff == f
        if te.sum() == 0:
            done.add(f)
            continue
        p = T.run_fold(XX, yy, ss, ff, f, float(TRUNK_LAM), seed, device, TRUNK_MODE,
                       target=tgt, l1=True, ext=ee, lam_ext=lam_ext)
        # External rows carry fold index -1, so `te` is False for every one of them and the
        # held-out block lies wholly inside the first n rows. Asserted, not assumed.
        if te.sum() != te[:n].sum():
            raise SystemExit("внешняя строка попала в отложенный фолд")
        pred[te[:n]] = p
        done.add(f)
        cache = load_pred_cache()
        cache[key] = {"pred": [[None if not np.isfinite(v) else float(v) for v in r]
                               for r in pred],
                      "done": sorted(done)}
        json.dump(cache, open(PRED, "w"))
        print(f"      фолд {f}: {time.time()-t0:.0f} c, записано в {PRED}", flush=True)
    return pred


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--device", default="mps")
    ap.add_argument("--lam-ext", type=float, default=LAM_EXT)
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    # ---- data ------------------------------------------------------------------------------
    X, y32, lo32, hi32, scr, smiles = T.load(T.BLOCKS)
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

    E = np.load(D + "ext_thirdhead.npz")
    Xe, ye = E["Xe"], E["ye"]
    n, m = len(X), len(Xe)
    if Xe.shape[1] != X.shape[1]:
        raise SystemExit(f"ширина не совпала: {Xe.shape[1]} против {X.shape[1]}")
    XX = np.vstack([X, Xe])
    yy = np.vstack([y32, np.full_like(ye, np.nan)])
    # THE SCREEN STAYS. This is the whole difference from trunkext.py, which put the external
    # labels in this block and so DROPPED the screen channel -- the reason its trunk scores
    # 0.5270 against the shipped 0.5966.
    ss = np.vstack([scr, np.full_like(ye, np.nan)])
    ee = np.vstack([np.full_like(scr, np.nan), ye])
    print(f"строк: наших {n}, внешних {m}, всего {len(XX)}", flush=True)
    print("внешних меток: " + ", ".join(f"{c} {int(np.isfinite(ye[:, e]).sum())}"
                                        for e, c in enumerate(CYPS)), flush=True)

    base_src = json.load(open(RES + "preds/trunk_twohead.json"))
    base_preds, base_meta = base_src["preds"], base_src.get("meta", {})

    out = json.load(open(LOG)) if _pl.Path(LOG).exists() else {}
    out.setdefault("meta", {}).update(
        lam_ext=a.lam_ext, lam=float(TRUNK_LAM), mode=TRUNK_MODE, blocks=T.BLOCKS,
        device=a.device, gate=GATE, n_ours=n, n_ext=m,
        arms={"A": "shipped ensemble, two-head trunk",
              "B": "three-head trunk in place, ext rows present, lam_ext=0",
              "C": f"three-head trunk in place, ext rows present, lam_ext={a.lam_ext}"})
    out.setdefault("seeds", {})

    hdr = f"  {'сид':>3s} {'арм':>26s} {'ранг':>8s} {'пара':>8s} {'Δранг':>9s} {'Δпара':>9s}"
    print(hdr, flush=True)
    for s in seeds:
        if str(s) in out["seeds"] and all(k in out["seeds"][str(s)] for k in "ABC"):
            r = out["seeds"][str(s)]
            print(f"  сид {s}: уже в {LOG}, пропуск "
                  f"(A {r['A'][0]:.4f} B {r['B'][0]:.4f} C {r['C'][0]:.4f})", flush=True)
            continue

        fold, _ = butina_folds(smiles, seed=s)
        d = fold_digest(fold)
        # Guard: the dead-zone target and the reference member must both have been computed on
        # THESE folds, or the target is not out-of-fold and the member is not comparable.
        for nm, src in (("trunk_twohead.json", base_meta),
                        ("trunk_twohead_dead.json",
                         json.load(open(RES + "preds/trunk_twohead_dead.json")).get("meta", {}))):
            rec = (src or {}).get("fold_digests") or {}
            if rec.get(str(s)) != d:
                raise SystemExit(f"сид {s}: дайджест {d}, а {nm} посчитан на {rec.get(str(s))}")

        key0 = f"{TRUNK_MODE}|{s}|{float(TRUNK_LAM)}"
        A0 = np.asarray(base_preds[key0], float)
        if not np.array_equal(np.isnan(A0), np.isnan(y)):
            raise SystemExit("маска NaN мишени не совпадает с маской меток")
        tgt = np.vstack([np.clip(A0, lo32, hi32).astype(y32.dtype),
                         np.full_like(ye, np.nan)])
        ff = np.concatenate([fold, np.full(m, -1, dtype=fold.dtype)])

        t0 = time.time()
        pB = train_arm(f"{s}|0.0", XX, yy, ss, ee, ff, tgt, 0.0, s, a.device, n)
        pC = train_arm(f"{s}|{a.lam_ext}", XX, yy, ss, ee, ff, tgt, a.lam_ext, s, a.device, n)
        print(f"    обучение сида {s}: {time.time()-t0:.0f} c", flush=True)

        # Controls that must succeed AND exercise the suspect clause -------------------------
        # 1. the external rows really changed the model (else arm B is arm A by another name);
        # 2. the channel really turned on (else arm C is arm B by another name).
        dBA = float(np.nanmax(np.abs(pB - A0)))
        dCB = float(np.nanmax(np.abs(pC - pB)))
        print(f"    контроль: |B - двухголовый| max {dBA:.3f} (>0 ожидается), "
              f"|C - B| max {dCB:.3f} (>0 ожидается)", flush=True)
        if dBA == 0.0 or dCB == 0.0:
            raise SystemExit("контроль не сработал: одна из рук совпала с другой побитово")

        t0 = time.time()
        pdz = members(s, XFULL, y, mask, fold, LO, HI)
        print(f"      члены: {time.time()-t0:.0f} c", flush=True)

        mem = lambda P: [_trunk_clip(P[mask[:, e], e], y[mask[:, e], e]) for e in range(len(CYPS))]
        res = out["seeds"].get(str(s), {})
        res["A"] = score(combine_keep(pdz), y, LO, HI, mask, fold)
        res["B"] = score(combine_replace(pdz, mem(pB)), y, LO, HI, mask, fold)
        res["C"] = score(combine_replace(pdz, mem(pC)), y, LO, HI, mask, fold)
        res["solo_ствол"] = score([P for k, P in pdz if k == "ствол"][0], y, LO, HI, mask, fold)
        res["solo_B"] = score(mem(pB), y, LO, HI, mask, fold)
        res["solo_C"] = score(mem(pC), y, LO, HI, mask, fold)
        res["controls"] = {"max_abs_B_minus_twohead": dBA, "max_abs_C_minus_B": dCB}

        base = res["A"][0]
        basep = res["A"][1]
        for nm in ["A", "B", "C", "solo_ствол", "solo_B", "solo_C"]:
            sc = res[nm]
            solo = nm.startswith("solo")
            dr = "" if solo else f"{sc[0]-base:+9.4f}"
            dp = "" if solo else f"{sc[1]-basep:+9.4f}"
            print(f"  {s:3d} {nm:>26s} {sc[0]:8.4f} {sc[1]:8.4f} {dr:>9s} {dp:>9s}", flush=True)

        out["seeds"][str(s)] = res
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
        print(f"    записано в {LOG}", flush=True)

    report(out)


def report(out):
    ks = sorted(out["seeds"], key=int)
    ks = [k for k in ks if all(x in out["seeds"][k] for x in "ABC")]
    if not ks:
        return
    print(f"\n  сидов: {len(ks)} ({', '.join(ks)})", flush=True)
    print(f"\n  {'арм':>34s} {'ранг ср.':>9s} {'ST-RAE ср.':>11s}", flush=True)
    for nm, lab in [("A", "A эталон (двухголовый ствол)"),
                    ("B", "B три головы, lam_ext=0"),
                    ("C", f"C три головы, lam_ext={out['meta']['lam_ext']}")]:
        r = np.array([out["seeds"][k][nm][0] for k in ks])
        p = np.array([out["seeds"][k][nm][1] for k in ks])
        print(f"  {lab:>34s} {r.mean():9.4f} {p.mean():11.4f}", flush=True)

    print(f"\n  {'разность':>34s} {'Δранг ср.':>10s} {'sd':>8s} {'знак>0':>8s} "
          f"{'Δпара ср.':>10s}   ворота {GATE}", flush=True)
    for lab, a_, b_ in [("C - A  (РЕШАЮЩЕЕ)", "C", "A"),
                        ("C - B  (канал)", "C", "B"),
                        ("B - A  (стоимость вставки строк)", "B", "A")]:
        d = np.array([out["seeds"][k][a_][0] - out["seeds"][k][b_][0] for k in ks])
        q = np.array([out["seeds"][k][a_][1] - out["seeds"][k][b_][1] for k in ks])
        sd = d.std(ddof=1) if len(d) > 1 else 0.0
        v = ""
        if lab.startswith("C - A"):
            v = ("ПРОХОДИТ" if (d.mean() >= GATE and (d > 0).sum() >= 3 and len(d) >= 4)
                 else "не проходит")
        print(f"  {lab:>34s} {d.mean():+10.4f} {sd:8.4f} {int((d > 0).sum()):5d}/{len(d):<2d} "
              f"{q.mean():+10.4f}   {v}", flush=True)

    print(f"\n  {'соло (член сам по себе)':>34s} {'ранг ср.':>9s}", flush=True)
    for nm, lab in [("solo_ствол", "ствол двухголовый (подаётся)"),
                    ("solo_B", "три головы lam_ext=0"),
                    ("solo_C", f"три головы lam_ext={out['meta']['lam_ext']}")]:
        kk = [k for k in ks if nm in out["seeds"][k]]
        if kk:
            r = np.array([out["seeds"][k][nm][0] for k in kk])
            print(f"  {lab:>34s} {r.mean():9.4f}", flush=True)


if __name__ == "__main__":
    main()
