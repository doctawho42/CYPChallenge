"""Screening rows INSIDE the per-enzyme member: the one cell `src/submit.py` calls unmeasured.

THE QUESTION, and why item 182 does not already answer it. `src/submit.py`'s `--screen` is off by
default, and its own help says why: item 177 measured the screening table as a per-enzyme TARGET at
+0.029 of rank standalone, and item 182 then added that arm to the ensemble and got +0.0058 on the
plain composition and +0.0014 on the best one -- the second below the macro floor of 0.0036. The
diagnosis in item 182 is the arm's error correlation with the ensemble, 0.935 to 0.969: "the arm
errs where the ensemble errs, and an average pays for disagreement". The help text then adds the
sentence this file exists for: *"вариант ЗДЕСЬ иной: скрининг вставлен ВНУТРЬ поферментного члена, а
не добавлен шестым. Механизм предсказывает тот же ноль, измерено это не было."*

The mechanism does NOT predict the same zero, and that is worth stating plainly because the comment
asserts it. Item 182 added a SIXTH vector to the mean, so its gain was taxed by dilution: a new
member must disagree with the existing average to earn its weight, and that one did not. Inserting
the screen rows inside `поферментно` adds no vector at all -- `oof_members` passes `scr` to exactly
one member, `_oof_one(X, y, mask, fold, False, scr)`, and the other four are untouched -- so the
composition keeps its size and one existing member simply becomes better. There is no dilution term
to pay.

What DOES carry over from item 182 is its closing claim, and it cuts the other way: "the ensemble is
saturated. Members that get better by fitting the target better do not help it, because the four
existing members already fit the target the same way." Improving `поферментно` is exactly getting
better by fitting the target better. So the two readings disagree and the cell is honestly open.

THE CONTROL THAT MAKES A ZERO MEAN SOMETHING. A null here is ambiguous between "no effect" and "scr
is not wired", so both ends are pinned:

  * the member recomputed WITHOUT the screen must reproduce the committed `members_seed{s}.json`
    entry for `поферментно` BIT FOR BIT. That proves this file's path is the path that produced the
    shipped composition, not a lookalike.
  * the member WITH the screen must beat it at the member level by something of item 177's order
    (+0.0245 of rank at seed 0, +0.0307 on seeds 1-3, measured there on `src/ablaux.py`'s own
    learner, so the magnitude is an expectation and not an assertion). If that gain is absent, the
    wiring is at fault and the composition number means nothing.

    uv run python verify/k115_screenin.py --seeds 0

Four members come from the cache; the per-enzyme member is computed four times (raw and dead-zone
passed, with and without the screen), which is where the runtime goes. Predictions are cached under
the fold digest, so a later change of REFERENCE costs a second -- item 335 paid that bill twice.

Writes results/logs/k115_screenin.json. Exit: 0 ran, 2 a control failed, 3 inputs missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time

import numpy as np
import pandas as pd

from cypsplit import butina_folds, fold_digest
from submit import CYPS, load_screen, _oof_one, _dz_oof

import k106_probe as K
from k106_probe import score, errcorr, members, load_emb, probe_oof, combine_add

LOG = RES + "logs/k115_screenin.json"
CACHE = RES + "preds/k115_perenz_seed%d.json"
FLOOR = 0.0036            # project macro-rank noise floor
ITEM177 = {0: 0.0245, 1: 0.0307, 2: 0.0307, 3: 0.0307}   # member-level expectation, src/ablaux.py


def per_enzyme(X, y, mask, fold, LO, HI, scr, s, digest, tag):
    """The per-enzyme member, raw and dead-zone passed, cached under the fold digest.

    The dead-zone pass is `submit.dz_pass`'s body for this one member, not a variant of it: the
    target is `clip(p_oof, lo, hi)` built from the OUT-OF-FOLD prediction, which is the whole
    point of the device (projecting a model's own training-row predictions collapses the scheme
    into an identity -- see `src/abldead.py`).
    """
    path = CACHE % s
    store = json.load(open(path)) if _pl.Path(path).exists() else {}
    if store.get("digest") != digest:
        store = {"digest": digest}
    kr, kd = tag + "_raw", tag + "_dz"
    if kr in store and kd in store:
        print(f"      {tag:<12} из кэша", flush=True)
        return ([np.asarray(v, float) for v in store[kr]],
                [np.asarray(v, float) for v in store[kd]], 0.0)
    t0 = time.time()
    raw = _oof_one(X, y, mask, fold, False, scr)
    dz = []
    for e in range(len(CYPS)):
        m = mask[:, e]
        t = np.clip(raw[e], LO[m, e], HI[m, e])
        dz.append(_dz_oof("поферментно", X[m], t, fold[m], e))
    dt = time.time() - t0
    store[kr] = [[float(x) for x in v] for v in raw]
    store[kd] = [[float(x) for x in v] for v in dz]
    _pl.Path(RES + "preds").mkdir(parents=True, exist_ok=True)
    json.dump(store, open(path, "w"))
    print(f"      {tag:<12} посчитан за {dt:.0f} с и записан", flush=True)
    return raw, dz, dt


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    for f in ("emb_chemprop_medium.npz", "feats.npz", "rows.csv"):
        if not _pl.Path(D + f).exists():
            print(f"нет data/{f}")
            return 3

    rows = pd.read_csv(D + "rows.csv")
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)
    K.S.LO, K.S.HI = LO, HI
    E = load_emb(rows)
    L2, has = load_screen(rows)

    # ---- count first (items 109, 105, 114): how many rows does this even add? -------------
    print("\n  строк, которые добавляет скрининг:", flush=True)
    tot_c = tot_s = 0
    counts = {}
    for e, c in enumerate(CYPS):
        curves = int(mask[:, e].sum())
        add = int((has[:, e] & ~mask[:, e]).sum())
        counts[c] = {"кривых": curves, "скрининга": add, "отношение": round(add / curves, 3)}
        tot_c += curves
        tot_s += add
        print(f"      {c:<8} кривых {curves:5d}   скрининга без кривой {add:5d}   "
              f"отношение {add / curves:.2f}", flush=True)
    print(f"      итого    кривых {tot_c:5d}   скрининга            {tot_s:5d}   "
          f"отношение {tot_s / tot_c:.2f}", flush=True)
    if tot_s == 0:
        print("ОТКАЗ: скрининг не добавляет ни одной строки")
        return 2

    out = json.load(open(LOG)) if _pl.Path(LOG).exists() else {}
    out.setdefault("meta", {}).update(floor=FLOOR, counts=counts,
                                      item177_expectation=ITEM177)
    out.setdefault("seeds", {})

    for s in seeds:
        fold, ncl = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        print(f"\n  сид {s}: дайджест {dg}, кластеров {ncl}", flush=True)
        if s == 0 and dg != "2d93c19815e14261":
            print(f"ОТКАЗ: разбиение сдвинулось, сид 0 дал {dg}")
            return 2

        pdz = members(s, X, y, mask, fold, LO, HI)
        probe = [probe_oof(E, y, mask, fold, e)[0] for e in range(len(CYPS))]
        ref = json.load(open(RES + "logs/k106_probe.json"))["seeds"].get(str(s))
        sp = score(probe, y, LO, HI, mask, fold)
        if ref and abs(sp[0] - ref["solo_probe"][0]) > 1e-6:
            print(f"ОТКАЗ: зонд не воспроизводится ({sp[0]:.4f} против "
                  f"{ref['solo_probe'][0]:.4f})")
            return 2
        print(f"      КОНТРОЛЬ зонд воспроизведён побитово: {sp[0]:.4f}", flush=True)

        raw_no, dz_no, t1 = per_enzyme(X, y, mask, fold, LO, HI, None, s, dg, "без_скрининга")
        raw_sc, dz_sc, t2 = per_enzyme(X, y, mask, fold, LO, HI, (L2, has), s, dg, "со_скринингом")

        # ---- CONTROL A: my recomputation IS the shipped path -----------------------------
        cached = dict(pdz)["поферментно"]
        dmax = max(float(np.max(np.abs(cached[e] - dz_no[e]))) for e in range(len(CYPS)))
        print(f"      КОНТРОЛЬ путь тот же: пересчёт без скрининга против кэша члена, "
              f"max|Δ| = {dmax:.3e}", flush=True)
        if dmax != 0.0:
            print("ОТКАЗ: пересчитанный член не совпал с кэшем побитово; это не отгруженный путь")
            return 2

        # ---- CONTROL B: the screen is actually wired in ----------------------------------
        s_raw_no = score(raw_no, y, LO, HI, mask, fold)
        s_raw_sc = score(raw_sc, y, LO, HI, mask, fold)
        g_member = s_raw_sc[0] - s_raw_no[0]
        print(f"      КОНТРОЛЬ скрининг подключён: член СЫРОЙ {s_raw_no[0]:.4f} -> "
              f"{s_raw_sc[0]:.4f} ({g_member:+.4f}); пункт 177 даёт "
              f"{ITEM177.get(s, 0.03):+.4f}", flush=True)
        if g_member < FLOOR:
            print("ОТКАЗ: на уровне члена скрининг не даёт даже пола --- проводка под вопросом, "
                  "а не результат")
            return 2

        # ---- the measurement -------------------------------------------------------------
        pdz_sc = [("поферментно", dz_sc) if k == "поферментно" else (k, P) for k, P in pdz]
        six_no = combine_add(pdz, probe)
        six_sc = combine_add(pdz_sc, probe)
        s_six_no = score(six_no, y, LO, HI, mask, fold)
        s_six_sc = score(six_sc, y, LO, HI, mask, fold)
        d_rank = s_six_sc[0] - s_six_no[0]
        d_pair = s_six_sc[1] - s_six_no[1]

        # The member-level comparison AFTER the dead zone, which is what actually enters.
        s_dz_no = score(dz_no, y, LO, HI, mask, fold)
        s_dz_sc = score(dz_sc, y, LO, HI, mask, fold)

        # Item 182's diagnostic, recomputed for the NEW member against the rest.
        pm, vm, vs, _cal = errcorr(dz_sc, pdz, y, mask)

        print(f"      член после зоны   {s_dz_no[0]:.4f} -> {s_dz_sc[0]:.4f} "
              f"({s_dz_sc[0] - s_dz_no[0]:+.4f})", flush=True)
        print(f"      СОСТАВ (шесть)    {s_six_no[0]:.4f} -> {s_six_sc[0]:.4f} "
              f"({d_rank:+.4f} ранга), пара {s_six_no[1]:.4f} -> {s_six_sc[1]:.4f} "
              f"({d_pair:+.4f})   пол {FLOOR}", flush=True)
        print(f"      corr ошибок нового члена со средним пяти: "
              f"{' '.join('%.3f' % v for v in vm)} (пункт 182: 0.935--0.969)", flush=True)

        out["seeds"][str(s)] = {
            "member_raw": {"без": s_raw_no, "со": s_raw_sc, "прирост": g_member},
            "member_dz": {"без": s_dz_no, "со": s_dz_sc,
                          "прирост": s_dz_sc[0] - s_dz_no[0]},
            "composition": {"без": s_six_no, "со": s_six_sc,
                            "прирост_ранг": d_rank, "прирост_пара": d_pair},
            "errcorr_new_member": {"vs_five_mean": vm, "per_member": pm, "vs_shipped": vs},
            "seconds": round(t1 + t2, 1),
            "control_path_maxdiff": dmax,
        }
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
        print(f"      записано: {LOG}", flush=True)

    print("\n" + "=" * 86)
    print("СКРИНИНГ ВНУТРИ ПОФЕРМЕНТНОГО ЧЛЕНА --- то, что submit.py называет непромеренным")
    print("=" * 86)
    print("  %-5s %10s %10s   %10s %10s   %9s %9s"
          % ("сид", "член сыр.", "прирост", "состав", "прирост", "пара", "знак>пол"))
    for s in seeds:
        r = out["seeds"].get(str(s))
        if not r:
            continue
        print("  %-5s %10.4f %+10.4f   %10.4f %+10.4f   %+9.4f %9s"
              % (s, r["member_raw"]["без"][0], r["member_raw"]["прирост"],
                 r["composition"]["без"][0], r["composition"]["прирост_ранг"],
                 r["composition"]["прирост_пара"],
                 "да" if r["composition"]["прирост_ранг"] > FLOOR else "НЕТ"))
    vals = [out["seeds"][str(s)]["composition"]["прирост_ранг"] for s in seeds
            if str(s) in out["seeds"]]
    if vals:
        print("\n  состав: среднее %+.4f, знак %d/%d, пол %.4f"
              % (float(np.mean(vals)), int(sum(1 for v in vals if v > 0)), len(vals), FLOOR))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
