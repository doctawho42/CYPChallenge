"""Score the fine-tuned encoders `src/finetune.py` wrote, in the PINNED environment.

The split of labour is forced, not stylistic. chemprop resolves `scikit-learn` to 1.9.1 and
CLAUDE.md records that 1.9.0 moves every prediction in this pipeline by ~0.135 pIC50 and doubles
the measured effect of the mechanistic block. So the encoder environment may train and may write an
array, and nothing else: every number that goes in the journal is computed here, under sklearn
1.3.2, by the same `score` and `errcorr` that item 324 used for the probe and items 341/342 used for
the graph nets. The contract between the two halves is a file, as `src/embed.py` established.

WHAT IS COMPARED, and why it is four arms and not two. The request was two checkpoints, with their
difference measuring what CYP pretraining adds. That works only if the two architectures are the
same, and they are not: `rdkit2d` carries d_h 2048 at depth 6 (8.71 M encoder parameters) and
`chemprop_medium` carries d_h 600 at depth 4 (0.82 M) --- a factor of ten. A straight
`medium - rdkit2d` would price CYP pretraining and a tenfold size change together. So each
checkpoint gets its own DISCRIMINATING CONTROL: the identical architecture with the published
weights discarded.

    rdkit2d  -  rand2048     что покупает предобучение на дескрипторах
    medium   -  rand600      что покупает предобучение на CYP/ADME
    rdkit2d  vs medium       сравнимы только через свои контроли, не напрямую

The second difference is the one the request was really about, and it is now free of the size
confound. What remains confounded, and is said rather than hidden: the two checkpoints differ in
pretraining TASK, in width, in depth and in activation, so `medium - rand600` prices "this
particular CYP/ADME checkpoint at its own size", not CYP pretraining in general.

REFERENCE POINTS, all measured on these same folds and already in the journal:

    зонд: гребневая на ЗАМОРОЖЕННОМ chemprop_medium        0.6261   пункт 324
    графовая сеть с нуля, скрининг + четыре сети           0.6079   пункт 342
    графовая сеть с нуля, одна                             0.5778   пункт 341
    слабейший член отгруженного состава (гребневая)        0.5778
    отгруженный состав из шести                            0.6568

The question this file answers is whether letting the encoder MOVE beats leaving it frozen --- that
is, whether any arm clears 0.6261 --- and whether any of them survives as a seventh member.

    uv run python verify/k120_finetune.py

Reads data/ft_*_seed0.npz. Writes results/logs/k120_finetune.json.
Exit: 0 ran, 2 a control failed, 3 no arms found.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json

import numpy as np
import pandas as pd

from cypsplit import butina_folds, fold_digest
from submit import CYPS

import k106_probe as K
from k106_probe import score, errcorr, members, probe_oof, combine_add, load_emb

LOG = RES + "logs/k120_finetune.json"
FLOOR = 0.0036
GOLDEN = "2d93c19815e14261"
# Pairs the differential rests on: (arm, its own random-weight control).
PAIRS = [("rdkit2d", "rand2048"), ("medium", "rand600")]
REF = {"зонд (замороженный, п.324)": 0.6261, "граф+скрининг x4 (п.342)": 0.6079,
       "граф с нуля (п.341)": 0.5778, "слабейший член": 0.5778, "состав из шести": 0.6568}


def load_arm(tag, seed):
    """Prefer the un-scaled file when it exists, and say which one was read.

    `src/finetune.py` built its `RegressionFFN` without an `output_transform`, so chemprop returned
    predictions in each fold's own StandardScaler space. Spearman is invariant to a monotone
    transform WITHIN a fold, but the five folds have different statistics (mu spread 0.012 to 0.039
    across enzymes), so even the rank over the concatenated vector was affected. The transform is
    linear and the split deterministic, so the repair is exact and needs no retraining ---
    `scratchpad/unscale.py` writes `*_unscaled.npz`.
    """
    u = _pl.Path(D + "ft_%s_seed%d_unscaled.npz" % (tag, seed))
    p = u if u.exists() else _pl.Path(D + "ft_%s_seed%d.npz" % (tag, seed))
    if not p.exists():
        return None, None
    if p is not u:
        print("      ВНИМАНИЕ: %s читается БЕЗ разворота масштаба, пара будет бессмысленной" % tag,
              flush=True)
    z = np.load(p, allow_pickle=True)
    dg = str(z["digest"])
    meta = json.loads(str(z["meta"]))
    return z["oof"], {"digest": dg, "meta": meta, "file": p.name}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)

    rows = pd.read_csv(D + "rows.csv")
    z = np.load(D + "feats.npz")
    XFULL = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)
    K.S.LO, K.S.HI = LO, HI

    fold, ncl = butina_folds(list(rows.SMILES), seed=a.seed)
    dg = fold_digest(fold)
    print("  сид %d: дайджест %s, кластеров %d" % (a.seed, dg, ncl), flush=True)
    if a.seed == 0 and dg != GOLDEN:
        print("ОТКАЗ: разбиение сдвинулось, сид 0 дал %s" % dg)
        return 2

    tags = [t for pair in PAIRS for t in pair]
    arms = {}
    for t in tags:
        P, info = load_arm(t, a.seed)
        if P is None:
            print("      %-9s файла нет, пропуск" % t, flush=True)
            continue
        # CONTROL: the other environment must have used OUR folds. Without this the numbers
        # describe a different split and nothing below is comparable to anything.
        if info["digest"] != dg:
            print("ОТКАЗ: %s собран на разбиении %s, а не %s" % (t, info["digest"], dg))
            return 2
        if not np.isfinite(P).all():
            nmiss = int((~np.isfinite(P)).any(1).sum())
            print("ОТКАЗ: %s оставил %d строк незаполненными --- не все фолды прошли"
                  % (t, nmiss))
            return 2
        arms[t] = (P, info)
    if not arms:
        print("ни одного арма не найдено: сначала прогоните src/finetune.py в окружении энкодера")
        return 3

    # ---- the shipped composition and the probe, as every other member measurement does ----
    pdz = members(a.seed, XFULL, Y, M, fold, LO, HI)
    E = load_emb(rows)
    probe = [probe_oof(E, Y, M, fold, e)[0] for e in range(len(CYPS))]
    five = [np.mean([Q[e] for k, Q in pdz if K._keep(e, k)], 0) for e in range(len(CYPS))]
    six = combine_add(pdz, probe)
    s_five, s_six = score(five, Y, LO, HI, M, fold), score(six, Y, LO, HI, M, fold)
    gain = s_six[0] - s_five[0]
    print("      КОНТРОЛЬ состав: пять %.4f -> шесть %.4f, +%.4f ранга (пункт 324: +0.0138)"
          % (s_five[0], s_six[0], gain), flush=True)
    if gain < 0.005:
        print("ОТКАЗ: добавление зонда не воспроизводит выигрыш пункта 324")
        return 2

    out = json.load(open(LOG)) if _pl.Path(LOG).exists() else {}
    out.setdefault("meta", {}).update(floor=FLOOR, refs=REF, pairs=PAIRS)
    res = {}
    for t, (P, info) in arms.items():
        cols = [P[M[:, e], e] for e in range(len(CYPS))]
        sc = score(cols, Y, LO, HI, M, fold)
        _pm, vm, _vs, _cal = errcorr(cols, pdz, Y, M)
        seven = [np.mean([Q[e] for k, Q in pdz if K._keep(e, k)] + [probe[e], cols[e]], 0)
                 for e in range(len(CYPS))]
        s7 = score(seven, Y, LO, HI, M, fold)
        ep = [d["epochs"] for d in info["meta"]["log"]]
        sec = sum(d["seconds"] for d in info["meta"]["log"])
        res[t] = {"solo": sc, "seventh": s7, "gain_rank": s7[0] - s_six[0],
                  "gain_pair": s7[1] - s_six[1], "errcorr_vs_five": vm,
                  "epochs": ep, "seconds": round(sec, 1),
                  "loaded_params": info["meta"]["log"][0].get("loaded_params", 0)}
        print("      %-9s ранг %.4f пара %.4f | corr %.3f | 7-м %+.4f | эпох %s | %.0f с"
              % (t, sc[0], sc[1], float(np.mean(vm)), s7[0] - s_six[0],
                 "/".join(str(x) for x in ep), sec), flush=True)

    out["seeds" if "seeds" in out else "seeds"] = out.get("seeds", {})
    out["seeds"][str(a.seed)] = {"arms": res, "five": s_five, "six": s_six}
    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)

    print("\n" + "=" * 92)
    print("ДООБУЧЕНИЕ ЭНКОДЕРА ПРОТИВ ЗАМОРОЗКИ, сид %d" % a.seed)
    print("=" * 92)
    print("  %-26s %9s %9s %8s %9s" % ("арм", "одиночн.", "пара", "corr", "Δ 7-м"))
    for name, v in REF.items():
        print("  %-26s %9.4f %9s %8s %9s" % (name, v, "-", "-", "-"))
    print("  " + "-" * 64)
    for t in tags:
        if t not in res:
            continue
        v = res[t]
        print("  %-26s %9.4f %9.4f %8.3f %+9.4f"
              % (t, v["solo"][0], v["solo"][1], float(np.mean(v["errcorr_vs_five"])),
                 v["gain_rank"]))
    print("\n  РАЗНОСТИ, каждая против СВОЕГО контроля той же архитектуры:")
    for arm, ctl in PAIRS:
        if arm in res and ctl in res:
            d = res[arm]["solo"][0] - res[ctl]["solo"][0]
            print("    %-10s - %-9s = %+.4f одиночного ранга  (пол %.4f)"
                  % (arm, ctl, d, FLOOR))
        else:
            print("    %-10s - %-9s = нет одного из армов" % (arm, ctl))
    if all(t in res for t in ("medium", "rand600", "rdkit2d", "rand2048")):
        a1 = res["medium"]["solo"][0] - res["rand600"]["solo"][0]
        a2 = res["rdkit2d"]["solo"][0] - res["rand2048"]["solo"][0]
        print("    и разность разностей: CYP/ADME сверх дескрипторного = %+.4f" % (a1 - a2))
        print("    (confounded: чекпойнты различаются и задачей предобучения, и шириной,")
        print("     и глубиной, и активацией --- это цена ЭТОГО чекпойнта, не CYP в общем)")
    print("\n  Решает: превышает ли какой-нибудь арм 0.6261 --- зонд на ЗАМОРОЖЕННОМ энкодере.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
