"""Does a NON-LINEAR head on the same 600 frozen dimensions decorrelate enough to be a member?

THE QUESTION, and why it is worth minutes rather than a dependency. Item 300 declares the
tabular foundation-model axis -- TabPFN, TabICL, Mitra -- the one axis this project has never
touched, and a re-grep 34 items later (313 items, 22..334) confirms it: `Mitra` and `autogluon`
return zero against controls `chemprop` 12 and `member` 429 that must and do hit, and every
TabPFN/TabICL hit resolves inside item 300 itself. In the code, every head ever fitted on
`data/emb_chemprop_medium.npz` is linear -- `RidgeCV` in `verify/k106_probe.py:probe_oof` and in
`src/submit.py:_oof_probe`. So the axis is open.

But a foundation head is a dependency and a build. Before paying for it, the SAME question can be
asked of a cheap non-linear head on the SAME columns, because what decides a new member is not its
architecture but whether its out-of-fold errors differ from what the ensemble already makes. That
is item 289's screen, and `k106_probe.py:errcorr` already computes it. This file swaps the head and
changes nothing else.

WHY THE GATE IS A LIVE INSTRUMENT AND NOT A SLOGAN. Read off the committed
`results/logs/k106_probe.json`, 8 seeds, correlation of out-of-fold ERROR against the five-member
mean:

    the five members themselves     0.966 -- 0.982      redundant by construction
    the SHIPPED probe               0.901               +0.0138 macro rank, sign 8/8, items 324-327
    item 289's gate                 0.93

The gate sits between the band the members occupy and the value the one member that ever cleared it
actually has. So 0.93 discriminates on these exact data, and this file re-derives both ends of that
spread as its controls rather than quoting them.

THE REFERENCE THE PROBE ITSELF DOES NOT NEED, and this file adds. The probe SHIPPED (items 325/327;
every column of the uploaded `activity_submission_best4.csv` descends from
`activity_submission_probe_base.csv`). A new head on the same 600 columns would therefore be a
SEVENTH member competing with a composition that already contains a ridge on that representation.
Correlation against the five-member mean is no longer the whole question: the head must also differ
from the PROBE. `vs_probe` is reported beside `vs_five_mean`, and the honest gate for this cell is
the larger of the two.

A caution the data carries and this file does not hide: the shipped probe's own correlation is
0.954--0.973 on CYP2D6 in six of eight seeds -- inside the member band -- and the ensemble gained
anyway. The gate is a macro screen, not a per-enzyme veto.

    uv run python verify/k114_gate93.py --seeds 0 --heads nystroem,mlp
    uv run python verify/k114_gate93.py --seeds 0 --heads gbm          # time it alone first

Writes results/logs/k114_gate93.json. Exit: 0 ran, 2 a control failed, 3 inputs missing.
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
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.kernel_approximation import Nystroem
from sklearn.linear_model import RidgeCV
from sklearn.neural_network import MLPRegressor

from cypsplit import butina_folds, fold_digest
from submit import CYPS

# The apparatus, imported rather than re-written: same folds, same score, same gate function,
# same member cache, same embedding loader. Re-implementing any of these would make the numbers
# incomparable to item 324's, which is the only thing this measurement is for.
import k106_probe as K
from k106_probe import (score, errcorr, members, load_emb, probe_oof,
                        combine_add, combine_keep)

LOG = RES + "logs/k114_gate93.json"
CACHE = RES + "preds/k114_heads_seed%d.json"
GATE_CORR = 0.93          # item 289's error-correlation screen
GATE_RANK_SLACK = 0.02    # item 289: standalone rank within 0.02 of the WEAKEST current member
WEAKEST = 0.5778          # `гребневая`, mean over 8 seeds of k106_probe.json's solo_гребневая
PROBE_SOLO = 0.6261       # the shipped probe, mean over the same 8 seeds
PROBE_CORR = 0.901        # the shipped probe's own vs_five_mean, mean over 32 cells
BAND_LO = 0.966           # the loosest of the five members (`ствол`), seed 0


def head_oof(E, y, mask, fold, e, kind, seed=0):
    """A non-linear head on the frozen embedding, out of fold, one enzyme.

    The per-fold standardisation and dead-column filter are lifted verbatim from
    `k106_probe.probe_oof`, so the only difference between this and the shipped probe is the
    estimator. Fitting the scaler on the union would be item 202's defect 2.
    """
    m = mask[:, e]
    Xm, ym, fm = E[m], y[m, e], fold[m]
    p = np.full(int(m.sum()), np.nan)
    for f in range(5):
        trn, te = fm != f, fm == f
        if te.sum() == 0:
            continue
        mu, sd = Xm[trn].mean(0), Xm[trn].std(0)
        keep = sd > 1e-6
        Z = (Xm[:, keep] - mu[keep]) / sd[keep]
        if kind == "nystroem":
            # RBF feature map then the SAME alpha grid as the probe: the smoothest possible
            # departure from linearity, at the probe's cost.
            ny = Nystroem(gamma=1.0 / Z.shape[1], n_components=600,
                          random_state=seed).fit(Z[trn])
            A, B = ny.transform(Z[trn]), ny.transform(Z[te])
            p[te] = RidgeCV(alphas=K.ALPHAS).fit(A, ym[trn]).predict(B)
        elif kind == "mlp":
            p[te] = MLPRegressor(hidden_layer_sizes=(256,), alpha=1e-3, max_iter=400,
                                 early_stopping=True, n_iter_no_change=15,
                                 random_state=seed).fit(Z[trn], ym[trn]).predict(Z[te])
        elif kind == "gbm":
            p[te] = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.06,
                                                  random_state=seed
                                                  ).fit(Z[trn], ym[trn]).predict(Z[te])
        else:
            raise SystemExit(f"неизвестная голова: {kind}")
    if not np.isfinite(p).all():
        raise SystemExit(f"{CYPS[e]}: голова {kind} оставила не заполненные строки")
    return p


def head_cached(E, y, mask, fold, kind, s, digest):
    """`head_oof` for all four enzymes, cached on disk under the fold digest.

    The reference a candidate member is judged against changed once already in this file --- from
    the five-member mean to the shipped six --- and re-deriving the heads to answer a question
    about the REFERENCE is paying twice for one fit. Predictions are valid only for the fold
    vector they were produced on, so the digest is stored beside them and a mismatch recomputes
    rather than silently returning the wrong vectors.
    """
    path = CACHE % s
    store = json.load(open(path)) if _pl.Path(path).exists() else {}
    if store.get("digest") == digest and kind in store.get("heads", {}):
        P = [np.asarray(v, float) for v in store["heads"][kind]]
        if all(len(P[e]) == int(mask[:, e].sum()) for e in range(len(CYPS))):
            print(f"      {kind:<9} из кэша ({path.split('/')[-1]})", flush=True)
            return P, 0.0
    t0 = time.time()
    P = [head_oof(E, y, mask, fold, e, kind, seed=s) for e in range(len(CYPS))]
    dt = time.time() - t0
    if store.get("digest") != digest:
        store = {"digest": digest, "heads": {}}
    store.setdefault("heads", {})[kind] = [[float(x) for x in v] for v in P]
    _pl.Path(RES + "preds").mkdir(parents=True, exist_ok=True)
    json.dump(store, open(path, "w"))
    return P, dt


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--heads", default="nystroem,mlp")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]
    heads = [h.strip() for h in a.heads.split(",") if h.strip()]

    if not _pl.Path(D + "emb_chemprop_medium.npz").exists():
        print("нет data/emb_chemprop_medium.npz (src/embed_medium.py)")
        return 3

    rows = pd.read_csv(D + "rows.csv")
    z = np.load(D + "feats.npz")
    XFULL = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)
    K.S.LO, K.S.HI = LO, HI          # score() reads the bands off the module
    E = load_emb(rows)

    out = json.load(open(LOG)) if _pl.Path(LOG).exists() else {}
    out.setdefault("meta", {}).update(gate_corr=GATE_CORR, weakest_member=WEAKEST,
                                      probe_solo=PROBE_SOLO, probe_corr=PROBE_CORR,
                                      band_loosest_member=BAND_LO, heads=heads)
    out.setdefault("seeds", {})

    for s in seeds:
        fold, ncl = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        print(f"\n  сид {s}: дайджест {dg}, кластеров {ncl}", flush=True)
        if s == 0 and dg != "2d93c19815e14261":
            print(f"ОТКАЗ: разбиение сдвинулось, сид 0 дал {dg}")
            return 2

        # ---- the linear probe, recomputed, and the member cache -------------------------
        t0 = time.time()
        probe = [probe_oof(E, y, mask, fold, e)[0] for e in range(len(CYPS))]
        t_probe = time.time() - t0
        pdz = members(s, XFULL, y, mask, fold, LO, HI)
        six = combine_add(pdz, probe)
        five = combine_keep(pdz)
        sp = score(probe, y, LO, HI, mask, fold)
        pm, vm, _vs, cal = errcorr(probe, pdz, y, mask)

        # ---- CONTROLS, each of which must succeed and exercises the suspect clause -------
        # 1. The probe I just recomputed must reproduce the COMMITTED log. If the apparatus
        #    moved, every number below is incomparable to item 324's and means nothing.
        ref = json.load(open(RES + "logs/k106_probe.json"))["seeds"].get(str(s))
        c_probe = None
        if ref:
            c_probe = abs(sp[0] - ref["solo_probe"][0])
            print(f"      КОНТРОЛЬ зонд воспроизведён: {sp[0]:.4f} против "
                  f"{ref['solo_probe'][0]:.4f} в логе, расхождение {c_probe:.2e}", flush=True)
            if c_probe > 1e-6:
                print("ОТКАЗ: зонд не воспроизводится, прибор сдвинулся")
                return 2
        # 2. The gate function must SEPARATE the two known ends on these same data: the member
        #    band near 0.98 and the shipped probe near 0.90. A measurement that could only
        #    return one of them would be a query that cannot fail.
        band = float(np.mean([np.mean(v) for v in cal.values()]))
        pcorr = float(np.mean(vm))
        print(f"      КОНТРОЛЬ прибор различает: члены против своего среднего {band:.3f}, "
              f"отгруженный зонд {pcorr:.3f}, разнос {band - pcorr:.3f}", flush=True)
        if band - pcorr < 0.05:
            print("ОТКАЗ: прибор не различает известные концы, ворота нечитаемы")
            return 2

        res = out["seeds"].get(str(s), {})
        res["probe"] = {"solo": sp, "vs_five_mean": vm, "per_member": pm,
                        "calibration": cal, "seconds": round(t_probe, 1),
                        "control_reproduces_log": c_probe}
        # CONTROL on the six-member reference: the composition WITH the probe must score better
        # than without it, by about item 324's +0.0138 of rank. If it does not, `six` is not the
        # shipped composition and every correlation against it is against the wrong vector.
        s_six, s_five = score(six, y, LO, HI, mask, fold), score(five, y, LO, HI, mask, fold)
        gain = s_six[0] - s_five[0]
        print(f"      КОНТРОЛЬ состав: пять {s_five[0]:.4f} -> шесть {s_six[0]:.4f}, "
              f"+{gain:.4f} ранга (пункт 324 даёт +0.0138)", flush=True)
        if gain < 0.005:
            print("ОТКАЗ: добавление зонда не воспроизводит выигрыш пункта 324; "
                  "состав собран не тот")
            return 2
        res["composition"] = {"five": s_five, "six": s_six, "gain": gain}
        res.setdefault("heads", {})

        for kind in heads:
            P, dt = head_cached(E, y, mask, fold, kind, s, dg)
            sc = score(P, y, LO, HI, mask, fold)
            hpm, hvm, _hvs, _hcal = errcorr(P, pdz, y, mask)
            # The reference the probe did not need: the probe is already shipped, so a new head
            # is a SEVENTH member and must differ from the probe too.
            vs_probe = [float(np.corrcoef(P[e] - y[mask[:, e], e],
                                          probe[e] - y[mask[:, e], e])[0, 1])
                        for e in range(len(CYPS))]
            # THE binding reference. The probe shipped (items 325/327), so the composition a
            # seventh member competes with is `combine_add` -- the kept members PLUS the probe,
            # equally weighted, exactly as src/submit.py:1374 builds it. Correlating against the
            # five alone would answer a question nobody is asking any more.
            vs_six = [float(np.corrcoef(P[e] - y[mask[:, e], e],
                                        six[e] - y[mask[:, e], e])[0, 1])
                      for e in range(len(CYPS))]
            # Controls on the head itself: not a constant, not a copy of the probe.
            sd = [float(np.std(P[e])) for e in range(len(CYPS))]
            dmin = min(float(np.min(np.abs(P[e] - probe[e]))) for e in range(len(CYPS)))
            corr_five, corr_probe = float(np.mean(hvm)), float(np.mean(vs_probe))
            corr_six = float(np.mean(vs_six))
            binding = corr_six
            passes = (binding < GATE_CORR) and (sc[0] > WEAKEST - GATE_RANK_SLACK)
            # The answer the gate only approximates: admit the head as a SEVENTH member,
            # equally weighted beside the kept five and the probe, and score the composition.
            seven = [np.mean([Q[e] for k, Q in pdz if K._keep(e, k)] + [probe[e], P[e]], axis=0)
                     for e in range(len(CYPS))]
            s7 = score(seven, y, LO, HI, mask, fold)
            d7 = s7[0] - s_six[0]
            print(f"        седьмым членом: ранг {s_six[0]:.4f} -> {s7[0]:.4f} "
                  f"({d7:+.4f}), пара {s_six[1]:.4f} -> {s7[1]:.4f} ({s7[1] - s_six[1]:+.4f})"
                  f"  пол 0.0036", flush=True)
            res["heads"][kind] = {
                "seventh": s7, "seventh_gain_rank": d7,
                "seventh_gain_pair": s7[1] - s_six[1],
                "solo": sc, "vs_five_mean": hvm, "vs_probe": vs_probe, "vs_six": vs_six,
                "per_member": hpm, "mean_vs_five": corr_five, "mean_vs_probe": corr_probe,
                "mean_vs_six": corr_six, "binding": binding, "passes_gate": bool(passes),
                "seconds": round(dt, 1), "sd": sd, "min_abs_minus_probe": dmin,
            }
            print(f"      {kind:<9} ранг {sc[0]:.4f} пара {sc[1]:.4f} | corr пять "
                  f"{corr_five:.3f} зонд {corr_probe:.3f} СОСТАВ {corr_six:.3f} | "
                  f"{'ПРОХОДИТ' if passes else 'НЕ проходит'} | {dt:.0f} с", flush=True)
            if min(sd) < 1e-9 or dmin < 1e-12:
                print(f"        ОТКАЗ: {kind} постоянна или совпала с зондом "
                      f"(sd={min(sd):.2e}, min|diff|={dmin:.2e})")
                return 2

        out["seeds"][str(s)] = res
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
        print(f"      записано: {LOG}", flush=True)

    # ---- the verdict table ------------------------------------------------------------
    print("\n" + "=" * 92)
    print("ВОРОТА 0.93 --- связывает корреляция против ОТГРУЖЕННОГО СОСТАВА (пять + зонд)")
    print("=" * 92)
    print("  %-10s %8s %8s   %8s %8s %8s  %s"
          % ("голова", "ранг", "пара", "vs пять", "vs зонд", "vs СОСТАВ", "вердикт"))
    print("  %-10s %8.4f %8s   %8.3f %8s %8s  %s"
          % ("ЗОНД(лин)", PROBE_SOLO, "-", PROBE_CORR, "-", "-", "отгружен, +0.0138"))
    print("  %-10s %8.4f %8s   %8.3f %8s %8s  %s"
          % ("члены", WEAKEST, "-", BAND_LO, "-", "-", "полоса избыточности"))
    for s in seeds:
        for kind, h in out["seeds"].get(str(s), {}).get("heads", {}).items():
            print("  %-10s %8.4f %8.4f   %8.3f %8.3f %8.3f  %s"
                  % (f"{kind}/с{s}", h["solo"][0], h["solo"][1], h["mean_vs_five"],
                     h["mean_vs_probe"], h.get("mean_vs_six", h["binding"]),
                     "ПРОХОДИТ" if h["passes_gate"] else "НЕ проходит"))
    print("\n  Ворота: связывающая корреляция < %.2f И одиночный ранг > %.4f "
          "(слабейший член %.4f минус запас %.2f)." % (GATE_CORR, WEAKEST - GATE_RANK_SLACK,
                                                       WEAKEST, GATE_RANK_SLACK))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
