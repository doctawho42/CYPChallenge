"""Do the three sub-floor members add up, or are they one signal in three costumes?

THE SITUATION. Three separate measurements have now produced a member that is consistently positive
and consistently under the macro-rank floor of 0.0036:

    пункт 335   gbm на 600 замороженных столбцах            +0.0023   знак 4/4
    пункт 342   графовая сеть, скрининг + четыре сети       +0.0006   знак 4/4
    пункт 345   второй чекпойнт, chemeleon, зондом          +0.0033   знак 4/4

Each was reported as not shipping, correctly, because the project's rule is that an effect smaller
than one seed's worth of noise is not a result. But three of them point the same way, and the
question nobody has asked is whether they are INDEPENDENT. If they are, equal-weight addition sums
their contributions and the total may clear the floor; if they are three reads of the same signal,
adding all three costs dilution and buys one of them.

This is the cheapest question left in the project: every prediction involved is already on disk.
`results/preds/k114_heads_seed0.json` holds the gbm head, `results/preds/k119_aux_seed0_lam1_n4.json`
holds the graph net, and the chemeleon probe is recomputed here through `probe_oof` unchanged and
cached beside them.

WHAT DECIDES IT, and it is not the total. The pairwise error correlation AMONG the three candidates
is the quantity that says whether addition can work, and it is reported before the compositions so
the numbers are read in the right order. Item 335 established that the 0.93 screen does not predict
a member's gain, so these correlations are offered as a MECHANISM for whatever the compositions
show, never as a substitute for measuring them.

Every subset is measured rather than the full set alone: six, six plus each candidate, six plus each
pair, six plus all three. Equal weights throughout, because item 316's fitted weights and item 280's
Caruana selector are both measured shut on this library and a weighted combiner would reopen a
settled question inside an unsettled one.

    uv run python verify/k122_stack.py --seeds 0,1,2,3

Writes results/logs/k122_stack.json. Exit: 0 ran, 2 a control failed, 3 an input is missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import itertools
import json

import numpy as np
import pandas as pd

from cypsplit import butina_folds, fold_digest
from submit import CYPS

import k106_probe as K
from k106_probe import score, members, probe_oof, combine_add, load_emb
from k121_chemeleon import load_any, NEW

LOG = RES + "logs/k122_stack.json"
CHEM_CACHE = RES + "preds/k122_chemeleon_probe_seed%d.json"
FLOOR = 0.0036
GOLDEN = "2d93c19815e14261"
CAND = ("chemeleon", "gbm", "граф")


def chemeleon_probe(rows, Y, M, fold, seed, dg):
    path = CHEM_CACHE % seed
    if _pl.Path(path).exists():
        st = json.load(open(path))
        if st.get("digest") == dg:
            print("      chemeleon из кэша", flush=True)
            return [np.asarray(v, float) for v in st["P"]]
    E, _named = load_any(NEW, rows)   # chemeleon carries names; the flag is for nameless files
    P = [probe_oof(E, Y, M, fold, e)[0] for e in range(len(CYPS))]
    json.dump({"digest": dg, "P": [[float(x) for x in v] for v in P]}, open(path, "w"))
    return P


def from_cache(path, key, M, dg, what, optional=False):
    if not _pl.Path(path).exists():
        if optional:
            print("      нет %s --- кандидат пропущен (%s)"
                  % (path.split("/")[-1], what), flush=True)
            return None
        print("нет %s --- %s" % (path, what))
        raise SystemExit(3)
    st = json.load(open(path))
    if st.get("digest") != dg:
        print("ОТКАЗ: %s собран на разбиении %s, а не %s" % (path, st.get("digest"), dg))
        raise SystemExit(2)
    raw = st["heads"][key] if key else st["P"]
    # Two shapes live in these caches and they must be told apart by STRUCTURE, not by which
    # file they came from: k114 stores four per-enzyme vectors of the masked lengths, k119 a
    # full (n_rows, 4) matrix. Guessing wrong silently misaligns every row.
    want = [int(M[:, e].sum()) for e in range(len(CYPS))]
    if len(raw) == len(CYPS) and [len(v) for v in raw] == want:
        return [np.asarray(v, float) for v in raw]
    A = np.asarray(raw, dtype=float)
    if A.shape == (len(M), len(CYPS)):
        return [A[M[:, e], e] for e in range(len(CYPS))]
    print("ОТКАЗ: %s имеет форму, которую я не узнаю: %s (ожидалось %s по ферментам "
          "или (%d, 4))" % (path, getattr(A, "shape", type(raw)), want, len(M)))
    raise SystemExit(2)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

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

    out = json.load(open(LOG)) if _pl.Path(LOG).exists() else {}
    out.setdefault("meta", {}).update(floor=FLOOR, candidates=list(CAND))
    out.setdefault("seeds", {})

    for s in seeds:
        fold, ncl = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        print("\n  сид %d: дайджест %s" % (s, dg), flush=True)
        if s == 0 and dg != GOLDEN:
            print("ОТКАЗ: разбиение сдвинулось, сид 0 дал %s" % dg)
            return 2

        pdz = members(s, XFULL, Y, M, fold, LO, HI)
        probe = [probe_oof(load_emb(rows), Y, M, fold, e)[0] for e in range(len(CYPS))]
        six = combine_add(pdz, probe)
        s_six = score(six, Y, LO, HI, M, fold)
        five = [np.mean([Q[e] for k, Q in pdz if K._keep(e, k)], 0) for e in range(len(CYPS))]
        s_five = score(five, Y, LO, HI, M, fold)
        print("      КОНТРОЛЬ состав: пять %.4f -> шесть %.4f, +%.4f (пункт 324: +0.0138)"
              % (s_five[0], s_six[0], s_six[0] - s_five[0]), flush=True)
        if s_six[0] - s_five[0] < 0.005:
            print("ОТКАЗ: зонд не воспроизводит выигрыш пункта 324")
            return 2

        cand = {
            "chemeleon": chemeleon_probe(rows, Y, M, fold, s, dg),
            "gbm": from_cache(RES + "preds/k114_heads_seed%d.json" % s, "gbm", M, dg,
                              "прогоните verify/k114_gate93.py --heads gbm"),
            "граф": from_cache(RES + "preds/k119_aux_seed0_lam1_n4.json" if s == 0 else
                               RES + "preds/k119_aux_seed%d_lam1_n4.json" % s, None, M, dg,
                               "прогоните verify/k119_dmpnnaux.py --lam 1 --nets 4",
                               optional=True),
        }
        cand = {k: v for k, v in cand.items() if v is not None}
        names_here = tuple(c for c in CAND if c in cand)
        print("      кандидатов на этом сиде: %s" % (", ".join(names_here),), flush=True)

        # --- the mechanism, printed BEFORE the compositions -----------------------------
        print("      корреляция ОШИБОК между кандидатами (среднее по ферментам):", flush=True)
        pair = {}
        for x, y2 in itertools.combinations(names_here, 2):
            v = [float(np.corrcoef(cand[x][e] - Y[M[:, e], e],
                                   cand[y2][e] - Y[M[:, e], e])[0, 1])
                 for e in range(len(CYPS))]
            pair["%s~%s" % (x, y2)] = v
            print("        %-20s %s   среднее %.3f"
                  % ("%s ~ %s" % (x, y2), " ".join("%.3f" % t for t in v),
                     float(np.mean(v))), flush=True)
        vs_six = {}
        for x in names_here:
            v = [float(np.corrcoef(cand[x][e] - Y[M[:, e], e], six[e] - Y[M[:, e], e])[0, 1])
                 for e in range(len(CYPS))]
            vs_six[x] = v
        print("      и каждого против СОСТАВА из шести: "
              + "  ".join("%s %.3f" % (x, float(np.mean(vs_six[x]))) for x in names_here), flush=True)

        # --- every subset, equal weights ------------------------------------------------
        kept = [[Q[e] for k, Q in pdz if K._keep(e, k)] for e in range(len(CYPS))]
        res = {}
        print("      составы (равные веса):", flush=True)
        for r in range(0, len(names_here) + 1):
            for sub in itertools.combinations(names_here, r):
                P = [np.mean(kept[e] + [probe[e]] + [cand[x][e] for x in sub], 0)
                     for e in range(len(CYPS))]
                sc = score(P, Y, LO, HI, M, fold)
                name = "шесть" if not sub else "шесть+" + "+".join(sub)
                res[name] = {"score": sc, "gain_rank": sc[0] - s_six[0],
                             "gain_pair": sc[1] - s_six[1]}
                print("        %-34s ранг %.4f (%+.4f)  пара %.4f (%+.4f)"
                      % (name, sc[0], sc[0] - s_six[0], sc[1], sc[1] - s_six[1]), flush=True)

        out["seeds"][str(s)] = {"six": s_six, "five": s_five, "subsets": res,
                                "pair_corr": pair, "vs_six": vs_six}
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)

    # --- summary across seeds ------------------------------------------------------------
    names = sorted({n for s in seeds if str(s) in out["seeds"]
                    for n in out["seeds"][str(s)]["subsets"]},
                   key=lambda n: (n.count("+"), n))
    print("\n" + "=" * 94)
    print("СКЛАДЫВАЮТСЯ ЛИ ТРИ ПОДПОРОГОВЫХ ЧЛЕНА, сиды %s" % seeds)
    print("=" * 94)
    print("  %-34s %9s %8s %7s %10s" % ("состав", "Δранг", "sd", "знак", "Δпара"))
    for n in names:
        have = [s for s in seeds if str(s) in out["seeds"]
                and n in out["seeds"][str(s)]["subsets"]]
        g = [out["seeds"][str(s)]["subsets"][n]["gain_rank"] for s in have]
        p = [out["seeds"][str(s)]["subsets"][n]["gain_pair"] for s in have]
        if not g:
            continue
        sd = float(np.std(g, ddof=1)) if len(g) > 1 else float("nan")
        print("  %-34s %+9.4f %8.4f %5d/%d %+10.4f"
              % (n, float(np.mean(g)), sd, int(sum(1 for x in g if x > 0)), len(g),
                 float(np.mean(p))))
    print("\n  пол ранга %.4f, пол ST-RAE 0.007" % FLOOR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
