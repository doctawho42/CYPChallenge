"""The probe's alpha grid is pinned at its CEILING on three enzymes of four, and item 324 said in
as many words that widening it is a different experiment. This is that experiment.

WHAT IS PINNED. `verify/k106_probe.py:ALPHAS` is `[100, 1000, 10000]`, and item 324 records the
reasoning: "Carried in, never searched. Small alphas overfit this embedding badly, and sklearn's
default 0.1/1/10 is NOT used. RidgeCV's own generalised cross-validation picks within this grid on
TRAINING rows only, so the choice is nested and no alpha is chosen against the outcome. Widening it
is a DIFFERENT experiment and must be labelled one."

Building the submission with the one reachable checkpoint (item 355) printed the thing that makes
this worth running: `alpha на полной выборке: CYP1A2 10000 (верх), CYP2C9 10000 (верх), CYP2D6 10000
(верх), CYP3A4 1000 (внутри)`. Three of four sit at the top of the grid. A regulariser at its
boundary is not a tuned parameter, it is a parameter the grid refused to tune --- and on 2048
standardised columns against 1285 to 2335 rows, the direction it is pushing is the one that should
have been allowed.

**This is not alpha chosen against the outcome.** RidgeCV's generalised cross-validation picks
within the grid using TRAINING rows of each fold only, exactly as before; what changes is which
values it may consider. The grid is extended upward only, by decades, and the arms are reported for
every grid so the comparison is a sweep rather than a pick.

WHAT IS MEASURED. The probe as a sixth member on the five `_keep` admits, item 324's device, scored
by `verify/k106_probe.py:score` after the per-fold affine pair, for each candidate grid. Four seeds.
The reference is the shipped grid, and the gain that matters is a member's gain over the five, so an
improvement here is an improvement to the one lever that has ever moved this project's board
ordering (item 344).

THE CONTROLS. The shipped grid must reproduce item 355's +0.0148 for this embedding, or the device
is not that device. The chosen alphas are printed per fold and enzyme for every grid, because an
optimum that moves to the new ceiling means the grid is still too narrow and the answer is another
decade rather than this one. Seed 0's Butina digest must be the golden `2d93c19815e14261`. And the
embedding's row order is checked by name with the rolled-row control, as item 345 requires.

    uv run python verify/k134_probealpha.py --emb data/emb_surrogate_adme_train.npz --seeds 0,1,2,3

Writes results/logs/k134_probealpha.json.
Exit: 0 ran, 2 a control failed, 3 an input is missing.
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
from submit import CYPS, _keep

import k106_probe as K
from k106_probe import score, probe_oof

LOG = RES + "logs/k134_probealpha.json"
GOLDEN = "2d93c19815e14261"
MEMBERS = ("поферментно", "пул", "GP", "гребневая", "ствол")
SHIPPED = [100.0, 1000.0, 10000.0]
GRIDS = {
    "отгружаемая 1e2..1e4": SHIPPED,
    "+1 декада  1e2..1e5": SHIPPED + [1e5],
    "+2 декады  1e2..1e6": SHIPPED + [1e5, 1e6],
    "+3 декады  1e2..1e7": SHIPPED + [1e5, 1e6, 1e7],
    "только верх 1e5..1e7": [1e5, 1e6, 1e7],
}
ITEM355 = 0.0148      # surrogate_adme as a sixth member, four seeds


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--emb", default=D + "emb_surrogate_adme_train.npz")
    ap.add_argument("--seeds", default="0,1,2,3")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    for f in ("rows.csv", "cyp-challenge-TRAIN_inhibition.csv"):
        if not _pl.Path(D + f).exists():
            print("нет data/%s --- make setup" % f, file=_sys.stderr)
            return 3
    if not _pl.Path(a.emb).exists():
        print("нет %s" % a.emb, file=_sys.stderr)
        return 3

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)
    K.S.LO, K.S.HI = LO, HI

    # the embedding, order-checked by name with the rolled control (item 345)
    z = np.load(a.emb, allow_pickle=True)
    E = z["train"].astype(np.float64)
    nm = np.asarray(z["molecule_name"], object)
    r = rows.Molecule_Name.to_numpy(object)
    ok, rolled = bool((nm == r).all()), bool((np.roll(nm, 1) == r).all())
    print("вложение %s %s, порядок совпал: %s; сдвинутый контроль (ожидается False): %s"
          % (_pl.Path(a.emb).name, E.shape, ok, rolled))
    if not ok or rolled:
        print("КОНТРОЛЬ ПРОВАЛЕН: порядок вложения", file=_sys.stderr)
        return 2

    out = {"emb": _pl.Path(a.emb).name, "grids": {k: v for k, v in GRIDS.items()},
           "seeds": seeds, "per_seed": {}}
    for s in seeds:
        fold, _n = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        if s == 0 and dg != GOLDEN:
            print("ОТКАЗ: сид 0 дал разбиение %s" % dg, file=_sys.stderr)
            return 2
        mem_p = _pl.Path(RES + f"preds/members_seed{s}.json")
        if not mem_p.exists():
            print("  сид %d: нет %s --- пропуск" % (s, mem_p.name))
            continue
        mem = json.load(open(mem_p))
        pdz = [(k, [np.asarray(mem[k][e], float) for e in range(len(CYPS))]) for k in MEMBERS]
        kept = [[Q[e] for k, Q in pdz if _keep(e, k)] for e in range(len(CYPS))]
        five = score([np.mean(kept[e], 0) for e in range(len(CYPS))], Y, LO, HI, M, fold)
        cell = {"digest": dg, "пять": five, "arms": {}}
        print("\n  сид %d (дайджест %s): пять = ранг %.4f, пара %.4f"
              % (s, dg, five[0], five[1]), flush=True)
        for name, grid in GRIDS.items():
            K.ALPHAS = np.asarray(grid, float)
            P, picked = [], []
            for e in range(len(CYPS)):
                pe, al, _nc = probe_oof(E, Y, M, fold, e)
                P.append(pe)
                picked.append([float(x) for x in al])
            six = score([np.mean(kept[e] + [P[e]], 0) for e in range(len(CYPS))],
                        Y, LO, HI, M, fold)
            solo = score(P, Y, LO, HI, M, fold)
            at_top = sum(1 for al in picked for x in al if x >= max(grid) - 1e-9)
            n_pick = sum(len(al) for al in picked)
            cell["arms"][name] = {"six": six, "solo": solo, "alphas": picked,
                                  "d_rank": six[0] - five[0], "d_pair": six[1] - five[1],
                                  "at_ceiling": [at_top, n_pick]}
            print("    %-22s шестым %+.4f ранга, %+.4f пары | одиночн. %.4f | на верху сетки %d/%d"
                  % (name, six[0] - five[0], six[1] - five[1], solo[0], at_top, n_pick),
                  flush=True)
        K.ALPHAS = np.asarray(SHIPPED, float)
        out["per_seed"][str(s)] = cell
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)

    have = [s for s in seeds if str(s) in out["per_seed"]]
    if not have:
        print("ни одного сида не посчитано", file=_sys.stderr)
        return 3
    print("\n" + "=" * 96)
    print("%-22s %9s %8s %7s %9s %9s %12s" % ("сетка alpha", "Δранг", "sd", "знак", "Δпара",
                                              "одиночн.", "на верху"))
    print("=" * 96)
    best = None
    for name in GRIDS:
        v = np.asarray([out["per_seed"][str(s)]["arms"][name]["d_rank"] for s in have], float)
        pr = np.asarray([out["per_seed"][str(s)]["arms"][name]["d_pair"] for s in have], float)
        so = np.mean([out["per_seed"][str(s)]["arms"][name]["solo"][0] for s in have])
        top = sum(out["per_seed"][str(s)]["arms"][name]["at_ceiling"][0] for s in have)
        tot = sum(out["per_seed"][str(s)]["arms"][name]["at_ceiling"][1] for s in have)
        sd = float(np.std(v, ddof=1)) if len(v) > 1 else float("nan")
        print("%-22s %+9.4f %8.4f %5d/%d %+9.4f %9.4f %9d/%d"
              % (name, v.mean(), sd, int((v > 0).sum()), len(v), pr.mean(), so, top, tot))
        out.setdefault("summary", {})[name] = {"d_rank": float(v.mean()), "sd": sd,
                                               "sign": [int((v > 0).sum()), len(v)],
                                               "d_pair": float(pr.mean()), "solo": float(so),
                                               "at_ceiling": [top, tot]}
        if best is None or v.mean() > out["summary"][best]["d_rank"]:
            best = name
    ship = out["summary"]["отгружаемая 1e2..1e4"]
    print("\nконтроль: отгружаемая сетка даёт %+.4f против %+.4f пункта 355 (|Δ| %.4f)"
          % (ship["d_rank"], ITEM355, abs(ship["d_rank"] - ITEM355)))
    if abs(ship["d_rank"] - ITEM355) > 0.003:
        print("КОНТРОЛЬ ПРОВАЛЕН: устройство не воспроизводит пункт 355", file=_sys.stderr)
        return 2
    gain = out["summary"][best]["d_rank"] - ship["d_rank"]
    print("лучшая сетка: %s, %+.4f ранга --- расширение стоит %+.4f (пол 0.0036)"
          % (best, out["summary"][best]["d_rank"], gain))
    if out["summary"][best]["at_ceiling"][0] > 0 and best != "отгружаемая 1e2..1e4":
        print("ВНИМАНИЕ: оптимум всё ещё задевает верх ЭТОЙ сетки --- нужна ещё декада")
    json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
    print("\nзапись: %s" % LOG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
