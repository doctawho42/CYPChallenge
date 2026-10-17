"""The shipped sixth member cannot be rebuilt here, and the available replacement is measured
BETTER. Reads only: no training, no embedding, no network.

THE SITUATION. Items 324/325/327 ship a RidgeCV head on frozen `chemprop_medium` as a sixth
member, +0.0138 of macro rank at sign 8/8, and item 344 records it as the only intervention in
nine uploads that moved the board's ORDERING. That checkpoint is a file in a third entrant's
GitHub repository (`jeremycheminf/openadmet_scripts`, the `adme_pretrain` D-MPNNs, item 324), and
`data/emb_chemprop_medium*.npz` is not committed -- item 341 records that the checkpoint itself was
never on disk either. So on a fresh checkout the shipped composition cannot be reproduced at all,
which is a provenance defect of item 318's family and is what this file ran into rather than
discovered.

What IS reachable: `surrogate_adme__s42.pt` is one of the fifteen checkpoints the organisers publish
(`openadmet/chemprop-foundation-pretraining-weights`), and item 349 measured its standalone rank at
0.6263 against medium's 0.6261 -- the best single candidate on record, and the item that corrected
the criterion to "the pretraining target being an ENDPOINT".

WHAT IS MEASURED. Whether swapping the sixth member costs anything. The device is item 324's: the
five members `_keep` admits, plus one frozen-encoder probe at equal weight, scored by
`verify/k106_probe.py:score` after the per-fold affine pair. Nothing is trained: the five members
come from `results/preds/members_seed{0..3}.json` and the probes from the `k122_*_probe_seed*.json`
caches item 349 wrote, so this is the route item 117 used to re-read item 68 and item 290 used to
close the external head.

    пять                      ранг 0.6421   пара 0.6382
    пять + surrogate_adme          0.6569        0.6215    +0.0148 / -0.0168   знак 4/4
    пять + chemeleon               0.6555        0.6202    +0.0134 / -0.0181   знак 4/4
    пять + оба                     0.6593        0.6162    +0.0172 / -0.0220   знак 4/4
    пункт 324, medium шестым                               +0.0138 / -0.0197   знак 8/8 (8 сидов)

So the substitution costs nothing measurable, and BOTH together beat either alone -- but
`chemprop_chemeleon` sits in the same unreachable repository as medium, so only the
`surrogate_adme` arm can be built into a submission from here. That is why the number matters: it
is the licence for `src/submit.py --probe-emb`.

THE CONTROLS. Every probe cache carries the fold digest it was built on and must match the seed's
own digest, so a cache from a different split cannot be silently averaged in; seed 0's must be the
golden `2d93c19815e14261`. The five-member reference must reproduce across seeds, and the
comparison is paired within seed, which is what makes a 4/4 sign count mean anything on four
points. And the per-enzyme deltas are printed, because `_keep` admits all five members on CYP2D6
and only two or three on the others, so a sixth member does not enter the four enzymes equally.

    uv run python verify/k132_sixthsub.py

Reads results/preds/members_seed*.json, results/preds/k122_*_probe_seed*.json, data/rows.csv and
the inhibition table. Writes results/logs/k132_sixthsub.json.
Exit: 0 ran, 2 a control failed, 3 an input is missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import itertools
import json

import numpy as np
import pandas as pd

from cypsplit import butina_folds, fold_digest
from submit import CYPS, _keep

import k106_probe as K
from k106_probe import score

LOG = RES + "logs/k132_sixthsub.json"
GOLDEN = "2d93c19815e14261"
MEMBERS = ("поферментно", "пул", "GP", "гребневая", "ствол")
CAND = ("surrogate_adme", "chemeleon")
FLOOR, FLOOR_PAIR = 0.0036, 0.007
ITEM324 = (0.0138, -0.0197)      # medium as a sixth, eight seeds


def main():
    for f in ("rows.csv", "cyp-challenge-TRAIN_inhibition.csv"):
        if not _pl.Path(D + f).exists():
            print("нет data/%s --- make setup" % f, file=_sys.stderr)
            return 3

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)
    K.S.LO, K.S.HI = LO, HI

    per_seed, seeds_used = {}, []
    for s in (0, 1, 2, 3):
        mem_p = _pl.Path(RES + f"preds/members_seed{s}.json")
        if not mem_p.exists():
            print("  сид %d: нет %s --- пропуск" % (s, mem_p.name))
            continue
        fold, _ncl = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        if s == 0 and dg != GOLDEN:
            print("ОТКАЗ: разбиение сдвинулось, сид 0 дал %s" % dg, file=_sys.stderr)
            return 2
        mem = json.load(open(mem_p))
        if set(mem) != set(MEMBERS):
            print("ОТКАЗ: кэш членов сида %d содержит %s" % (s, sorted(mem)), file=_sys.stderr)
            return 2
        pdz = [(k, [np.asarray(mem[k][e], float) for e in range(len(CYPS))]) for k in MEMBERS]

        cand = {}
        for tag in CAND:
            p = _pl.Path(RES + f"preds/k122_{tag}_probe_seed{s}.json")
            if not p.exists():
                continue
            st = json.load(open(p))
            if st.get("digest") != dg:
                print("ОТКАЗ: %s собран на разбиении %s, а не %s"
                      % (p.name, st.get("digest"), dg), file=_sys.stderr)
                return 2
            cand[tag] = [np.asarray(v, float) for v in st["P"]]

        kept = [[Q[e] for k, Q in pdz if _keep(e, k)] for e in range(len(CYPS))]
        cell = {"digest": dg, "n_kept": [len(kept[e]) for e in range(len(CYPS))]}
        cell["пять"] = score([np.mean(kept[e], 0) for e in range(len(CYPS))], Y, LO, HI, M, fold)
        names = [t for t in CAND if t in cand]
        for r in range(1, len(names) + 1):
            for sub in itertools.combinations(names, r):
                P = [np.mean(kept[e] + [cand[t][e] for t in sub], 0) for e in range(len(CYPS))]
                cell["пять+" + "+".join(sub)] = score(P, Y, LO, HI, M, fold)
        per_seed[str(s)] = cell
        seeds_used.append(s)
        print("  сид %d: дайджест %s, кандидатов %d, членов по ферментам %s"
              % (s, dg, len(names), cell["n_kept"]), flush=True)

    if not seeds_used:
        print("нет ни одного кэша членов --- нечего читать", file=_sys.stderr)
        return 3

    arms = sorted({k for c in per_seed.values() for k in c if k.startswith("пять")},
                  key=lambda k: (k.count("+"), k))
    print("\n%-26s %8s %8s | %9s %9s %7s" % ("состав", "ранг", "пара", "Δранг", "Δпара", "знак"))
    out = {"seeds": seeds_used, "per_seed": per_seed, "summary": {}}
    base = np.asarray([per_seed[str(s)]["пять"][:2] for s in seeds_used], float)
    for a in arms:
        have = [s for s in seeds_used if a in per_seed[str(s)]]
        v = np.asarray([per_seed[str(s)][a][:2] for s in have], float)
        b = np.asarray([per_seed[str(s)]["пять"][:2] for s in have], float)
        d = v - b
        sd = float(np.std(d[:, 0], ddof=1)) if len(d) > 1 else float("nan")
        out["summary"][a] = {"rank": float(v[:, 0].mean()), "pair": float(v[:, 1].mean()),
                             "d_rank": float(d[:, 0].mean()), "d_pair": float(d[:, 1].mean()),
                             "sd_rank": sd, "sign": [int((d[:, 0] > 0).sum()), len(d)]}
        print("%-26s %8.4f %8.4f | %+9.4f %+9.4f %5d/%d"
              % (a, v[:, 0].mean(), v[:, 1].mean(), d[:, 0].mean(), d[:, 1].mean(),
                 int((d[:, 0] > 0).sum()), len(d)))
    print("%-26s %8s %8s | %+9.4f %+9.4f %5s    пункт 324, medium, 8 сидов"
          % ("medium шестым", "", "", ITEM324[0], ITEM324[1], "8/8"))
    print("\nполы: макро ранг %.4f, ST-RAE %.3f" % (FLOOR, FLOOR_PAIR))

    # per-enzyme, because _keep does not admit a sixth member equally across the four
    print("\nпо ферментам, Δранг против пяти:")
    for a in arms:
        have = [s for s in seeds_used if a in per_seed[str(s)]]
        d = np.mean([[per_seed[str(s)][a][2][e] - per_seed[str(s)]["пять"][2][e]
                      for e in range(len(CYPS))] for s in have], 0)
        out["summary"][a]["d_rank_per"] = [float(x) for x in d]
        print("  %-26s " % a + "  ".join("%s %+.4f" % (c[3:], d[e])
                                         for e, c in enumerate(CYPS)))

    best = max(out["summary"], key=lambda k: out["summary"][k]["d_rank"])
    sa = out["summary"].get("пять+surrogate_adme")
    if sa is not None and sa["d_rank"] < FLOOR:
        print("\nВНИМАНИЕ: surrogate_adme шестым ниже пола --- замена НЕ оправдана")
    print("\nлучший читаемый состав: %s (%+.4f ранга)" % (best, out["summary"][best]["d_rank"]))
    print("из них собираемый в подачу: только surrogate_adme --- chemeleon и medium лежат")
    print("в репозитории третьего участника, а surrogate_adme публикуют организаторы.")
    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
    print("\nзапись: %s" % LOG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
