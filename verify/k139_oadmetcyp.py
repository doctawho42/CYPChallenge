"""The organisers' own CYP-trained model as a member: P5-P8 all fail, and the frozen embedding
lands BELOW a random-weight control of the same architecture.

WHAT IS MEASURED AND WHY. Item 350 named the remaining form of the transfer axis as a warm-started
multi-head fine-tune over challenge+ChEMBL+qHTS and deferred it on GPU cost. Half of it is already
trained and published: `openadmet/cyp1a2-cyp2d6-cyp3a4-cyp2c9-chemeleon-v1`, a multitask chemprop
D-MPNN warm-started from CheMeleon with four outputs on exactly our four enzymes, fitted to 8068
ChEMBL 37 pIC50 structures with no validation split. `src/embed.py` refused it once -- "a model
already trained on CYP pIC50 data would make the comparison circular" -- which is true of the
comparison that file makes and says nothing about whether the checkpoint is a member.

The source is pinned in a way no previous external arm managed: **`data/chemeleon_X_train.csv` is
byte-for-byte this model's training set** -- the same 8068 SMILES in the same row order with the
same label counts (4800 / 2583 / 2523 / 1361). So the same data has now been routed three ways:
as training ROWS through the trunk's third head (item 319, channel +0.0043 against a +0.0105 gate,
closed), as THEIR trained model (here), and as our own learner on their rows
(`verify/k140_extmodel.py`). The source is held fixed and only the role changes, which is exactly
the question `docs/transfer_axis.md` was built to answer.

THE ARMS, and the first one has no precedent: the model emits four direct estimates of our four
targets, so it enters the composition with ZERO parameters fitted on our labels. All twenty-seven
previous candidates were encoders with no head on our endpoint.

    арм                        одиночный ранг   шестым    седьмым   корр. с surrogate_adme
    прямые (0 параметров)             0.2881   -0.0138    -0.0080            0.778
    зонд(эмб 2048)                    0.4372   -0.0162    -0.0139            0.870
    зонд(голова 1024)                 0.3415   -0.0106    -0.0073            0.819
    зонд(4 предсказания)              0.3126   -0.0067    -0.0043            0.814
                                                            (сид 0; сводка по сидам ниже)

**Every pre-registered prediction fails.** P5 wanted item 348's standalone threshold of 0.60 and the
best arm reaches 0.4372. P6 wanted it above `surrogate_adme`'s 0.6287. P7 wanted a positive sixth
member and all four arms are negative. P8 wanted a positive seventh and an error correlation below
0.90, and the seventh slot is negative on every arm.

**And the number that makes this a finding rather than another closed door:** item 343's
same-device random-weight control is 0.4833 (`rand2048`), so the best arm here is **below random
weights**. A model fine-tuned on our own endpoint, by the organisers, for 100 epochs, is a worse
frozen feature extractor for our labels than untrained weights of the same shape. Item 343's
"freeze, don't fine-tune" now has its sharpest instance, and the comparison is run again on THIS
architecture rather than quoted (`--random 0` in `src/embed_cypmodel.py`, which inverts control 1:
random weights must FAIL to reproduce the model's training labels, and they do, at ρ 0.1353).

WHY, MEASURED RATHER THAN ARGUED, and it refutes the pre-registration's own reasoning. The
pre-registration argued that ChEMBL's defect is LEVEL (+0.22 to +0.60, items 60-63), that Spearman
is invariant to level exactly, and so the defect that closed the external-rows axis does not touch
the quantity item 353 says we are behind on. The source agrees with ours on ordering well enough to
support that: on the 64 shared structures ρ is 0.821 / 0.667 / 0.426 / 0.750 with offsets +0.215 to
+0.603, reproducing items 60-63's table. So the labels are not the problem. **The model is: it
reproduces its own 8068 training labels at ρ 0.9987 and orders our 4905 unselected library at ρ
0.2881.** A no-split model trained to convergence on a publication-selected sample is fitted to that
sample's distribution, and our library is not it. The defect is generalisation off the source's
distribution, not the metric and not the offset.

WHAT IT SAYS ABOUT FIRST PLACE, which is the reason this was worth a run. Item 341 records that the
best readable entrant trains on **challenge+ChEMBL+qHTS together**. It is now measured that ChEMBL
ALONE transfers at ρ 0.29 as a model and +0.0043 as rows, and item 356 measured that the qHTS panel
dies over the ensemble. So what works in that configuration is the UNION with our own distribution,
not either external source on its own -- which is precisely the joint multi-head fine-tune item 350
named and deferred, and the one thing on this axis that genuinely needs a GPU.

    uv run python verify/k139_oadmetcyp.py --seeds 0,1

Needs data/emb_oadmetcyp_train.npz from `src/embed_cypmodel.py`, run out of process in the chemprop
environment, and optionally data/emb_oadmetcyp_rand_train.npz for the random-weight arm. Reads the
member caches and item 349's probe caches. Writes results/logs/k139_oadmetcyp.json.

Controls: row order against rows.csv with a rolled negative; the five-member reference must
reproduce k132/k134/k135's 0.6434 at seed 0 to 1e-4; the fold digest must be the golden one; the
model's output columns are mapped to our enzymes BY NAME rather than by position; and the
random-weight arm is read if present, with item 343's rand2048 quoted when it is not.
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
import time

import numpy as np
import pandas as pd

from cypsplit import butina_folds, fold_digest
from submit import CYPS, _keep

import k106_probe as K
from k106_probe import score, probe_oof

LOG = RES + "logs/k139_oadmetcyp.json"
GOLDEN = "2d93c19815e14261"
MEMBERS = ("поферментно", "пул", "GP", "гребневая", "ствол")
EMB = D + "emb_oadmetcyp_train.npz"
RND = D + "emb_oadmetcyp_rand_train.npz"
FLOOR, FLOOR_PAIR = 0.0036, 0.007
SCREEN = 0.60            # item 348's standalone admission threshold
ITEM355 = 0.6287         # surrogate_adme's standalone rank (items 355/360), P6's bar
ITEM343_RAND = 0.4833    # rand2048, item 343's same-architecture random-weight control
ITEM361_CORR = 0.907     # the lowest error correlation on record (item 361)
# The published model's task order, which is NOT ours.
TASKS = ["OPENADMET_LOGAC50_CYP3A4", "OPENADMET_LOGAC50_CYP2D6",
         "OPENADMET_LOGAC50_CYP2C9", "OPENADMET_LOGAC50_CYP1A2"]


def taskmap(stored):
    """Column of the model's output for each of OUR enzymes, by name rather than by position."""
    names = [str(x) for x in stored]
    out = []
    for c in CYPS:
        hit = [j for j, n in enumerate(names) if n.upper().endswith(c.upper())]
        if len(hit) != 1:
            raise SystemExit("колонка модели для %s не опознана среди %r" % (c, names))
        out.append(hit[0])
    return out


def errcorr(a, b, Y, M):
    """Mean over enzymes of the correlation between the two members' residuals."""
    v = []
    for e in range(len(CYPS)):
        m = M[:, e]
        v.append(float(np.corrcoef(a[e] - Y[m, e], b[e] - Y[m, e])[0, 1]))
    return float(np.mean(v)), v


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    for f in ("rows.csv", "cyp-challenge-TRAIN_inhibition.csv"):
        if not _pl.Path(D + f).exists():
            print("нет data/%s --- make setup" % f, file=_sys.stderr)
            return 3
    if not _pl.Path(EMB).exists():
        print("нет %s --- сначала src/embed_cypmodel.py в /tmp/encenv" % EMB, file=_sys.stderr)
        return 3

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)
    K.S.LO, K.S.HI = LO, HI

    z = np.load(EMB, allow_pickle=True)
    PRED, E2048, F1024 = z["pred"], z["emb"], z["ffn"]
    col = taskmap(z["tasks"])
    print("вложение: pred %s, emb %s, ffn %s" % (PRED.shape, E2048.shape, F1024.shape))
    print("колонки модели -> наши ферменты: " +
          ", ".join("%s<-%d" % (c, col[e]) for e, c in enumerate(CYPS)))

    # ---- control 1: row order, by SMILES, with a deliberately shifted negative ----------
    smi = [str(x) for x in z["smiles"]]
    ok = smi == list(rows.SMILES)
    shifted = smi[1:] + smi[:1] == list(rows.SMILES)
    print("контроль 1: порядок строк совпал: %s; сдвинутый контроль (ожидается False): %s"
          % (ok, shifted))
    if not ok or shifted:
        print("КОНТРОЛЬ 1 ПРОВАЛЕН", file=_sys.stderr)
        return 2
    if len(E2048) != len(rows):
        print("КОНТРОЛЬ: %d строк против %d" % (len(E2048), len(rows)), file=_sys.stderr)
        return 2

    SETS = {"зонд(эмб 2048)": E2048, "зонд(голова 1024)": F1024,
            "зонд(4 предсказания)": PRED.astype(np.float64)}
    if _pl.Path(RND).exists():
        zr = np.load(RND, allow_pickle=True)
        if [str(x) for x in zr["smiles"]] != list(rows.SMILES):
            print("КОНТРОЛЬ: случайное вложение в другом порядке", file=_sys.stderr)
            return 2
        SETS["зонд(случайные веса)"] = zr["emb"]
        print("контроль 4: случайные веса той же архитектуры прочитаны из %s"
              % _pl.Path(RND).name)
    else:
        print("контроль 4: %s нет --- эталон берётся из пункта 343, rand2048 %.4f"
              % (_pl.Path(RND).name, ITEM343_RAND))

    per_seed, used = {}, []
    for s in seeds:
        mem_p = _pl.Path(RES + f"preds/members_seed{s}.json")
        if not mem_p.exists():
            print("  сид %d: нет %s --- пропуск" % (s, mem_p.name))
            continue
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        if s == 0 and dg != GOLDEN:
            print("ОТКАЗ: разбиение сдвинулось, сид 0 дал %s" % dg, file=_sys.stderr)
            return 2
        mem = json.load(open(mem_p))
        if set(mem) != set(MEMBERS):
            print("ОТКАЗ: кэш членов сида %d содержит %s" % (s, sorted(mem)), file=_sys.stderr)
            return 2
        pdz = [(k, [np.asarray(mem[k][e], float) for e in range(len(CYPS))]) for k in MEMBERS]
        kept = [[Q[e] for k, Q in pdz if _keep(e, k)] for e in range(len(CYPS))]
        five = [np.mean(kept[e], 0) for e in range(len(CYPS))]
        cell = {"digest": dg, "пять": score(five, Y, LO, HI, M, fold)}
        print("\n  сид %d: дайджест %s, пять = ранг %.4f, пара %.4f"
              % (s, dg, cell["пять"][0], cell["пять"][1]), flush=True)

        # the shipped sixth member we can actually build, as the seventh-slot baseline
        sixth = None
        p_sa = _pl.Path(RES + f"preds/k122_surrogate_adme_probe_seed{s}.json")
        if p_sa.exists():
            st = json.load(open(p_sa))
            if st.get("digest") == dg:
                sixth = [np.asarray(v, float) for v in st["P"]]
                cell["пять+sa"] = score([np.mean(kept[e] + [sixth[e]], 0)
                                         for e in range(len(CYPS))], Y, LO, HI, M, fold)
                print("    пять+surrogate_adme = ранг %.4f (%+.4f)"
                      % (cell["пять+sa"][0], cell["пять+sa"][0] - cell["пять"][0]), flush=True)

        # the model's own four outputs, as a member with ZERO parameters fitted on our labels
        direct = [PRED[M[:, e], col[e]].astype(np.float64) for e in range(len(CYPS))]
        arms = {"прямые (0 параметров)": direct}
        for tag, Xs in SETS.items():
            t0 = time.time()
            arms[tag] = [probe_oof(np.asarray(Xs, float), Y, M, fold, e)[0]
                         for e in range(len(CYPS))]
            print("    %-24s %4.0f с" % (tag, time.time() - t0), flush=True)

        for tag, P in arms.items():
            solo = score(P, Y, LO, HI, M, fold)
            six = score([np.mean(kept[e] + [P[e]], 0) for e in range(len(CYPS))],
                        Y, LO, HI, M, fold)
            rec = {"solo": solo, "sixth": six}
            if sixth is not None:
                rec["seventh"] = score(
                    [np.mean(kept[e] + [sixth[e], P[e]], 0) for e in range(len(CYPS))],
                    Y, LO, HI, M, fold)
                rec["corr_sa"] = errcorr(P, sixth, Y, M)
                rec["corr_five"] = errcorr(P, five, Y, M)
            cell[tag] = rec
            print("    %-24s одиночн. %.4f | шестым %+.4f | седьмым %+.4f | корр. с sa %s"
                  % (tag, solo[0], six[0] - cell["пять"][0],
                     (rec["seventh"][0] - cell.get("пять+sa", cell["пять"])[0])
                     if sixth is not None else float("nan"),
                     ("%.3f" % rec["corr_sa"][0]) if sixth is not None else "—"), flush=True)
        per_seed[str(s)] = cell
        used.append(s)

    if not used:
        print("нет ни одного кэша членов", file=_sys.stderr)
        return 3

    # ---- control 2: the five-member reference must reproduce k132/k134/k135 -------------
    if 0 in used and abs(per_seed["0"]["пять"][0] - 0.6434) > 1e-4:
        print("КОНТРОЛЬ 2 ПРОВАЛЕН: пять на сиде 0 дали %.4f, а k132/k134/k135 дали 0.6434"
              % per_seed["0"]["пять"][0], file=_sys.stderr)
        return 2
    if 0 in used:
        print("\nконтроль 2: пять членов на сиде 0 дают ранг %.4f (k132/k134/k135: 0.6434)"
              % per_seed["0"]["пять"][0])

    tags = [t for t in per_seed[str(used[0])] if t not in ("digest", "пять", "пять+sa")]
    print("\n%s" % ("=" * 104))
    print("%-24s %9s %9s | %9s %7s | %9s %7s | %7s"
          % ("арм", "одиночн.", "Δпара", "Δ шестым", "знак", "Δ седьмым", "знак", "корр."))
    print("=" * 104)
    out = {"seeds": used, "per_seed": per_seed, "summary": {}}
    for t in tags:
        sol = np.array([per_seed[str(s)][t]["solo"][0] for s in used])
        solp = np.array([per_seed[str(s)][t]["solo"][1] for s in used])
        d6 = np.array([per_seed[str(s)][t]["sixth"][0] - per_seed[str(s)]["пять"][0]
                       for s in used])
        has7 = all("seventh" in per_seed[str(s)][t] for s in used)
        d7 = (np.array([per_seed[str(s)][t]["seventh"][0] - per_seed[str(s)]["пять+sa"][0]
                        for s in used]) if has7 else np.array([np.nan]))
        cr = (float(np.mean([per_seed[str(s)][t]["corr_sa"][0] for s in used]))
              if has7 else float("nan"))
        out["summary"][t] = {"solo_rank": float(sol.mean()), "solo_pair": float(solp.mean()),
                             "d6": float(d6.mean()), "sign6": [int((d6 > 0).sum()), len(d6)],
                             "d7": float(np.nanmean(d7)),
                             "sign7": [int(np.nansum(d7 > 0)), len(used)] if has7 else None,
                             "corr_sa": cr}
        print("%-24s %9.4f %9.4f | %+9.4f %5d/%d | %+9.4f %5s | %7s"
              % (t, sol.mean(), solp.mean(), d6.mean(), int((d6 > 0).sum()), len(d6),
                 np.nanmean(d7), ("%d/%d" % (int(np.nansum(d7 > 0)), len(used))) if has7 else "—",
                 ("%.3f" % cr) if has7 else "—"))
    print("\nполы: макро ранг %.4f, ST-RAE %.3f; экран пункта 348 --- одиночный ранг %.2f"
          % (FLOOR, FLOOR_PAIR, SCREEN))
    print("эталоны: surrogate_adme одиночный %.4f (P6), случайные веса пункта 343 %.4f,"
          % (ITEM355, ITEM343_RAND))
    print("         самая низкая корреляция ошибок на записи %.3f (пункт 361)" % ITEM361_CORR)

    best = max(out["summary"], key=lambda k: out["summary"][k]["solo_rank"])
    b = out["summary"][best]
    print("\nпредсказания, вынесенные ДО запуска (docs/transfer_axis.md):")
    print("  P5 одиночный ранг >= %.2f:      %s  (%.4f, арм %s)"
          % (SCREEN, "ДЕРЖИТСЯ" if b["solo_rank"] >= SCREEN else "ПРОВАЛИЛОСЬ",
             b["solo_rank"], best))
    print("  P6 одиночный > %.4f:           %s  (%+.4f)"
          % (ITEM355, "ДЕРЖИТСЯ" if b["solo_rank"] > ITEM355 else "ПРОВАЛИЛОСЬ",
             b["solo_rank"] - ITEM355))
    s6 = b["sign6"]
    print("  P7 шестым > пола при 2/2:       %s  (%+.4f, знак %d/%d)"
          % ("ДЕРЖИТСЯ" if (b["d6"] > FLOOR and s6[0] == s6[1]) else "ПРОВАЛИЛОСЬ",
             b["d6"], s6[0], s6[1]))
    if b["sign7"] is not None:
        s7 = b["sign7"]
        p8 = b["d7"] > FLOOR and s7[0] == s7[1] and b["corr_sa"] < 0.90
        print("  P8 седьмым > пола И корр < 0.90: %s  (%+.4f, знак %d/%d, корр %.3f)"
              % ("ДЕРЖИТСЯ" if p8 else "ПРОВАЛИЛОСЬ", b["d7"], s7[0], s7[1], b["corr_sa"]))
    out["best_arm"] = best
    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
    print("\nзапись: %s" % LOG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
