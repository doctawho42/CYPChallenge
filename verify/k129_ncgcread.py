"""The two NCGC panel arms that were measured, cached, and never reported --- read out and nested.

WHY THIS EXISTS. `src/ablncgc.py` wrote four caches into `results/preds/`. The journal reports two
of them: `oof_ncgc.json` at item 153 and `oof_ncgc_w.json` at item 157. It reports neither
`oof_ncgc_w123.json` --- the `--seeds 1,2,3` replication item 157 called "queued and what settles
it" --- nor `oof_ncgc_per.json`, the per-enzyme arm that item 187 recorded as having "died instantly
three separate times" on a parser defect, "twice inside a queue whose wrapper reported success".
The parser was repaired in that same item. The repaired run's output has been on disk since.

Read straight out of those caches, in the currency `src/ablncgc.py` already used:

    пул + NCGC w0.1, 4 сида      1A2 +0.0141 знак 4/4 (пол 0.0061)   MACRO +0.0008 (пол 0.0036)
    поферментно + NCGC w0.1      MACRO +0.0112 ранга, знак 2/2,  -0.0087 пары, знак 2/2
      его перемешанный контроль  MACRO -0.0022 ранга,            +0.0037 пары  (в обе стороны мимо)
      несут 2C9 +0.0243 и 3A4 +0.0123, а НЕ 1A2 (+0.0060)

which inverts item 157's two conclusions --- "whatever the panel is worth, it is worth it on CYP1A2
alone" and "+0.0048 macro is below the 0.007 floor ... a repair, not a result" --- both of which were
statements about the POOLED model. `docs/transfer_axis.md` has the full reading and why the
per-enzyme fit is where a foreign protocol can enter without being confused for an enzyme contrast
(items 157 and 188 both measured that mechanism already).

THE CURRENCY NEEDS NO CONVERSION, which is the first thing this file establishes rather than
assumes. `src/ablncgc.py:269-278` scores per enzyme on that enzyme's masked vector --- `spearmanr(y,
p)` for rank and `strae` after `fit_apply(p, lo, hi, fold, u)` for the pair, macro as the mean of
four --- and that is `verify/k106_probe.py:score()` exactly. So the numbers above are in the same
units as item 324's +0.0138 and item 346's +0.0047. This file recomputes them through `score()` and
requires the agreement, so the claim is checked and not argued.

WHAT IS ACTUALLY DECIDED HERE, because the numbers above are member-level and `docs/teammate_asks.md`
states the ceiling that eats member-level gains: a delta sinks over the ensemble by averaging
(~delta/N) and again by redundancy, 3 to 22x on the fused enzymes. Items 290 and 319 are the worked
example of a real channel dying on exactly that. So the deciding measurement is the panel arm nested
over the shipped composition --- the same test item 290 applied and item 324's probe passed --- and
it needs no retraining, because the arm's out-of-fold predictions are in the cache. Three nestings
are measured, not one:

    ADDED      шесть + плечо              does it add as a seventh member
    SUBSTITUTED  шесть, плечо вместо «поферментно»   does it improve the member it is a version of
    STANDALONE   ранг самого плеча        item 348's admission screen, threshold 0.60

Item 345 measured that pair of questions for its own candidate and the answers differed in sign, so
both are reported here too.

THE CONTROLS THAT MUST SUCCEED.
  1. ALIGNMENT. For every (seed, arm) in the caches, the stored `preds` must reproduce the stored
     `table` row to 1e-4 in all ten published cells. A cache whose predictions do not reproduce its
     own table is misaligned against `rows.csv` and nothing below it means anything.
  2. THE SPLIT. Seed 0's Butina fold digest must be the golden `2d93c19815e14261`, as every verifier
     in this project requires.
  3. THE COMPOSITION. The six-member reference must reproduce item 324's +0.0138 for the probe, which
     is `verify/k122_stack.py`'s own control and is what shows the device here is that device.
  4. THE SHUFFLE. The shuffled-panel arm must not out-score the real one, in either currency. It is
     the arm `src/ablncgc.py` already writes --- identical rows, identical target marginal, the
     molecule-to-measurement link broken --- and item 242 is the precedent for a block whose shuffle
     beat it.

    uv run python verify/k129_ncgcread.py --seeds 0,1
    uv run python verify/k129_ncgcread.py --seeds 0,1 --no-nest     (caches only, no refit)

Reads results/preds/oof_ncgc{,_w,_w123,_per}.json plus data/feats.npz, data/rows.csv and the
challenge inhibition table. Writes results/logs/k129_ncgcread.json.
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
from k106_probe import score, members, probe_oof, combine_add, load_emb

LOG = RES + "logs/k129_ncgcread.json"
GOLDEN = "2d93c19815e14261"
FLOOR_MACRO = 0.0036
FLOOR_PER = {"CYP1A2": 0.0061, "CYP2C9": 0.0071, "CYP2D6": 0.0049, "CYP3A4": 0.0033}
SCREEN = 0.60          # item 348's standalone-rank admission threshold
PROBE_GAIN = 0.005     # the floor k122 uses for "the probe reproduced item 324"

# (cache, base arm, panel arm, shuffled arm or None). The per-enzyme pair is the unread one;
# the pooled pair is here so the two are read on one device and the comparison is not across files.
FAMILIES = [
    ("поферментно", ["oof_ncgc_per.json"],
     "поферментно", "поферментно+NCGC w0.1", "поферментно+NCGC w0.1, перемешанный"),
    ("пул", ["oof_ncgc_w.json", "oof_ncgc_w123.json"],
     "база", "+NCGC активные w0.1", None),
]


def load_caches(names):
    """Merge the tables and prediction dicts of several caches; later files do not overwrite."""
    table, preds, seen = [], {}, []
    for n in names:
        p = _pl.Path(RES + "preds/" + n)
        if not p.exists():
            print("нет %s" % p)
            raise SystemExit(3)
        d = json.load(open(p))
        table += d["table"]
        for k, v in d["preds"].items():
            preds.setdefault(k, v)
        seen.append(n)
    return table, preds, seen


def row_of(table, seed, arm):
    for r in table:
        if r["seed"] == seed and r["рука"] == arm:
            return r
    return None


def vectors(preds, seed, arm, M):
    """The four masked out-of-fold vectors for one (seed, arm), length-checked per enzyme."""
    out = []
    for e, c in enumerate(CYPS):
        k = "%d|%s|%s" % (seed, arm, c)
        if k not in preds:
            return None
        v = np.asarray(preds[k], float)
        if len(v) != int(M[:, e].sum()):
            print("ОТКАЗ: %r имеет %d значений, а маска %s даёт %d"
                  % (k, len(v), c, int(M[:, e].sum())))
            raise SystemExit(2)
        out.append(v)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1")
    ap.add_argument("--probe-emb", default=None,
                    help="какое замороженное вложение читает зонд шестым членом. По умолчанию "
                         "k106_probe.EMB, то есть chemprop_medium. На свежем клоне его нет "
                         "(пункт 355): чекпойнт лежит в репозитории третьего участника. "
                         "surrogate_adme публикуют организаторы и k132 мерит замену шестым "
                         "членом в +0.0148 ранга против +0.0138 у medium, знак 4/4, так что "
                         "состав ниже --- «пять плюс зонд», а какой именно зонд, печатается.")
    ap.add_argument("--no-nest", action="store_true",
                    help="only re-read the caches; skip the composition, which refits members")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    for f in ("rows.csv", "feats.npz", "cyp-challenge-TRAIN_inhibition.csv"):
        if not _pl.Path(D + f).exists():
            print("нет data/%s --- соберите данные и признаки (make setup)" % f)
            return 3

    if a.probe_emb:
        if not _pl.Path(a.probe_emb).exists():
            print("нет %s" % a.probe_emb, file=_sys.stderr)
            return 3
        K.EMB = a.probe_emb
    if not a.no_nest and not _pl.Path(K.EMB).exists():
        print("нет %s --- вложение зонда. Передайте --probe-emb или --no-nest." % K.EMB,
              file=_sys.stderr)
        return 3
    print("зонд читает: %s" % _pl.Path(K.EMB).name)

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

    out = {"meta": {"floor_macro": FLOOR_MACRO, "floor_per": FLOOR_PER, "screen": SCREEN,
                    "seeds": seeds, "nested": not a.no_nest}, "families": {}}

    for fam, files, base_arm, panel_arm, shuf_arm in FAMILIES:
        table, preds, seen = load_caches(files)
        print("\n" + "=" * 92)
        print("СЕМЕЙСТВО %s   (%s)" % (fam, ", ".join(seen)))
        print("=" * 92)
        res = {"files": seen, "base": base_arm, "panel": panel_arm, "shuffle": shuf_arm,
               "seeds": {}}

        for s in seeds:
            if row_of(table, s, panel_arm) is None:
                print("  сид %d: этой руки в кэше нет --- пропускаю" % s)
                continue
            fold, _ncl = butina_folds(list(rows.SMILES), seed=s)
            dg = fold_digest(fold)
            if s == 0 and dg != GOLDEN:
                print("ОТКАЗ: разбиение сдвинулось, сид 0 дал %s" % dg)
                return 2
            print("\n  сид %d: дайджест %s" % (s, dg), flush=True)

            arms, S = {}, {}
            for tag, arm in (("base", base_arm), ("panel", panel_arm), ("shuf", shuf_arm)):
                if arm is None:
                    continue
                P = vectors(preds, s, arm, M)
                if P is None:
                    continue
                arms[tag] = P
                S[tag] = score(P, Y, LO, HI, M, fold)
                # ---- control 1: the cache's own table must come back out of its own preds ----
                r = row_of(table, s, arm)
                worst = 0.0
                for e, c in enumerate(CYPS):
                    worst = max(worst, abs(S[tag][2][e] - r["%s rho" % c]),
                                abs(S[tag][3][e] - r["%s пара" % c]))
                worst = max(worst, abs(S[tag][0] - r["MACRO rho"]),
                            abs(S[tag][1] - r["MACRO пара"]))
                if worst > 1e-4:
                    print("КОНТРОЛЬ 1 ПРОВАЛЕН: %r не воспроизводит свою таблицу, max|Δ| %.2e"
                          % (arm, worst))
                    return 2
                print("      контроль 1: %-38s таблица воспроизведена, max|Δ| %.1e"
                      % (arm, worst), flush=True)

            if "base" not in arms or "panel" not in arms:
                print("  сид %d: нет пары база/плечо --- пропускаю" % s)
                continue

            d_rank = S["panel"][0] - S["base"][0]
            d_pair = S["panel"][1] - S["base"][1]
            print("      плечо против базы:  ранг %+.4f (пол %.4f)   пара %+.4f (пол 0.007)"
                  % (d_rank, FLOOR_MACRO, d_pair))
            print("        по ферментам: " + "  ".join(
                "%s %+.4f" % (c[3:], S["panel"][2][e] - S["base"][2][e])
                for e, c in enumerate(CYPS)))
            cell = {"base": S["base"], "panel": S["panel"], "d_rank": d_rank, "d_pair": d_pair,
                    "d_rank_per": [S["panel"][2][e] - S["base"][2][e] for e in range(len(CYPS))],
                    "d_pair_per": [S["panel"][3][e] - S["base"][3][e] for e in range(len(CYPS))]}

            # ---- control 4: the shuffled panel must lose, in both currencies ----
            if "shuf" in arms:
                sr = S["shuf"][0] - S["base"][0]
                sp = S["shuf"][1] - S["base"][1]
                print("      контроль 4: перемешанная панель  ранг %+.4f  пара %+.4f" % (sr, sp))
                if sr >= d_rank or sp <= d_pair:
                    print("КОНТРОЛЬ 4 ПРОВАЛЕН: перемешанная панель не хуже настоящей")
                    return 2
                cell.update(shuf=S["shuf"], d_rank_shuf=sr, d_pair_shuf=sp)

            # ---- the screen, and the nesting ----
            print("      одиночный ранг плеча %.4f против порога %.2f пункта 348 --- %s"
                  % (S["panel"][0], SCREEN, "проходит" if S["panel"][0] >= SCREEN else "НЕ проходит"))
            cell["screen_pass"] = bool(S["panel"][0] >= SCREEN)

            if not a.no_nest:
                pdz = members(s, XFULL, Y, M, fold, LO, HI)
                probe = [probe_oof(load_emb(rows), Y, M, fold, e)[0] for e in range(len(CYPS))]
                six = combine_add(pdz, probe)
                s_six = score(six, Y, LO, HI, M, fold)
                five = [np.mean([Q[e] for k, Q in pdz if _keep(e, k)], 0) for e in range(len(CYPS))]
                s_five = score(five, Y, LO, HI, M, fold)
                print("      КОНТРОЛЬ 3 состав: пять %.4f -> шесть %.4f, %+.4f (пункт 324: +0.0138)"
                      % (s_five[0], s_six[0], s_six[0] - s_five[0]), flush=True)
                if s_six[0] - s_five[0] < PROBE_GAIN:
                    print("КОНТРОЛЬ 3 ПРОВАЛЕН: зонд не воспроизводит выигрыш пункта 324")
                    return 2

                kept = [[Q[e] for k, Q in pdz if _keep(e, k)] for e in range(len(CYPS))]
                added = [np.mean(kept[e] + [probe[e]] + [arms["panel"][e]], 0)
                         for e in range(len(CYPS))]
                s_add = score(added, Y, LO, HI, M, fold)
                subst = [np.mean([Q[e] for k, Q in pdz if _keep(e, k) and k != "поферментно"]
                                 + [probe[e]]
                                 + ([arms["panel"][e]] if _keep(e, "поферментно") else []), 0)
                         for e in range(len(CYPS))]
                s_sub = score(subst, Y, LO, HI, M, fold)
                print("      ВЛОЖЕНО: седьмым  ранг %+.4f  пара %+.4f"
                      % (s_add[0] - s_six[0], s_add[1] - s_six[1]))
                print("               заменой  ранг %+.4f  пара %+.4f"
                      % (s_sub[0] - s_six[0], s_sub[1] - s_six[1]))
                corr = [float(np.corrcoef(arms["panel"][e] - Y[M[:, e], e],
                                          six[e] - Y[M[:, e], e])[0, 1])
                        for e in range(len(CYPS))]
                print("               корреляция ошибок с составом %.3f (ворота пункта 289: 0.93)"
                      % float(np.mean(corr)))
                cell.update(six=s_six, five=s_five, added=s_add, substituted=s_sub,
                            gain_added=s_add[0] - s_six[0], gain_subst=s_sub[0] - s_six[0],
                            pair_added=s_add[1] - s_six[1], pair_subst=s_sub[1] - s_six[1],
                            err_corr_vs_six=corr)

            res["seeds"][str(s)] = cell

        out["families"][fam] = res
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)

    # ---- summary over seeds, in the house format: mean, sd, sign ----
    print("\n" + "=" * 92)
    print("СВОДКА. Пол макро-ранга %.4f, пол ST-RAE 0.007." % FLOOR_MACRO)
    print("=" * 92)
    for fam, res in out["families"].items():
        cells = list(res["seeds"].values())
        if not cells:
            continue
        def agg(key):
            v = [c[key] for c in cells if key in c]
            if not v:
                return None
            sd = float(np.std(v, ddof=1)) if len(v) > 1 else float("nan")
            return float(np.mean(v)), sd, sum(1 for x in v if x > 0), len(v)
        for label, key, better_up in (("плечо против базы, ранг", "d_rank", True),
                                      ("плечо против базы, пара", "d_pair", False),
                                      ("вложено седьмым, ранг", "gain_added", True),
                                      ("вложено заменой, ранг", "gain_subst", True)):
            g = agg(key)
            if g is None:
                continue
            m, sd, pos, n = g
            sign = pos if better_up else n - pos
            print("  %-12s %-26s %+8.4f  sd %7.4f  знак %d/%d"
                  % (fam, label, m, sd, sign, n))
        per = [c["d_rank_per"] for c in cells if "d_rank_per" in c]
        if per:
            mp = np.mean(np.asarray(per, float), 0)
            print("  %-12s %-26s %s" % (fam, "по ферментам, ранг", "  ".join(
                "%s %+.4f (пол %.4f)" % (c[3:], mp[e], FLOOR_PER[c])
                for e, c in enumerate(CYPS))))
    print("\nзапись: %s" % LOG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
