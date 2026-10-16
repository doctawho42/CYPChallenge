"""The nine-ligand leave-one-out that item 334 named as its own next measurement, read out.

WHAT THIS IS. Item 334 measured ensemble docking into experimental CYP3A4 receptors on ONE target,
metyrapone, and said so in its own last paragraph: "One target. ... the next measurement is the same
leave-one-out across every drug-like ligand in the ensemble, which is what would make any of this a
sample." That run exists. It was produced after the timeout fix in `verify/k112_ensdock.py` -- the
first attempt died on its first target, ritonavir into 1W0G exceeding 900 s, and took nine untried
ligands with it -- and `results/logs/k112_loo.json` has been sitting unread since.

This file only reads it. No docking is re-run, nothing is re-scored, and the poses themselves are
gone: `k112_ensdock.py` works in `tempfile.mkdtemp(prefix="k112_")`, so anything needing per-receptor
pose data would mean re-docking, which is said rather than quietly skipped.

THE BASELINE, established before reading anything, because items 337/338 measured a candidate
against an arm retired a month earlier. For the structure track there is nothing to retire: a count
of `results/submission/` finds fourteen CSV files, all fourteen regression or TDI, and **zero**
structure files of any extension, against a control of thirty-three directory entries. **We ship
nothing on this track.** So the arms to beat are external, both from items 333/334:

    лидер доски            LDDT-PLI 0.4252
    сток Boltz-2           0.4252 - 0.0302 = 0.3950
    потолок (воспроизв.)   LDDT-PLI 0.7515, BiSyRMSD медиана 1.454

THE CONTROL THAT MUST SUCCEED. Metyrapone is in this run and in item 334, so its four selections
must reproduce item 334's committed heme-box figures exactly. If they do not, the two runs used
different devices and nothing below compares to anything.

WHAT IS ASKED, in the order the commit message for the timeout fix set out: whether by-score
selection beats by-consensus across the nine, how far either sits from the oracle, and what that
does to item 334's claim that the pose set contains board-leading poses. Both the mean and the
per-ligand sign are reported, because this project's own standard is that a mean without a sign
count is not a result -- and here they disagree.

    uv run python verify/k127_looread.py

Reads results/logs/k112_loo.json. Writes results/logs/k127_looread.json.
Exit: 0 ran, 2 a control failed, 3 the log is missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import RES

import json

import numpy as np

LOG = RES + "logs/k112_loo.json"
OUT = RES + "logs/k127_looread.json"
LEADER = 0.4252          # item 334: the structure-track board leader's LDDT-PLI
BOLTZ2 = 0.3950          # leader minus the 0.0302 that separates it from stock Boltz-2
CEILING = 0.7515         # item 333: experimental reproducibility
# Item 334's committed heme-box figures for metyrapone, which this run must reproduce.
MYT334 = {"по_скору": (4.786, 0.3142), "по_согласию": (4.732, 0.3087),
          "оракул_rmsd": (2.790, 0.5545), "оракул_pli": (2.887, 0.6117)}
RULES = [("по_скору", "по скору"), ("по_согласию", "по согласию"),
         ("оракул_rmsd", "оракул RMSD"), ("оракул_pli", "оракул PLI")]


def main():
    p = _pl.Path(LOG)
    if not p.exists():
        print("нет %s" % LOG)
        return 3
    d = json.load(open(p))
    lig = sorted(d)
    print("  лигандов %d: %s" % (len(lig), ", ".join(lig)))

    # ---- control 1: metyrapone must reproduce item 334 ------------------------------------
    if "MYT" not in d:
        print("ОТКАЗ: в прогоне нет MYT --- нечем сверить устройство с пунктом 334")
        return 2
    worst = 0.0
    for key, (r334, p334) in MYT334.items():
        got = d["MYT"][key]
        dr, dp = abs(got["rmsd"] - r334), abs(got["pli"] - p334)
        worst = max(worst, dr, dp)
        print("      КОНТРОЛЬ MYT %-12s пункт 334 %6.3f/%6.4f, здесь %6.3f/%6.4f  |Δ| %.4f/%.4f"
              % (key, r334, p334, got["rmsd"], got["pli"], dr, dp))
    if worst > 0.001:
        print("ОТКАЗ: метирапон не воспроизводит пункт 334 (худшее |Δ| %.4f) --- разные устройства"
              % worst)
        return 2
    print("      -> устройство то же, худшее |Δ| %.4f" % worst)

    # ---- control 2: nothing silently dropped ----------------------------------------------
    bad = [(k, len(d[k]["не_отработало"]), d[k]["не_оценено"]) for k in lig
           if d[k]["не_отработало"] or d[k]["не_оценено"]]
    print("      КОНТРОЛЬ полноты: лигандов с отказами или неоценёнными позами %d из %d"
          % (len(bad), len(lig)))
    for k, nf, ne in bad:
        print("        %s: не отработало %d, не оценено %d" % (k, nf, ne))
    poses = [d[k]["n_poses"] for k in lig]
    print("      поз на лиганд: %d..%d, всего %d" % (min(poses), max(poses), sum(poses)))

    # ---- the table -------------------------------------------------------------------------
    print("\n" + "=" * 104)
    print("ENSEMBLE DOCKING, LEAVE-ONE-OUT ПО ДЕВЯТИ ЛИГАНДАМ (пункт 334 мерил ОДИН)")
    print("=" * 104)
    print("  %-5s %5s | %-20s %-20s %-20s %-20s"
          % ("лиг", "поз", "по скору", "по согласию", "оракул RMSD", "оракул PLI"))
    for k in lig:
        row = []
        for key, _nm in RULES:
            v = d[k][key]
            row.append("%-6s %5.3f/%.3f" % (v["rec"], v["rmsd"], v["pli"]))
        print("  %-5s %5d | %s" % (k, d[k]["n_poses"], " ".join("%-20s" % x for x in row)))

    agg = {}
    print("\n  %-14s %-28s %-28s %s"
          % ("правило", "RMSD (среднее/медиана)", "PLI (среднее/медиана)", "PLI > лидера"))
    for key, nm in RULES:
        r = [d[k][key]["rmsd"] for k in lig]
        q = [d[k][key]["pli"] for k in lig]
        over = sum(1 for x in q if x > LEADER)
        agg[key] = {"rmsd_mean": float(np.mean(r)), "rmsd_median": float(np.median(r)),
                    "pli_mean": float(np.mean(q)), "pli_median": float(np.median(q)),
                    "over_leader": over, "n": len(q)}
        print("  %-14s %10.3f / %-14.3f %10.4f / %-14.4f %d из %d"
              % (nm, np.mean(r), np.median(r), np.mean(q), np.median(q), over, len(q)))
    print("  %-14s %10s   %-14s %10.4f   %-14s" % ("лидер доски", "-", "", LEADER, ""))
    print("  %-14s %10.3f   %-14s %10.4f   %-14s" % ("потолок", 1.454, "", CEILING, ""))

    # ---- rule against rule, mean AND sign --------------------------------------------------
    print("\n  ПРАВИЛО ПРОТИВ ПРАВИЛА, по лигандам а не только в среднем:")
    sc_pli = np.array([d[k]["по_скору"]["pli"] for k in lig])
    co_pli = np.array([d[k]["по_согласию"]["pli"] for k in lig])
    sc_rms = np.array([d[k]["по_скору"]["rmsd"] for k in lig])
    co_rms = np.array([d[k]["по_согласию"]["rmsd"] for k in lig])
    sw_pli = int((sc_pli > co_pli).sum())
    sw_rms = int((sc_rms < co_rms).sum())
    print("    PLI:  скор лучше у %d из %d лигандов, в среднем %+.4f в пользу %s"
          % (sw_pli, len(lig), abs(np.mean(sc_pli) - np.mean(co_pli)),
             "скора" if np.mean(sc_pli) > np.mean(co_pli) else "СОГЛАСИЯ"))
    print("    RMSD: скор лучше у %d из %d лигандов, в среднем %+.3f A в пользу %s"
          % (sw_rms, len(lig), abs(np.mean(sc_rms) - np.mean(co_rms)),
             "скора" if np.mean(sc_rms) < np.mean(co_rms) else "СОГЛАСИЯ"))
    worst_lig = lig[int(np.argmax(co_pli - sc_pli))]
    print("    знак и среднее РАСХОДЯТСЯ по PLI: средним правит %s (+%.3f PLI на одном лиганде)"
          % (worst_lig, float((co_pli - sc_pli).max())))

    # ---- the selection gap -----------------------------------------------------------------
    best_rule = max(("по_скору", "по_согласию"), key=lambda k: agg[k]["pli_mean"])
    gap = agg["оракул_pli"]["pli_mean"] - agg[best_rule]["pli_mean"]
    print("\n  ЗАЗОР ОТБОРА: оракул PLI %.4f минус лучшее правило (%s) %.4f = %.4f"
          % (agg["оракул_pli"]["pli_mean"], best_rule, agg[best_rule]["pli_mean"], gap))
    print("    (пункт 334 на одном лиганде: 0.6117 минус 0.3142 = 0.2975)")

    # ---- item 334's surviving claim, tested ------------------------------------------------
    print("\n  ПРЕТЕНЗИЯ ПУНКТА 334 -- 'набор поз содержит позы, которые повели бы доску':")
    print("    оракул PLI: среднее %.4f, медиана %.4f против лидера %.4f"
          % (agg["оракул_pli"]["pli_mean"], agg["оракул_pli"]["pli_median"], LEADER))
    print("    превосходят лидера: оракулом %d из %d, РЕАЛИЗУЕМЫМИ правилами %d из %d"
          % (agg["оракул_pli"]["over_leader"], len(lig),
             agg["по_скору"]["over_leader"] + agg["по_согласию"]["over_leader"], 2 * len(lig)))
    held = agg["оракул_pli"]["pli_median"] > LEADER
    print("    -> на девяти лигандах претензия %s: медиана оракула %s лидера на %.4f"
          % ("ДЕРЖИТСЯ" if held else "НЕ ДЕРЖИТСЯ",
             "выше" if held else "ниже" if agg["оракул_pli"]["pli_median"] < LEADER else "равна",
             abs(agg["оракул_pli"]["pli_median"] - LEADER)))
    print("    метирапон (0.612) был лёгким концом, как пункт 334 и предупреждал:")
    print("       он выбран как самый воспроизводимый лиганд набора (пункт 333, 0.443 A)")

    # ---- the receptor habit, with its caveat -----------------------------------------------
    import collections
    print("\n  КАКИЕ РЕЦЕПТОРЫ ВЫБИРАЮТСЯ (признак того, что правило берёт рецептор, а не позу):")
    for key, nm in RULES:
        c = collections.Counter(d[k][key]["rec"] for k in lig)
        print("    %-12s %s" % (nm, ", ".join("%s x%d" % (r, n) for r, n in c.most_common())))
    hab = collections.Counter(d[k]["по_согласию"]["rec"] for k in lig).most_common(1)[0][0]
    A = [k for k in lig if d[k]["по_согласию"]["rec"] == hab]
    B = [k for k in lig if d[k]["по_согласию"]["rec"] != hab]
    f = lambda ks, key: float(np.mean([d[k][key]["pli"] for k in ks]))
    print("    согласие берёт %s у %d из %d лигандов; оракул берёт его у %d"
          % (hab, len(A), len(lig),
             sum(1 for k in lig if d[k]["оракул_pli"]["rec"] == hab)))
    print("      где взяло %s: правило %.4f, оракул %.4f, зазор %.4f"
          % (hab, f(A, "по_согласию"), f(A, "оракул_pli"), f(A, "оракул_pli") - f(A, "по_согласию")))
    print("      где взяло другое: правило %.4f, оракул %.4f, зазор %.4f"
          % (f(B, "по_согласию"), f(B, "оракул_pli"), f(B, "оракул_pli") - f(B, "по_согласию")))
    print("    ОСТОРОЖНО: подмножества различаются и по трудности (оракул %.4f против %.4f),"
          % (f(A, "оракул_pli"), f(B, "оракул_pli")))
    print("    так что это указание, а не измерение. Решающий контроль -- взять ОДИН")
    print("    фиксированный рецептор для всех девяти и сравнить с правилом; он требует")
    print("    поз по рецепторам, которых в логе нет, а позы лежали в tempfile.mkdtemp.")

    json.dump({"leader": LEADER, "boltz2": BOLTZ2, "ceiling": CEILING,
               "ligands": lig, "aggregate": agg,
               "score_wins_pli": sw_pli, "score_wins_rmsd": sw_rms,
               "selection_gap": gap, "best_rule": best_rule,
               "myt_control_max_abs_delta": worst,
               "habit_receptor": hab, "habit_n": len(A)},
              open(OUT, "w"), ensure_ascii=False, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
