"""Item 305's own measurement: the two clash policies side by side, and what the difference costs.

Written BLIND, before either arm of k97_dock.py finished, so that the shape of the comparison
could not be chosen to suit the answer. That is the same discipline item 305 applies to the
amendment itself: the rule is fixed before the number exists.

WHY THIS EXISTS AT ALL. `verify/k97_dock.py` prints a verdict for whichever policy it was run
under, and nothing prints the thing item 305 actually cares about -- the DIFFERENCE between
`--clash keep` (item 298 literally, 202 positive affinities left in place) and `--clash zero`
(the amendment, positives clipped to 0.0). If the two arms agree, the clash scores never mattered
and the amendment was insurance. If the literal arm passes where the clipped one does not, the
gain rode on outliers rather than on pocket complementarity, and item 305 says in advance that
this is to be reported as such rather than as a pass.

INDEPENDENCE, and it is deliberate. Every number below is recomputed from the per-seed records in
`seeds`, never read from the `итог` section k97 wrote. A summary is a claim like any other, and
this file is the second opinion on it -- twice in one day a wrong status survived because nothing
recomputed it (item 304's review recorded a finding as confirmed on zero votes, then killed a true
one). The script prints k97's own stored verdict beside the recomputed one and flags a mismatch.

Reads results/logs/k97_dock.json and results/logs/k97_dock_zero.json. Writes nothing.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import RES

import json
import numpy as np

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
# Per-enzyme rank floors, item 165. A gain under its own floor is not a result.
FLOOR = {"CYP1A2": 0.0061, "CYP2C9": 0.0071, "CYP2D6": 0.0049, "CYP3A4": 0.0033}
# After items 282-285 the per-enzyme member is kept only here; the other two are a harness
# control that must read exactly 0.000000 (item 298, condition 4).
LIVE = ("CYP1A2", "CYP2D6")
ARM_A, ARM_B = "A целевая", "B неверная изоформа"
ARM_C = "C ненаправленная"
MIN_SEEDS = 3          # item 298 fixed "at least 3 of the 4 seeds"; see item 304 on 0.75*len(d)

POLICIES = [("keep", "logs/k97_dock.json", "буквальный пункт 298"),
            ("zero", "logs/k97_dock_zero.json", "поправка пункта 305")]


def deltas(blob, arm, enz):
    """Per-seed Δrank for one (arm, enzyme), recomputed from `seeds` rather than from `итог`."""
    out = []
    for _seed, rec in sorted(blob.get("seeds", {}).items()):
        if "эталон" not in rec or arm not in rec:
            continue
        if enz not in rec[arm] or enz not in rec["эталон"]:
            continue
        out.append(float(rec[arm][enz]) - float(rec["эталон"][enz]))
    return np.array(out, float)


def verdict(d, enz):
    """Item 298's conditions 1 and 2, with item 304's fix to the sign rule."""
    if len(d) < MIN_SEEDS:
        return {"n": len(d), "mean": float("nan"), "sd": float("nan"), "sign": 0,
                "passes": False, "why": f"сидов {len(d)} < {MIN_SEEDS}"}
    sign = int((d > 0).sum())
    ok = d.mean() > FLOOR[enz] and sign >= MIN_SEEDS
    return {"n": len(d), "mean": float(d.mean()),
            "sd": float(d.std(ddof=1)) if len(d) > 1 else 0.0,
            "sign": sign, "passes": bool(ok),
            "why": "" if ok else ("ниже порога" if d.mean() <= FLOOR[enz] else "знак не держится")}


def main():
    blobs = {}
    for tag, rel, _ in POLICIES:
        p = RES + rel
        if not _pl.Path(p).exists():
            print(f"нет {p} --- плечо '{tag}' ещё не досчитало либо не запускалось", flush=True)
            return
        blobs[tag] = json.load(open(p))

    print("=== условия сборки, оба плеча ===", flush=True)
    for tag, _, human in POLICIES:
        b = blobs[tag]
        cl, im = b.get("клинчи", {}), b.get("импутировано", {})
        print(f"  {tag:>4s} ({human}): клинчей {cl.get('значений', '?')}, "
              f"политика '{cl.get('политика', '?')}', импутировано "
              f"{im.get('строк', '?')} строк ({im.get('доля', float('nan')):.3%}), "
              f"сидов {len(b.get('seeds', {}))}", flush=True)

    # The volume confound is scored on untouched raw in both arms by construction, so a
    # difference here would mean the two runs did not see the same block.
    v = [blobs[t].get("объёмный конфаунд", {}) for t, _, _ in POLICIES]
    same = all(abs(v[0].get(c, float("nan")) - v[1].get(c, float("nan"))) < 1e-12 for c in CYPS)
    print(f"  объёмный конфаунд совпадает в обоих плечах: {same}"
          + ("" if same else "  <-- РАСХОЖДЕНИЕ, плечи видели разные данные"), flush=True)

    print("\n=== Δранг по плечам и политикам (пересчитано из посидовых записей) ===", flush=True)
    print(f"  {'политика':>8s} {'арм':>20s} {'фермент':>8s} {'Δср.':>9s} {'sd':>7s} "
          f"{'знак':>6s} {'пол':>7s} {'вердикт':>9s}", flush=True)
    summ = {}
    for tag, _, _ in POLICIES:
        summ[tag] = {}
        for arm in (ARM_A, ARM_B, ARM_C):
            summ[tag][arm] = {}
            for enz in CYPS:
                d = deltas(blobs[tag], arm, enz)
                if not len(d):
                    continue
                r = verdict(d, enz)
                summ[tag][arm][enz] = r
                mark = "ПРОХОДИТ" if r["passes"] else "нет"
                print(f"  {tag:>8s} {arm:>20s} {enz:>8s} {r['mean']:+9.4f} {r['sd']:7.4f} "
                      f"{r['sign']}/{r['n']:<4d} {FLOOR[enz]:7.4f} {mark:>9s}", flush=True)

    print("\n=== условие 4: контроль стенда (обязан быть 0.000000) ===", flush=True)
    bad4 = []
    for tag, _, _ in POLICIES:
        for enz in ("CYP2C9", "CYP3A4"):
            r = summ[tag].get(ARM_A, {}).get(enz)
            if r is None:
                continue
            clean = abs(r["mean"]) < 5e-7 and r["sd"] < 5e-7
            if not clean:
                bad4.append((tag, enz, r["mean"]))
            print(f"  {tag:>4s} {enz}: Δ {r['mean']:+.6f} sd {r['sd']:.6f} -> "
                  f"{'чисто' if clean else 'НАРУШЕНО'}", flush=True)
    if bad4:
        print("  условие 4 нарушено --- поферментный член не должен входить в эти клетки;"
              "\n  любой ненулевой Δ здесь означает дефект стенда, а не результат.", flush=True)

    print("\n=== условия пункта 298, по живым клеткам и по каждому плечу отдельно ===", flush=True)
    for enz in LIVE:
        for tag, _, human in POLICIES:
            A = summ[tag].get(ARM_A, {}).get(enz)
            B = summ[tag].get(ARM_B, {}).get(enz)
            if A is None or B is None:
                continue
            beats = A["mean"] > B["mean"]
            adopted = A["passes"] and beats and not bad4
            print(f"  {enz} / {tag:>4s}: A {A['mean']:+.4f} (знак {A['sign']}/{A['n']}) против "
                  f"B {B['mean']:+.4f} -> целевая {'бьёт' if beats else 'НЕ бьёт'} контроль; "
                  f"порог {'взят' if A['passes'] else 'не взят'}"
                  + (f" ({A['why']})" if A['why'] else "")
                  + f"  =>  {'ПРИНЯТО' if adopted else 'ОТКЛОНЕНО'}", flush=True)

    print("\n=== ИЗМЕРЕНИЕ ПУНКТА 305: чего стоили клинчи ===", flush=True)
    print(f"  {'фермент':>8s} {'арм':>20s} {'keep':>9s} {'zero':>9s} {'keep-zero':>10s} "
          f"{'пол':>7s} {'больше пола':>12s}", flush=True)
    for enz in LIVE:
        for arm in (ARM_A, ARM_B, ARM_C):
            k = summ["keep"].get(arm, {}).get(enz)
            z = summ["zero"].get(arm, {}).get(enz)
            if k is None or z is None:
                continue
            diff = k["mean"] - z["mean"]
            print(f"  {enz:>8s} {arm:>20s} {k['mean']:+9.4f} {z['mean']:+9.4f} {diff:+10.4f} "
                  f"{FLOOR[enz]:7.4f} {'ДА' if abs(diff) > FLOOR[enz] else 'нет':>12s}",
                  flush=True)

    print("\n  Как это читать, записано до появления чисел:", flush=True)
    print("    -- оба плеча отклонены и разница меньше пола: клинчи ничего не решали,", flush=True)
    print("       поправка была страховкой, и вывод пункта 298 стоит на обоих прочтениях;", flush=True)
    print("    -- keep принято, zero отклонено: выигрыш держался на выбросах, и пункт 305", flush=True)
    print("       заранее велит докладывать это как отказ, а не как проход;", flush=True)
    print("    -- zero принято, keep отклонено: клинчи маскировали сигнал, и тогда результат", flush=True)
    print("       принадлежит поправке, а не буквальной предрегистрации --- говорить надо так.", flush=True)

    print("\n=== сверка с тем, что k97 записал сам (сводка --- такое же утверждение) ===", flush=True)
    for tag, _, _ in POLICIES:
        stored = blobs[tag].get("итог")
        if not stored:
            print(f"  {tag:>4s}: раздела 'итог' нет", flush=True)
            continue
        bad = []
        for arm, cells in stored.items():
            for enz, rec in cells.items():
                mine = summ[tag].get(arm, {}).get(enz)
                if mine is None:
                    continue
                if abs(float(rec.get("mean", np.nan)) - mine["mean"]) > 1e-9:
                    bad.append(f"{arm}/{enz}")
        print(f"  {tag:>4s}: расхождений со сводкой {len(bad)}"
              + (f" -> {bad}" if bad else " --- пересчёт согласен"), flush=True)


if __name__ == "__main__":
    main()
