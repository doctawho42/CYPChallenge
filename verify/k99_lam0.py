"""The identity control items 174 and 175 rest on, supplied for the two-site trunk arms.

Written BLIND, before either control run finished, so the pass rule could not be chosen to fit
the answer -- the same discipline as verify/k98_clash_arms.py and item 305's amendment.

WHY IT WAS MISSING. Item 307 closed the parametric leaf partly on `results/preds/trunk_caltwo.json`
and `trunk_caltwo3a4.json`, two runs committed on 2 September and never written up. Both carry
lambda 0.3 and 3.0 only. Items 174 and 175 both rest on the lambda = 0 arm showing that the mode
flag does nothing without a screening term, so as committed those runs could not have been adopted
even had they been positive. This supplies the missing arm.

WHAT THE CONTROL ASSERTS, and why it is exact rather than statistical. a comment in `src/trunk.py`'s `main` says
it outright: "lambda = 0 does not touch g_of_pi at all, so that arm is unaffected by construction
and serves as the leak check." At lambda = 0 the loss has no screening term, `g_of_pi` is never
reached, and the instrument constants cannot enter. So a two-site mode at lambda = 0 must produce
predictions IDENTICAL to any other mode at lambda = 0 -- not close, identical. `calshift` already
carries lambda = 0 for seeds 0-3 on the same device (mps) and the same torch (2.13.0), so the
reference is on disk and no new reference run is needed.

PASS RULE, fixed here before the numbers exist:
  1. every (seed, compound, enzyme) prediction of caltwo|s|0.0 equals calshift|s|0.0 exactly,
     NaN mask included -- max |difference| must be 0.0, not merely small;
  2. the same for caltwo3a4|s|0.0;
  3. all four seeds present in both new files.
A non-zero difference means the mode flag leaks outside the screening channel, which would make
every mode comparison in items 174-176, 307 and this file's own tables suspect. A pass means the
two-site arms of item 307 are now controlled exactly as items 174 and 175 require, and item 307's
conclusion stands on a complete record.

Reads results/preds/trunk_{calshift,caltwo_lam0,caltwo3a4_lam0}.json. Writes nothing.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import RES

import json
import numpy as np

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
SEEDS = [0, 1, 2, 3]
REF = "calshift"
ARMS = [("caltwo", "preds/trunk_caltwo_lam0.json"),
        ("caltwo3a4", "preds/trunk_caltwo3a4_lam0.json")]


def load(rel):
    p = RES + rel
    if not _pl.Path(p).exists():
        return None, p
    return json.load(open(p)), p


def main():
    ref, refp = load("preds/trunk_calshift.json")
    if ref is None:
        print(f"нет эталона {refp}", flush=True)
        return
    print(f"эталон: {refp}  (устройство {ref['meta'].get('device')}, "
          f"torch {ref['meta'].get('torch')})", flush=True)

    blobs = {}
    for name, rel in ARMS:
        b, p = load(rel)
        if b is None:
            print(f"нет {p} --- контрольный прогон '{name}' ещё не завершён", flush=True)
            return
        blobs[name] = b
        print(f"  {name}: {p}, устройство {b['meta'].get('device')}, "
              f"torch {b['meta'].get('torch')}, сиды {b['meta'].get('seeds')}, "
              f"лямбды {b['meta'].get('lams')}", flush=True)

    # Same device and torch, or a bit-for-bit comparison is testing the hardware, not the flag.
    bad_env = [n for n, b in blobs.items()
               if b['meta'].get('device') != ref['meta'].get('device')
               or b['meta'].get('torch') != ref['meta'].get('torch')]
    if bad_env:
        print(f"\n  ВНИМАНИЕ: {bad_env} посчитаны на другом устройстве или другом torch, "
              f"чем эталон --- точное совпадение проверять нельзя", flush=True)

    print(f"\n  {'арм':>11s} {'сид':>4s} {'ячеек':>7s} {'макс |разность|':>16s} "
          f"{'маски совпали':>14s} {'вердикт':>9s}", flush=True)
    fails = []
    for name, _ in ARMS:
        for s in SEEDS:
            k_new, k_ref = f"{name}|{s}|0.0", f"{REF}|{s}|0.0"
            if k_new not in blobs[name]["preds"]:
                print(f"  {name:>11s} {s:>4d}   нет ключа {k_new}", flush=True)
                fails.append(f"{name} сид {s}: нет ключа"); continue
            if k_ref not in ref["preds"]:
                print(f"  {name:>11s} {s:>4d}   нет эталонного ключа {k_ref}", flush=True)
                fails.append(f"{name} сид {s}: нет эталона"); continue
            A = np.asarray(blobs[name]["preds"][k_new], float)
            B = np.asarray(ref["preds"][k_ref], float)
            if A.shape != B.shape:
                print(f"  {name:>11s} {s:>4d}   форма {A.shape} против {B.shape}", flush=True)
                fails.append(f"{name} сид {s}: форма"); continue
            same_mask = bool(np.array_equal(np.isnan(A), np.isnan(B)))
            m = ~np.isnan(A) & ~np.isnan(B)
            d = float(np.max(np.abs(A[m] - B[m]))) if m.any() else float("nan")
            ok = same_mask and d == 0.0
            print(f"  {name:>11s} {s:>4d} {int(m.sum()):>7d} {d:16.3e} "
                  f"{str(same_mask):>14s} {'СОВПАЛО' if ok else 'НЕТ':>9s}", flush=True)
            if not ok:
                fails.append(f"{name} сид {s}: max|d|={d:.3e}, маски {same_mask}")

    print("\n  === вердикт, по правилу, записанному до прогона ===", flush=True)
    if fails:
        print("  КОНТРОЛЬ НЕ ПРОЙДЕН:", flush=True)
        for f in fails:
            print(f"    {f}", flush=True)
        print("  Ненулевая разность при lambda = 0 означает, что флаг режима протекает мимо\n"
              "  скринингового канала. Тогда под сомнением все сравнения режимов --- пункты\n"
              "  174-176, 307 и таблицы в них, --- а не только это плечо.", flush=True)
    else:
        print("  ПРОЙДЕН: при lambda = 0 оба двухсайтовых режима воспроизводят calshift\n"
              "  побитово на всех четырёх сидах. Флаг режима действует только через\n"
              "  скрининговый канал, и двухсайтовые плечи пункта 307 теперь\n"
              "  проконтролированы ровно так, как требуют пункты 174 и 175.", flush=True)

    # Descriptive, not part of the pass rule: what the identity arm scores.
    print(f"\n  для записи --- метрики плеча тождества (все режимы обязаны совпасть):", flush=True)
    print(f"  {'источник':>13s} {'сид':>4s} {'макро пара':>11s} {'макро ранг':>11s}", flush=True)
    for s in SEEDS:
        for src, b in [(REF, ref)] + [(n, blobs[n]) for n, _ in ARMS]:
            r = next((r for r in b["table"] if r.get("seed") == s and r.get("lambda") == 0.0), None)
            if r:
                print(f"  {src:>13s} {s:>4d} {r['MACRO']:11.4f} {r['MACRO_rho']:11.3f}", flush=True)


if __name__ == "__main__":
    main()
