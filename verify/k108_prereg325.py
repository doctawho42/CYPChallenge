"""Item 325's prediction, scored against the board. Rules fixed and committed BEFORE the upload.

Item 325 shipped the six-member composition of item 324 (the shipped five plus a ridge head on a
frozen third-party chemprop encoder, +0.0138 of macro rank out of fold, sign 8/8). It is the first
candidate this project has placed on the board whose rank is required to MOVE: every submission
since item 308 was a monotone re-placement of one vector, so Spearman was pinned by construction
and a rank rule could only confirm that nothing had broken.

  1. RANK MUST RISE. Macro Spearman strictly above 0.6935.
  2. RANK MAGNITUDE. Macro Spearman inside [0.6995, 0.7150] --- the out-of-fold +0.0138 with a
     band from half to one and a half times it.
  3. METRIC. Macro MA-ST-RAE strictly below 0.5535.
  Prediction: macro rho 0.7073, macro ST-RAE 0.5370, place 53 of 219.

THE CONTROL, in both directions. k104_prereg320.py could only guard one: written before its upload
existed, it could refuse a board still showing the PRE-upload file but could not name the stamp its
verdict belonged to. Here both stamps are known, so the guard is symmetric --- the run refuses if
the board still carries the pre-upload submission, and warns loudly if a stamp LATER than the one
this verdict was recorded against has replaced it, because then the numbers describe a different
file. That asymmetry was flagged in k104's own docstring as something to fix once the upload landed.

WHY THE PLACE IS DECOMPOSED RATHER THAN SUBTRACTED. Item 323's defect was pricing a competitor's
move as 70/230 -> 55/219 = 15 places, comparing ranks across two snapshots of different field size.
The same trap is live here: the field grew from 219 to 227 between the prediction and the verdict.
So our previous macro is re-ranked inside the CURRENT snapshot, and the difference between that and
our new place is the part that is ours.

    uv run python verify/k108_prereg325.py --pull          # fetch a dated snapshot and score
    uv run python verify/k108_prereg325.py --snap <path>   # score a snapshot already on disk

Exit: 0 all three held, 1 at least one failed, 2 refused (pre-upload board), 4 our row absent.
Runtime: seconds. Writes only the snapshot it pulls, and only when asked to pull.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import RES

import argparse
import json

sys = _sys
_v = _pl.Path(__file__).resolve().parent
_sys.path.insert(0, str(_v))
from k104_prereg320 import NAMES, n_entrants, our_row, pull, utcstamp   # one definition, not a copy

RHO = "MA-Spearman's ρ"
TAU = "MA-Kendall's τ"
STR = "MA-ST-RAE"

PRE_UPLOAD_STAMP = "2026-09-25 01:49 UTC"   # the half-step candidate item 321 scored
VERDICT_STAMP = "2026-09-25 17:33 UTC"      # the probe candidate these rules were written for

# Our row as it stood BEFORE this upload, from results/leaderboard_2026-09-25T0558Z.json.
BEFORE = {"rho": 0.6935, "strae": 0.5535, "rank": 58, "field": 219}
PRED = {"rho": 0.7073, "strae": 0.5370, "rank": 53, "field": 219}

RULES = (
    ("1 РАНГ ДОЛЖЕН ВЫРАСТИ", "макро ро строго выше 0.6935",
     lambda d: d[RHO] > 0.6935),
    ("2 ВЕЛИЧИНА РАНГА", "макро ро внутри [0.6995, 0.7150]",
     lambda d: 0.6995 <= d[RHO] <= 0.7150),
    ("3 МЕТРИКА", "макро MA-ST-RAE строго ниже 0.5535",
     lambda d: d[STR] < 0.5535),
)


def control(snap):
    """Refuse a pre-upload board; warn if a newer file has replaced the one scored."""
    row = our_row(snap, "partial_9")
    if row is None:
        print(f"\nНАШЕЙ СТРОКИ НА ДОСКЕ НЕТ. Искал, в порядке: {', '.join(NAMES)}.")
        print("Доска ранее уже теряла строки на время (пункт 314), так что отсутствие ---")
        print("не снятие. Повторить позже, прежде чем считать это событием.")
        return None, 4
    now = row["Submitted"]
    print("\nКОНТРОЛЬ: та ли это подача, для которой писались правила\n")
    print(f"    до загрузки зонда (пункт 321)   {PRE_UPLOAD_STAMP}")
    print(f"    правила писались под            {VERDICT_STAMP}")
    print(f"    на доске сейчас                 {now}")
    if now == PRE_UPLOAD_STAMP:
        print("    -> подача НЕ СМЕНИЛАСЬ: это ещё файл пункта 321, сверять нечего")
        return row, 2
    if now != VERDICT_STAMP:
        print("    -> ВНИМАНИЕ: штамп НЕ ТОТ, под который писались правила.")
        print("       Числа ниже описывают ДРУГОЙ файл --- вероятно, более свежую загрузку.")
        print("       Вердикт пункта 325 к нему не относится.")
        return row, 1
    print("    -> совпал, сверка относится именно к предрегистрированному файлу")
    return row, 0


def decompose(snap, row):
    """Our place, and where our OLD macro would sit in THIS field. Item 323's lesson."""
    b = snap["boards"]["partial_9"]
    vals = sorted(float(dict(zip(b["headers"], r))[STR]) for r in b["data"] if len(r) > 1)
    field = n_entrants(snap)          # shared with k104 rather than counted again here
    assert len(vals) == field, f"строк с метрикой {len(vals)}, участников {field}"
    old_place = sum(1 for v in vals if v < BEFORE["strae"]) + 1
    return field, old_place


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--pull", action="store_true")
    ap.add_argument("--snap")
    a = ap.parse_args(argv)
    if a.pull:
        stamp = utcstamp()
        snap = pull(stamp)
        out = RES + "leaderboard_" + stamp[:10] + "T" + stamp[11:13] + stamp[14:16] + "Z.json"
        json.dump(snap, open(out, "w"), ensure_ascii=False, indent=1)
        print(f"снимок сохранён: {out}, вкладок {len(snap['boards'])}, тянут {stamp}")
    elif a.snap:
        snap = json.load(open(a.snap))
        print(f"снимок с диска: {a.snap}")
    else:
        raise SystemExit("нужен --pull или --snap <path>")

    row, code = control(snap)
    if code in (2, 4):
        return code

    print("\nПРАВИЛА ПУНКТА 325, зафиксированные ДО загрузки\n")
    held = 0
    for name, text, test in RULES:
        ok = bool(test(row))
        held += ok
        print(f"    правило {name}\n        {text}")
        print(f"        факт: макро ро {row[RHO]}, макро ST-RAE {row[STR]}"
              f"  ->  {'ВЫПОЛНИЛОСЬ' if ok else 'ПРОВАЛИЛОСЬ'}")
    print(f"\n    выполнилось {held} из {len(RULES)}")

    print("\nПРЕДСКАЗАНИЕ ПРОТИВ ФАКТА\n")
    print(f"    {'величина':22s}{'предсказано':>13s}{'факт':>10s}{'промах':>10s}")
    print(f"    {'макро ро':22s}{PRED['rho']:>13.4f}{row[RHO]:>10.4f}"
          f"{row[RHO]-PRED['rho']:>+10.4f}")
    print(f"    {'макро ST-RAE':22s}{PRED['strae']:>13.4f}{row[STR]:>10.4f}"
          f"{row[STR]-PRED['strae']:>+10.4f}")

    field, old_place = decompose(snap, row)
    print(f"\nМЕСТО, РАЗЛОЖЕННОЕ ВНУТРИ ОДНОГО СНИМКА (пункт 323)\n")
    print(f"    было                       {BEFORE['rank']} из {BEFORE['field']}"
          f"   (макро {BEFORE['strae']})")
    print(f"    стало                      {row['Rank']} из {field}"
          f"   (макро {row[STR]})")
    print(f"    наш прежний макро в ЭТОМ поле  {old_place} из {field}")
    ours = old_place - int(row["Rank"])
    drift = int(BEFORE["rank"]) - old_place
    print(f"    -> наших мест {ours:+d}, дрейф поля {drift:+d}"
          f"   (предсказано было {PRED['rank']} из {PRED['field']})")
    print("    Разность мест между двумя снимками разного размера мерой не является ---")
    print("    поле выросло с 219 до %d, и это единственная причина так считать." % field)

    print(f"\nИТОГ: выполнилось {held} из {len(RULES)} предрегистрированных правил "
          f"(пункт 325, закоммичены ДО загрузки; расширению не подлежат)")
    if code == 1:
        print("      НО штамп доски не тот --- см. контроль выше.")
        return 1
    return 0 if held == len(RULES) else 1


if __name__ == "__main__":
    raise SystemExit(main())
