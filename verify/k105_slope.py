"""Item 323. Recompute item 322's rank-to-ST-RAE slope from the committed boards, and
separate a competitor's realised move from the field churn underneath it.

Item 322 corrected the slope from -0.990 to -1.62 and gave two routes: a competitor's own
before-and-after (-1.625) and a cross-field fit over ranks 30-90 (-1.619). Both recompute
here exactly, so the correction itself is sound and -0.990 is dead. Two things about it are
not, and both are the same shape as defects this file already records.

FIRST, the second route is not a second route. Its 95 per cent interval is [-2.196, -1.042],
which contains -0.990, so it cannot distinguish the old value from the new one at all; and it
moves with the window far more than with anything physical -- the same width-61 window slid
twenty places up or down the board reads anywhere from -1.76 to -0.08. The window where the
two routes agree to 0.006 is one of several, and it is the one that was chosen. That is a
selected coincidence, not corroboration, and it matters because it is the only thing standing
between the causal route and being a single observation with no error bar.

SECOND, "their +0.0280 of rank bought 15 board places" counts churn as gain. Rank is a
position in a field and the field changed size: 230 entrants on 24 September, 219 on the 25th.
Scoring their OLD and NEW scores against ONE field gives 63 and 55 -- eight places. The other
seven arrived between two snapshots in which their score did not move by a digit and the board
shed twenty rows. Item 322's own last paragraph is the finding that the board transiently drops
rows; its headline number is priced across exactly such a drop.

The controls, which are the point. Every number below is read from files in results/, so the
questions that could fail loudly are: does our own row carry the numbers the scoreboard claims
(0.5535 / 0.6935 at rank 58), and can the row-matching find a name that IS there? The board
wraps usernames in an HTML anchor, so an equality test against a bare name returns nothing and
would read as "the competitor is absent" -- the same false negative item 316 recorded and item
322 corrected. Matching is therefore a case-insensitive substring over the raw cell, and the
run refuses if either row is missing rather than reporting a slope computed from one of them.

    uv run python verify/k105_slope.py

Exit: 0 both of item 322's routes reproduce, 1 either does not, 4 a row is missing.
Runtime: under a second. Fits nothing, trains nothing, writes nothing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import RES

import glob
import json

import numpy as np

# The two snapshots that bracket the competitor's 24 September upload, and the current board.
BEFORE = "leaderboard_2026-09-24T1219Z.json"
AFTER = "leaderboard_2026-09-25T0644Z.json"
MID = "leaderboard_2026-09-25T0150Z.json"

RIVAL = "jacksonburns"                 # substring, case-insensitive: the cell is an <a> tag
OURS = ("Lizard Wizard Gizzard", "Wizard Lizard Gizzard")

# Item 322's numbers, to be reproduced rather than trusted.
CLAIM_CAUSAL = -1.625
CLAIM_FIELD = -1.619
CLAIM_PLACES = 15


def board(fn):
    """The macro regression board: (pulled_utc, rows). partial_9 carries MA-* columns."""
    with open(RES + fn) as fh:
        d = json.load(fh)
    return d["pulled_utc"], d["boards"]["partial_9"]["data"]


def find(rows, frag):
    """Rows whose username CELL contains frag, case-insensitively.

    Not equality: the board returns the username wrapped in an anchor tag, so `cell == name`
    is a query that cannot succeed and would report every entrant as withdrawn.
    """
    return [r for r in rows if frag.lower() in str(r[1]).lower()]


def main():
    tb, Bb = board(BEFORE)
    tm, Bm = board(MID)
    ta, Ba = board(AFTER)

    # --- control: our own row must be there and must carry the scoreboard's numbers -------
    ours = [r for r in Ba if any(n.lower() in str(r[1]).lower() for n in OURS)]
    if not ours:
        print("ОТКАЗ: нашей строки нет на доске 09-25; сравнивать не с чем")
        return 4
    o = ours[0]
    print(f"контроль, наша строка {ta}: место {o[0]} из {len(Ba)}, "
          f"ST-RAE {o[4]}, ро {o[7]}")
    if abs(o[4] - 0.5535) > 5e-4 or abs(o[7] - 0.6935) > 5e-4:
        print("  ВНИМАНИЕ: не совпало с табло в verify/README.md (0.5535 / 0.6935)")

    rb, ra = find(Bb, RIVAL), find(Ba, RIVAL)
    rb = [r for r in rb if "burns" in str(r[1]).lower()]
    ra = [r for r in ra if "burns" in str(r[1]).lower()]
    if not rb or not ra:
        print("ОТКАЗ: строка соперника отсутствует в одном из снимков")
        return 4
    rb, ra = rb[0], ra[0]

    # --- route 1: the causal before-and-after ---------------------------------------------
    dst, drho = ra[4] - rb[4], ra[7] - rb[7]
    s_causal = dst / drho
    print()
    print(f"1. причинный маршрут (ход соперника)")
    print(f"   {tb}: ST-RAE {rb[4]:.4f}  ро {rb[7]:.4f}  место {rb[0]} из {len(Bb)}")
    print(f"   {ta}: ST-RAE {ra[4]:.4f}  ро {ra[7]:.4f}  место {ra[0]} из {len(Ba)}")
    print(f"   dST-RAE {dst:+.4f} / dро {drho:+.4f} = наклон {s_causal:+.4f} "
          f"(пункт 322: {CLAIM_CAUSAL})")

    # --- route 2: the cross-field fit, and its dependence on the window -------------------
    st = np.array([r[4] for r in Ba], float)
    rho = np.array([r[7] for r in Ba], float)
    rank = np.array([r[0] for r in Ba], int)

    def fit(lo, hi):
        m = (rank >= lo) & (rank <= hi)
        return int(m.sum()), float(np.polyfit(rho[m], st[m], 1)[0])

    n30, s30 = fit(30, 90)
    print()
    print(f"2. поперёк поля, ранги 30-90 (n={n30}): наклон {s30:+.4f} "
          f"(пункт 322: {CLAIM_FIELD})")
    m = (rank >= 30) & (rank <= 90)
    res = st[m] - np.polyval(np.polyfit(rho[m], st[m], 1), rho[m])
    sx = np.sum((rho[m] - rho[m].mean()) ** 2)
    se = float(np.sqrt(np.sum(res ** 2) / (m.sum() - 2) / sx))
    print(f"   se {se:.4f}, 95% [{s30 - 1.96 * se:+.4f}, {s30 + 1.96 * se:+.4f}] "
          f"-- содержит ли -0.990: {'ДА' if s30 - 1.96*se <= -0.990 <= s30 + 1.96*se else 'нет'}")
    print("   то же окно шириной 61, сдвинутое по доске:")
    print("   " + " ".join(f"{lo:>6d}-{lo+60}" for lo in range(10, 130, 10)))
    print("   " + " ".join(f"{fit(lo, lo + 60)[1]:>9.2f}" for lo in range(10, 130, 10)))

    # --- the places, scored against ONE field ---------------------------------------------
    others = np.array([r[4] for r in Ba if r[1] != ra[1]], float)
    place = lambda x: int((others < x).sum()) + 1
    p_old, p_new = place(rb[4]), place(ra[4])
    print()
    print(f"3. места, оба счёта на ОДНОМ поле ({ta}, n={len(Ba)})")
    print(f"   старый счёт {rb[4]:.4f} -> место {p_old}")
    print(f"   новый счёт  {ra[4]:.4f} -> место {p_new}")
    print(f"   приписывается ходу: {p_old - p_new} мест (пункт 322: {CLAIM_PLACES})")
    rm = [r for r in find(Bm, RIVAL) if "burns" in str(r[1]).lower()]
    if rm:
        rm = rm[0]
        print(f"   контроль churn: {tm} ST-RAE {rm[4]:.4f} ро {rm[7]:.4f} место {rm[0]} из {len(Bm)}"
              f" -- счёт тот же, поле {len(Bm)} -> {len(Ba)}, место {rm[0]} -> {ra[0]}")

    # --- what one place costs us, read off the field rather than through a slope -----------
    mine = np.array([r[4] for r in Ba if abs(r[4] - o[4]) > 1e-12], float)
    here = int((mine < o[4]).sum()) + 1
    above = np.sort(mine)[here - 2]
    cost = o[4] - above
    print()
    print(f"4. одно место отсюда: сосед сверху {above:.4f}, нужно {cost:+.4f} ST-RAE")
    print(f"   {'наклон':>28s} {'одно место в ро':>18s}")
    for nm, sl in (("-1.62 (причинный)", 1.62), ("-0.97 (локальный, 43-73)", 0.97),
                   ("-0.990 (старый)", 0.990)):
        print(f"   {nm:>28s} {cost / sl:18.4f}")
    print("   пол макро-ранга 0.0036 -- ниже него результат не результат")

    ok = abs(s_causal - CLAIM_CAUSAL) < 5e-3 and abs(s30 - CLAIM_FIELD) < 5e-3
    print()
    print("оба маршрута пункта 322 воспроизвелись" if ok else "РАСХОЖДЕНИЕ с пунктом 322")
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
