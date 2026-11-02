"""What first place actually costs, decomposed out of the committed board snapshots.

WHY THIS EXISTS. The project has nineteen leaderboard snapshots in `results/` and has used them
for exactly two things: reading our own row's movement (items 313, 319, 344) and reading three
entrants' published code (item 316). Nobody has used the other 272 rows as a SAMPLE. They are one,
and they answer a question no out-of-fold measurement can: of the distance between us and first
place, how much is ordering and how much is where we put the predictions inside the bands?

That matters because this project has spent most of its effort on the second. The dead zone, the
affine pair, the delta tilt, the monotone placement budget, the threshold rules -- items 128, 148,
256, 277, 282-285, 294, 313 and 322 are all placement. If placement is already better than the
teams above us, every further idea in that family is worth approximately nothing and the entire
remaining gap is the ordering axis, which is where representation lives.

HOW THE DECOMPOSITION WORKS. For each entrant the board publishes both a ST-RAE and a Spearman rho
on the same predictions. Rho sees only the ordering; ST-RAE sees the ordering AND the placement,
because a prediction inside a compound's published band scores zero whatever its rank. So regress
ST-RAE on rho across the competitive band and read each team's residual: a negative residual is a
team getting more out of its ordering than the field does, which is placement skill. The slope of
that regression is NOT the project's own -1.62 (items 322, 323) and must not be confused with it --
that one is a within-team transform, measured on one before-and-after; this one is cross-sectional
and is contaminated by teams that are good at both. Both are reported, labelled, and the gap is
converted through each.

THE CONTROLS THAT MUST SUCCEED, and the first one is the reason this file can name enzymes at all.
The board's five regression tabs are called `partial_9` to `partial_13` and nothing in them says
which enzyme is which. Item 313 published our four per-enzyme ST-RAE values for the 21 September
pull -- CYP1A2 0.5889, CYP2C9 0.4624, CYP3A4 0.4822, CYP2D6 0.7707, macro 0.5760 -- so the 21
September snapshot identifies every tab EXACTLY, by matching our own row against numbers committed
to the journal before this file existed. All five must match to 1e-6 or the identification is a
guess and this exits 2. Second: our row in the latest snapshot must carry item 344's last recorded
upload, ST-RAE 0.5357 and rho 0.7062, or the row found is not ours. Third, item 313's own control:
a username that is certainly absent must return zero rows, because a substring matcher that
silently matches three usernames is how item 313's reader broke.

    uv run python verify/k128_boardgap.py        (stdlib only; plain python3 also runs it)

Reads results/leaderboard_*.json. Writes results/logs/k128_boardgap.json.
Exit: 0 ran, 2 a control failed, 3 the snapshots are missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))

import glob
import json
import os
import re

ROOT = _pl.Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "logs" / "k128_boardgap.json"

US = "Lizard Wizard Gizzard"
ABSENT = "этой команды на доске нет точно"

# Item 344's last recorded upload, and item 313's committed per-enzyme table for 21.09.
US_LAST = {"strae": 0.5357, "rho": 0.7062}
ITEM313 = {"CYP1A2": 0.5889, "CYP2C9": 0.4624, "CYP3A4": 0.4822, "CYP2D6": 0.7707}
ITEM313_MACRO = 0.5760
SLOPE_WITHIN = -1.62        # item 322, corrected from -0.990 by item 323
PLACE_PRICE = 0.0031        # item 346: ST-RAE per board place
BEST_LEVER = 0.0127         # item 327: the frozen probe's REALISED board rho gain


def clean(x):
    return re.sub(r"<[^>]+>", "", str(x)).strip()


def numeric(rows, si, ri):
    out = []
    for r in rows:
        if isinstance(r[si], (int, float)) and isinstance(r[ri], (int, float)):
            out.append({"rank": r[0], "name": clean(r[1]), "strae": r[si], "rho": r[ri]})
    return out


def regression_boards(boards):
    """The tabs that carry a ST-RAE column, with the indices of ST-RAE and rho."""
    out = {}
    for bk, b in boards.items():
        h = [str(x) for x in b["headers"]]
        si = next((i for i, x in enumerate(h) if x.endswith("ST-RAE")), None)
        ri = next((i for i, x in enumerate(h) if "Spearman" in x), None)
        if si is not None and ri is not None:
            out[bk] = (b["data"], si, ri)
    return out


def ols(pts):
    n = len(pts)
    mx = sum(p["rho"] for p in pts) / n
    my = sum(p["strae"] for p in pts) / n
    sxx = sum((p["rho"] - mx) ** 2 for p in pts)
    b = sum((p["rho"] - mx) * (p["strae"] - my) for p in pts) / sxx
    a = my - b * mx
    r = [p["strae"] - (a + b * p["rho"]) for p in pts]
    sd = (sum(x * x for x in r) / (n - 2)) ** 0.5
    return a, b, sd, n


def main():
    snaps = sorted(glob.glob(str(ROOT / "results" / "leaderboard_*.json")))
    if not snaps:
        print("нет снимков доски в results/", file=_sys.stderr)
        return 3
    latest, pull21 = snaps[-1], None
    for s in snaps:
        if "2026-09-21T1658Z" in s:
            pull21 = s
    if pull21 is None:
        print("нет снимка 21.09, по которому опознаются вкладки", file=_sys.stderr)
        return 3

    print("снимки: %d, последний %s" % (len(snaps), os.path.basename(latest)))

    # ---- control 1: identify every per-enzyme tab against item 313's committed table ----
    tabs = regression_boards(json.load(open(pull21))["boards"])
    ident, macro_tab = {}, None
    for bk, (rows, si, ri) in tabs.items():
        ours = [r for r in rows if clean(r[1]) == US]
        if not ours:
            continue
        v = ours[0][si]
        for e, x in ITEM313.items():
            if isinstance(v, (int, float)) and abs(v - x) < 1e-6:
                ident[bk] = e
        if isinstance(v, (int, float)) and abs(v - ITEM313_MACRO) < 1e-6:
            macro_tab = bk
    if sorted(ident.values()) != sorted(ITEM313) or macro_tab is None:
        print("КОНТРОЛЬ 1 ПРОВАЛЕН: вкладки не опознаются по пункту 313: %r, макро %r"
              % (ident, macro_tab), file=_sys.stderr)
        return 2
    print("контроль 1: вкладки опознаны по пункту 313 -- " +
          ", ".join("%s=%s" % (k, v) for k, v in sorted(ident.items())) +
          ", макро=%s" % macro_tab)

    # ---- control 3: a name that is certainly absent must match nothing ----
    tabs_l = regression_boards(json.load(open(latest))["boards"])
    rows_m, si, ri = tabs_l[macro_tab]
    if [r for r in rows_m if clean(r[1]) == ABSENT]:
        print("КОНТРОЛЬ 3 ПРОВАЛЕН: отсутствующее имя нашлось", file=_sys.stderr)
        return 2
    pts = numeric(rows_m, si, ri)
    us = [p for p in pts if p["name"] == US]
    if len(us) != 1:
        print("КОНТРОЛЬ: наша строка найдена %d раз" % len(us), file=_sys.stderr)
        return 2
    us = us[0]
    print("контроль 3: отсутствующее имя даёт 0 строк; наша строка найдена ровно одна")

    # ---- control 2: our row carries item 344's last upload ----
    if abs(us["strae"] - US_LAST["strae"]) > 1e-6 or abs(us["rho"] - US_LAST["rho"]) > 1e-6:
        print("КОНТРОЛЬ 2 ПРОВАЛЕН: строка %r не несёт числа пункта 344 %r"
              % (us, US_LAST), file=_sys.stderr)
        return 2
    print("контроль 2: наша строка несёт ST-RAE %.4f и ро %.4f пункта 344"
          % (us["strae"], us["rho"]))

    first = min(pts, key=lambda p: p["rank"])
    band = [p for p in pts if p["strae"] < 0.9 and p["rho"] > 0.5]
    a, b, sd, n = ols(band)
    resid = lambda p: p["strae"] - (a + b * p["rho"])

    print("\nМАКРО: нас %d-е из %d; ST-RAE %.4f, ро %.4f.  Первое место (%s): %.4f / %.4f"
          % (us["rank"], len(pts), us["strae"], us["rho"], first["name"], first["strae"], first["rho"]))
    print("  отрыв: ST-RAE %+.4f, ро %+.4f" % (us["strae"] - first["strae"], first["rho"] - us["rho"]))

    print("\nРАЗЛОЖЕНИЕ. Поперечная связь по конкурентной полосе (n=%d):" % n)
    print("  ST-RAE = %.4f %+.4f * ро,  sd остатка %.4f" % (a, b, sd))
    print("  НЕ путать с наклоном %.2f пункта 322: тот внутрикомандный, этот поперечный"
          % SLOPE_WITHIN)
    print("  остаток (ниже = лучше кладём внутрь полос):  мы %+.4f (%.2f sd),  первое место %+.4f"
          % (resid(us), resid(us) / sd, resid(first)))
    ordering = b * (first["rho"] - us["rho"])
    print("  объясняется ПОРЯДКОМ (ро %+.4f на поперечном наклоне): %+.4f"
          % (first["rho"] - us["rho"], ordering))
    print("  остаётся на РАСКЛАДКУ:                                 %+.4f"
          % ((first["strae"] - us["strae"]) - ordering))
    need_x = (first["strae"] - resid(us) - a) / b - us["rho"]
    need_w = (us["strae"] - first["strae"]) / abs(SLOPE_WITHIN)
    print("\nСКОЛЬКО РО НУЖНО, чтобы при НАШЕЙ раскладке выйти на ST-RAE первого места:")
    print("  поперечным наклоном  %+.4f ро" % need_x)
    print("  нашим наклоном (322) %+.4f ро" % need_w)
    print("  в рычагах размера пункта 327 (+%.4f ро на доске): от %.1f до %.1f"
          % (BEST_LEVER, need_x / BEST_LEVER, need_w / BEST_LEVER))

    # ---- per enzyme ----
    print("\nПО ФЕРМЕНТАМ (последний снимок):")
    per, gaps = {}, {}
    for bk, e in sorted(ident.items(), key=lambda kv: kv[1]):
        rows, s2, r2 = tabs_l[bk]
        p2 = numeric(rows, s2, r2)
        o = next(p for p in p2 if p["name"] == US)
        f1 = min(p2, key=lambda p: p["rank"])
        bestrho = max(p2, key=lambda p: p["rho"])
        per[e] = {"ours_rank": o["rank"], "ours_strae": o["strae"], "ours_rho": o["rho"],
                  "first_name": f1["name"], "first_strae": f1["strae"], "first_rho": f1["rho"],
                  "best_rho": bestrho["rho"],
                  "gap_strae": round(o["strae"] - f1["strae"], 4),
                  "gap_rho": round(bestrho["rho"] - o["rho"], 4), "n": len(p2)}
        gaps[e] = per[e]["gap_strae"]
        print("  %-7s нас #%3d  ST-RAE %.4f ро %.4f  |  1-е %.4f ро %.4f  (лучшее ро %.4f)"
              % (e, o["rank"], o["strae"], o["rho"], f1["strae"], f1["rho"], bestrho["rho"]))
        print("          отрыв ST-RAE %+.4f, отрыв ро до лучшего %+.4f"
              % (per[e]["gap_strae"], per[e]["gap_rho"]))
    tot = sum(gaps.values())
    print("\n  доля макро-отрыва по ферментам (среднее отрывов %.4f против макро %.4f):"
          % (tot / 4, us["strae"] - first["strae"]))
    for e, g in sorted(gaps.items(), key=lambda kv: -kv[1]):
        print("    %-7s %.4f  = %4.1f %% отрыва" % (e, g, 100 * g / tot))

    # ---- is first place a moving target? all snapshots, same tab ----
    print("\nТРАЕКТОРИЯ ПОЛЯ (вкладка %s):" % macro_tab)
    print("  %-20s %5s %8s %8s %8s %8s" % ("снимок", "n", "1-е", "2-е", "10-е", "мы"))
    traj = []
    for s_ in snaps:
        try:
            bb = json.load(open(s_))["boards"]
        except Exception:
            continue
        if macro_tab not in bb:
            continue
        rr, s3, r3 = regression_boards(bb)[macro_tab]
        pp = numeric(rr, s3, r3)
        if not pp:
            continue
        byr = sorted(pp, key=lambda p: p["rank"])
        mine = [p for p in pp if p["name"] == US]
        tag = re.search(r"leaderboard_(.+)\.json", s_).group(1)
        row = {"snapshot": tag, "n": len(pp), "first": byr[0]["strae"],
               "second": byr[1]["strae"] if len(byr) > 1 else None,
               "tenth": byr[9]["strae"] if len(byr) > 9 else None,
               "us": mine[0]["strae"] if mine else None}
        traj.append(row)
        print("  %-20s %5d %8.4f %8s %8s %8s"
              % (tag, row["n"], row["first"],
                 "%.4f" % row["second"] if row["second"] is not None else "—",
                 "%.4f" % row["tenth"] if row["tenth"] is not None else "—",
                 "%.4f" % row["us"] if row["us"] is not None else "—"))
    if len(traj) > 1:
        d1 = traj[-1]["first"] - traj[0]["first"]
        print("  первое место сдвинулось на %+.4f ST-RAE за %d снимков -- %s"
              % (d1, len(traj), "цель почти неподвижна" if abs(d1) < 0.05 else "цель движется"))

    # ---- a control on the price per place, against item 346 ----
    near = sorted([p for p in pts if abs(p["strae"] - us["strae"]) < 0.06], key=lambda p: p["strae"])
    price = None
    if len(near) > 3:
        dp = near[-1]["rank"] - near[0]["rank"]
        price = (near[-1]["strae"] - near[0]["strae"]) / max(dp, 1)
        print("\nцена места рядом с нами: %.4f ST-RAE (пункт 346 писал %.4f)" % (price, PLACE_PRICE))

    json.dump({"snapshot": os.path.basename(latest), "tabs": ident, "macro_tab": macro_tab,
               "us": us, "first": first, "n_board": len(pts),
               "fit": {"intercept": a, "slope_cross": b, "resid_sd": sd, "n": n,
                       "slope_within_item322": SLOPE_WITHIN},
               "resid_us": resid(us), "resid_first": resid(first),
               "gap_strae": us["strae"] - first["strae"], "gap_rho": first["rho"] - us["rho"],
               "gap_from_ordering": ordering,
               "gap_from_placement": (first["strae"] - us["strae"]) - ordering,
               "rho_needed_cross": need_x, "rho_needed_within": need_w,
               "levers_item327": [need_x / BEST_LEVER, need_w / BEST_LEVER],
               "per_enzyme": per, "price_per_place": price, "trajectory": traj},
              open(OUT, "w"), ensure_ascii=False, indent=1)
    print("\nзапись: %s" % OUT.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
