"""Item 320's prediction, scored against the board. Written BEFORE the result exists.

Item 320 re-placed CYP2D6 alone, on the finding that item 317's "the blind bands are negligible"
was refuted by the board's own arithmetic -- 69 entrants post a strictly worse CYP2D6 MAE than ours
and a strictly better ST-RAE, which no band-free metric permits. The three rules below were fixed in
`verify/README.md` and committed before the upload, and they are not to be widened afterwards.

  1. RANK (sharp). The transform is affine with b > 0, so every Spearman and Kendall must be
     UNCHANGED. Measured before upload at rho = 1.000000000000 against the shipped file on all four
     enzymes, so this rule tests the BOARD rather than the arithmetic: a moved rank means something
     other than the placement changed, and that falsifies the submission rather than the prediction.
     Tolerance 0.0005 exists only because the board prints four places.
  2. EXISTENCE. CYP2D6 ST-RAE <= 0.6882 -- no worse than a strictly weaker model already achieves
     on that cell (JacksonBurns, behind us on MAE, R2, Spearman AND Kendall simultaneously). This is
     the gate that licensed the upload; failing it means the band reading was wrong in kind.
  3. MAGNITUDE. CYP2D6 ST-RAE inside [0.6146, 0.6623] AND macro strictly below 0.5760. The band is
     the world family's spread, not a hedge added afterwards.

WHAT THIS RUN ALSO MEASURES, and it is the scientific content rather than the verdict. The blind
band width is the one quantity nobody could measure; the placement's sensitivity to it was mapped
before upload. The observed CYP2D6 ST-RAE therefore IDENTIFIES that multiplier, whichever way the
rules fall. A failure that lands at 0.75x is a different fact from a failure that lands at 0.3x.

THE CONTROL, without which the rest is decoration. If the candidate is not up yet the board still
carries the old submission, and every number below would "confirm" the state before the change. So
the run refuses while `Submitted` still equals the pre-upload stamp.

It does NOT guard the other direction the way `k101_prereg308.py` came to. That file could name the
submission its rules were written for, because that submission already existed when the guard was
added. This one is written BEFORE the upload, so the timestamp it should match is not yet known and
cannot be hard-coded without inventing it. Instead the control PRINTS the stamp it is scoring, and a
reader who sees a date later than the upload knows the verdict is about a different file. That is
weaker than `k101`'s guard, and it is said here rather than papered over; once the candidate is up,
its stamp can be pinned and the guard made symmetric.

    uv run python verify/k104_prereg320.py --pull          # fetch a dated snapshot and score
    uv run python verify/k104_prereg320.py --snap <path>   # score a snapshot already on disk

Exit: 0 all three held, 1 at least one failed, 2 refused (nothing uploaded yet),
4 our row is not on the board at all -- a different event from a failed prediction.
Runtime: seconds. Writes only the snapshot it pulls, and only when asked to pull.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import RES

import argparse
import datetime
import json
import urllib.request

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]

# Every name this team has submitted under, newest first. The board stores the username AS IT WAS
# AT SUBMISSION TIME (item 314), so one team can appear under two names at once and a single
# hard-coded name silently stops matching -- which is exactly how the 21 September pull read like a
# disqualification.
NAMES = ("Lizard Wizard Gizzard", "Wizard Lizard Gizzard")
BASE = "https://openadmet-cyp-challenge.hf.space/gradio_api/call"
TABS = dict(zip(CYPS, ["partial_10", "partial_11", "partial_12", "partial_13"]))

# The board as it stood on the submission item 320 replaces, pulled 2026-09-24T12:19Z and committed
# as results/leaderboard_2026-09-24T1219Z.json.
PRE_SUBMITTED = "2026-09-24 08:12 UTC"
PRE = {
    "CYP1A2": {"ST-RAE": 0.5889, "MAE": 0.7906, "R2": 0.5093, "rho": 0.7362, "tau": 0.5329},
    "CYP2C9": {"ST-RAE": 0.4624, "MAE": 0.4694, "R2": 0.6274, "rho": 0.7976, "tau": 0.6089},
    "CYP2D6": {"ST-RAE": 0.7707, "MAE": 0.9222, "R2": 0.3859, "rho": 0.4375, "tau": 0.3011},
    "CYP3A4": {"ST-RAE": 0.4822, "MAE": 0.5358, "R2": 0.6595, "rho": 0.8026, "tau": 0.6170},
}
PRE_MACRO = {"ST-RAE": 0.5760, "rho": 0.6935, "rank": 68, "of": 230}

# Pre-registered in item 320. Not to be edited after the result arrives.
GATE_EXISTENCE = 0.6882           # rule 2: a strictly weaker model already reaches this
BAND_2D6 = (0.6146, 0.6623)       # rule 3, the world family's spread
MACRO_CEILING = 0.5760            # rule 3: strictly below
PRED_2D6 = 0.6246
PRED_MACRO = 0.5395
PRED_RANK = 60
RANK_TOL = 0.0005

# The placement's sensitivity to the blind band width, mapped BEFORE the upload. Used only to read
# the implied multiplier off the result; it changes no verdict.
SENS = [(0.4, 0.8701), (0.6, 0.7653), (0.8, 0.6709), (1.0, 0.6246), (1.5, 0.6046)]


def utcstamp():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def pull(stamp):
    """Fetch every board tab through the Space's own Gradio endpoints.

    Stamped by the caller, because a pulled board is EVIDENCE and not a scratch file.
    """
    snap = {"pulled_utc": stamp, "source": BASE + "/partial_N", "boards": {}}
    for n in range(9, 18):
        ep = f"partial_{n}"
        req = urllib.request.Request(f"{BASE}/{ep}", data=json.dumps({"data": []}).encode(),
                                     headers={"Content-Type": "application/json"})
        with urllib.request.urlopen(req, timeout=60) as r:
            eid = json.load(r)["event_id"]
        body = None
        with urllib.request.urlopen(f"{BASE}/{ep}/{eid}", timeout=180) as r:
            for line in r:
                line = line.decode().strip()
                if line.startswith("data: "):
                    d = json.loads(line[6:])
                    body = d[0] if isinstance(d, list) and d else d
                    break
        if isinstance(body, dict) and "headers" in body:
            snap["boards"][ep] = {"headers": body["headers"], "data": body.get("data") or []}
    return snap


def our_row(snap, ep):
    """Our row on one tab, under any name this team has submitted under.

    The failure branch DIAGNOSES rather than merely refusing: a bare "row not found" once read like
    a disqualification when the row was there under a transposed name.
    """
    b = snap["boards"][ep]
    rows = [r for r in b["data"] if len(r) > 1]
    for name in NAMES:
        for r in rows:
            if isinstance(r[1], str) and name in r[1]:
                return dict(zip(b["headers"], r))
    near = sorted({r[1] for r in rows if isinstance(r[1], str)
                   and set(r[1].lower().split()) & {w for n in NAMES for w in n.lower().split()}})
    # Exit 4, not 1. "We could not find ourselves" and "the prediction failed" are different
    # events and must not share a code: on 21 September a bare row-not-found read like a
    # disqualification, and a caller that cannot tell them apart repeats that reading.
    print("\n".join([
        f"{ep}: ни одно из известных имён не найдено -- сверять нечего.",
        f"    искали: {', '.join(NAMES)}",
        f"    строк на вкладке: {len(rows)}",
        (f"    похожие имена на доске: {', '.join(near)}" if near else
         "    похожих имён на доске нет"),
        "    Если команда переименована, добавьте новое имя ПЕРВЫМ в NAMES.",
    ]), file=_sys.stderr)
    _sys.exit(4)


def n_entrants(snap, ep="partial_9"):
    """How many ENTRANTS the board lists, which is the field size a place is out of.

    Named n_scored until item 327, where it returned 227 and was read as a compound count --
    the live reveal scores 375 compounds (item 294 computes the rank band at that n), so a
    function called n_scored returning the number of rows is a trap in this file's own units.
    """
    return len([x for x in snap["boards"][ep]["data"] if len(x) > 1])


def implied_multiplier(observed):
    """Read the blind band multiplier off the observed CYP2D6 cell by linear interpolation.

    The mapping is monotone DECREASING in the multiplier (wider bands forgive more), so the table is
    walked from the wide end. Outside the mapped range the answer is reported as a bound rather than
    extrapolated -- an extrapolated multiplier would be a number with no measurement behind it.
    """
    pts = sorted(SENS, key=lambda t: t[1])           # ascending ST-RAE
    if observed <= pts[0][1]:
        return f"> {pts[0][0]:.1f}x (вне размеченного диапазона)"
    if observed >= pts[-1][1]:
        return f"< {pts[-1][0]:.1f}x (вне размеченного диапазона)"
    for (m0, s0), (m1, s1) in zip(pts, pts[1:]):
        if s0 <= observed <= s1:
            f = (observed - s0) / (s1 - s0) if s1 != s0 else 0.0
            return f"~ {m0 + f * (m1 - m0):.2f}x"
    return "не определён"


def control(snap):
    """Is the board showing the submission these rules were written for?"""
    now = our_row(snap, "partial_9")["Submitted"]
    print("КОНТРОЛЬ: та ли это подача, для которой писались правила\n")
    print(f"    до перекалибровки (пункт 320)   {PRE_SUBMITTED}")
    print(f"    на доске сейчас                 {now}")
    if now == PRE_SUBMITTED:
        print("\n    ОТКАЗ: на доске всё ещё ТА ЖЕ подача. Ниже стояли бы числа состояния ДО")
        print("    изменения, и любое 'подтверждение' было бы подтверждением того, что ничего")
        print("    не менялось. Это не вердикт. Загрузите файл и прогоните снова.")
        return "рано"
    print("    -> подача сменилась, сверка осмысленна\n")
    return "ок"


def verdict(snap):
    m9 = our_row(snap, "partial_9")
    rows = {c: our_row(snap, TABS[c]) for c in CYPS}

    print("ПРАВИЛО 1 --- РАНГ. Преобразование аффинно с b > 0, ро и тау двигаться НЕ МОГУТ.\n")
    print("    %-8s %9s %9s %11s   %9s %9s %11s" %
          ("фермент", "ро до", "ро сейчас", "разность", "тау до", "тау сейчас", "разность"))
    rank_ok = True
    for c in CYPS:
        r = rows[c]
        dr = r["Spearman's ρ"] - PRE[c]["rho"]
        dt = r["Kendall's τ"] - PRE[c]["tau"]
        ok = abs(dr) <= RANK_TOL and abs(dt) <= RANK_TOL
        rank_ok &= ok
        print("    %-8s %9.4f %9.4f %+11.4f   %9.4f %9.4f %+11.4f %s" %
              (c, PRE[c]["rho"], r["Spearman's ρ"], dr,
               PRE[c]["tau"], r["Kendall's τ"], dt, "" if ok else "<- СДВИНУЛСЯ"))
    dm = m9["MA-Spearman's ρ"] - PRE_MACRO["rho"]
    rank_ok &= abs(dm) <= RANK_TOL
    print("    %-8s %9.4f %9.4f %+11.4f" % ("макро", PRE_MACRO["rho"], m9["MA-Spearman's ρ"], dm))
    print(f"\n    правило 1 -> {'ВЫПОЛНИЛОСЬ' if rank_ok else 'ПРОВАЛИЛОСЬ'}")
    if not rank_ok:
        print("    Сдвиг ранга означает, что изменилось НЕ ТОЛЬКО размещение. Это опровергает")
        print("    подачу, а не предсказание.")

    got = rows["CYP2D6"]["ST-RAE"]
    print("\n\nПРАВИЛО 2 --- СУЩЕСТВОВАНИЕ. CYP2D6 ST-RAE <= %.4f,\n"
          "    то есть не хуже того, чего уже достигла заведомо более слабая модель.\n"
          % GATE_EXISTENCE)
    ex_ok = got <= GATE_EXISTENCE
    print(f"    было {PRE['CYP2D6']['ST-RAE']:.4f}  ->  стало {got:.4f}   порог {GATE_EXISTENCE}")
    print(f"    правило 2 -> {'ВЫПОЛНИЛОСЬ' if ex_ok else 'ПРОВАЛИЛОСЬ'} "
          f"(запас {GATE_EXISTENCE - got:+.4f})")

    print("\n\nПРАВИЛО 3 --- ВЕЛИЧИНА.\n")
    in_band = BAND_2D6[0] <= got <= BAND_2D6[1]
    macro = m9["MA-ST-RAE"]
    below = macro < MACRO_CEILING
    mag_ok = in_band and below
    print(f"    CYP2D6 внутри [{BAND_2D6[0]}, {BAND_2D6[1]}]: факт {got:.4f} -> "
          f"{'ДА' if in_band else 'НЕТ'}   (предсказано {PRED_2D6})")
    print(f"    макро строго ниже {MACRO_CEILING}:        факт {macro:.4f} -> "
          f"{'ДА' if below else 'НЕТ'}   (предсказано {PRED_MACRO})")
    print(f"    правило 3 -> {'ВЫПОЛНИЛОСЬ' if mag_ok else 'ПРОВАЛИЛОСЬ'}")

    print("\n\nНАБЛЮДЕНО, НО НЕ ПРАВИЛО --- три остальных фермента не трогались,\n"
          "    так что их сдвиг измеряет ДРЕЙФ ДОСКИ, а не нас.\n")
    print("    %-8s %10s %11s %11s" % ("фермент", "ST-RAE до", "ST-RAE факт", "разность"))
    for c in CYPS:
        d = rows[c]["ST-RAE"] - PRE[c]["ST-RAE"]
        note = "  <- менялся" if c == "CYP2D6" else ""
        print("    %-8s %10.4f %11.4f %+11.4f%s" % (c, PRE[c]["ST-RAE"], rows[c]["ST-RAE"], d, note))

    n = n_entrants(snap)
    print(f"\n    место: было {PRE_MACRO['rank']} из {PRE_MACRO['of']}, "
          f"стало {m9['Rank']} из {n}  ({m9['Rank'] - PRE_MACRO['rank']:+d}, "
          f"предсказано {PRED_RANK})")
    print(f"    MAE и R2 на CYP2D6 обязаны УХУДШИТЬСЯ -- это цена размещения, а не дефект: "
          f"MAE {PRE['CYP2D6']['MAE']:.4f} -> {rows['CYP2D6']['MAE']:.4f}, "
          f"R2 {PRE['CYP2D6']['R2']:.4f} -> {rows['CYP2D6']['R²']:.4f}")

    print("\n\nЧТО ЭТОТ ПРОГОН ИЗМЕРЯЕТ ПОМИМО ВЕРДИКТА\n")
    print("    Ширина слепых полос была единственной величиной, которую нельзя было измерить.")
    print("    Наблюдённая CYP2D6 ST-RAE её ОПОЗНАЁТ, как бы ни легли правила:")
    print(f"        подразумеваемый множитель полос: {implied_multiplier(got)}")
    print("    (предрегистрированная карта: " +
          ", ".join(f"{m:.1f}x->{s:.4f}" for m, s in SENS) + ")")

    held = sum([rank_ok, ex_ok, mag_ok])
    print(f"\n\nИТОГ: выполнилось {held} из 3 предрегистрированных правил")
    print("    (1 ранг, 2 существование, 3 величина --- зафиксированы в verify/README.md,")
    print("     пункт 320, и закоммичены ДО загрузки; расширению не подлежат)")
    return held == 3


def main(argv=None):
    """Returns the exit code rather than calling sys.exit, so the guards are testable offline."""
    ap = argparse.ArgumentParser()
    ap.add_argument("--snap", default="", help="снимок доски на диске")
    ap.add_argument("--pull", action="store_true", help="стянуть свежий снимок и сохранить")
    ap.add_argument("--out", default="", help="куда писать; по умолчанию с датой и временем пулла")
    ap.add_argument("--force", action="store_true", help="разрешить перезапись существующего снимка")
    a = ap.parse_args(argv)
    if not a.snap and not a.pull:
        raise SystemExit("нужен либо --snap <путь>, либо --pull")
    if a.pull:
        # Resolve and guard the destination BEFORE the network call. A snapshot is evidence, and a
        # guard that costs a round trip to exercise is a guard nobody exercises.
        stamp = utcstamp()
        auto = RES + "leaderboard_" + stamp[:10] + "T" + stamp[11:13] + stamp[14:16] + "Z.json"
        out = a.out or auto
        if _pl.Path(out).exists() and not a.force:
            raise SystemExit(
                f"{out} уже существует, и перезаписан не будет.\n"
                f"Снимок доски --- свидетельство, а не черновик. Укажите --out другим именем,\n"
                f"или --force, если перезапись действительно нужна.")
        snap = pull(stamp)
        json.dump(snap, open(out, "w"), ensure_ascii=False, indent=1)
        print(f"снимок сохранён: {out}, вкладок {len(snap['boards'])}, "
              f"тянут {snap['pulled_utc']}\n")
    else:
        snap = json.load(open(a.snap))
        print(f"снимок: {a.snap}, тянут {snap.get('pulled_utc')}\n")
    c = control(snap)
    if c == "рано":
        return 2
    return 0 if verdict(snap) else 1


if __name__ == "__main__":
    _sys.exit(main())
