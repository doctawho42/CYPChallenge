"""Item 308's prediction, scored against the board. Written BEFORE the result exists.

Item 308 recalibrated the regression submission by an affine map per enzyme, fixed before any
result was seen, and made a closed-form prediction. This file exists so that prediction can be
scored by a rule that was also fixed before the number arrived -- the same order as `k99_lam0.py`,
which was committed while its run was still going.

The prediction, from item 308:

  * the transform is monotone in p, so EVERY per-enzyme Spearman must be UNCHANGED. This is the
    sharp test: rank is this project's currency and the only thing the transform cannot move.
  * with the centre matched, R2 = 2*r*k' - k'^2 in closed form, so macro R2 must rise from 0.0985
    to about 0.50 (per enzyme 0.5396 / 0.6307 / 0.1911 / 0.6392).

The pass rules below are FIXED HERE, before the file goes to the board, and are not to be widened
afterwards. Their tolerances and their reasons:

  1. RANK (sharp). Every per-enzyme Spearman within 0.0005 of the 15 September value. Zero is
     expected exactly; the tolerance exists only because the board rounds to four places and a tie
     in the prediction vector could in principle move the last digit. A larger movement means
     something other than the transform changed, and that falsifies the whole submission, not just
     the arithmetic.
  2. DIRECTION. Macro R2 must rise by at least +0.30 from 0.0985. Anything less and the scale
     diagnosis was wrong.
  3. MAGNITUDE. Macro R2 must land inside [0.40, 0.60]. The band is wide on purpose: the closed
     form gives 0.5001 on the inferred spreads and 0.5031 on briford's probed ones, but item 308
     records that the MAE-to-RMSE assumption underneath is strictly VIOLATED on CYP1A2 and CYP2D6
     (the implied ST-RAE/RAE forgiveness exceeds one, which is impossible), so the per-cell numbers
     on those two are not trustworthy to two decimals. The direction and the order of magnitude are
     what was claimed; this band is what that claim deserves.

ST-RAE is reported but NOT predicted: with no credible bands on the test set there was no closed
form for it, and inventing a target afterwards would be exactly the tuning item 308 forbade.

THE CONTROL, without which this script is worthless. If the recalibrated file has not been uploaded
yet, the board still carries the OLD submission, and every number below would "confirm" the state
before the change. So the script first checks our row's submission timestamp against the recorded
pre-recalibration one and REFUSES to render a verdict if it has not moved. A query that could not
have failed is not a negative answer.

    uv run python verify/k101_prereg308.py --pull          # fetch a fresh dated snapshot and score
    uv run python verify/k101_prereg308.py --snap <path>   # score a snapshot already on disk

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

# Every name this team has submitted under, newest first. Not a convenience: the board stores
# the username AS IT WAS AT SUBMISSION TIME, so one team can appear under two names at once and
# a single hard-coded name silently stops matching. Measured on 21 September (item 314) --- the
# regression tabs carried "Lizard Wizard Gizzard" from the 19 September upload while the TDI
# tabs, untouched since 15 September, still carried "Wizard Lizard Gizzard", with byte-identical
# metrics. A live name lookup would have moved both.
NAMES = ("Lizard Wizard Gizzard", "Wizard Lizard Gizzard")
BASE = "https://openadmet-cyp-challenge.hf.space/gradio_api/call"
TABS = dict(zip(CYPS, ["partial_10", "partial_11", "partial_12", "partial_13"]))

# The board as it stood on the submission item 308 replaces. Anchors every comparison below.
PRE_SUBMITTED = "2026-09-15 10:01 UTC"
# The submission these three rules were written to score, and the ONLY one they score. Item 313
# later shipped a further affine step (2026-09-19 18:53 UTC); rules 2 and 3 are about item 308's
# R2 correction and do not transfer to it. Scoring a later upload against them would be a verdict
# on the wrong object, which is the failure this file exists to prevent, not commit.
SCORED_SUBMITTED = "2026-09-16 13:15 UTC"
PRE = {
    "CYP1A2": {"ST-RAE": 0.8853, "MAE": 1.1034, "R2": 0.0960, "rho": 0.7362},
    "CYP2C9": {"ST-RAE": 0.4862, "MAE": 0.5179, "R2": 0.5511, "rho": 0.7976},
    "CYP2D6": {"ST-RAE": 1.4043, "MAE": 1.7170, "R2": -0.8702, "rho": 0.4375},
    "CYP3A4": {"ST-RAE": 0.4837, "MAE": 0.5624, "R2": 0.6173, "rho": 0.8026},
}
PRE_MACRO = {"ST-RAE": 0.8149, "MAE": 0.9752, "R2": 0.0985, "rho": 0.6935, "rank": 103, "of": 164}

# PREDICTED by item 308's closed form. Not to be edited after the result arrives.
PRED_R2 = {"CYP1A2": 0.5396, "CYP2C9": 0.6307, "CYP2D6": 0.1911, "CYP3A4": 0.6392}
PRED_MAE = {"CYP1A2": 0.7875, "CYP2C9": 0.4698, "CYP2D6": 1.1292, "CYP3A4": 0.5461}
PRED_MACRO_R2 = 0.5001

# Pass rules, fixed before the result exists.
RANK_TOL = 0.0005
MIN_R2_RISE = 0.30
R2_BAND = (0.40, 0.60)


def utcstamp():
    return datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def pull(stamp):
    """Fetch every board tab through the Space's own Gradio endpoints.

    Stamped by the caller, because a pulled board is EVIDENCE and not a scratch file: the boards
    move on their own -- item 308 records our position sliding from 102/163 to 103/164 inside one
    afternoon while not one number of ours changed -- so a snapshot without a time is not evidence
    of anything. The caller stamps it so the filename and the record inside cannot disagree.
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

    The failure branch DIAGNOSES instead of just refusing. The 21 September pull died on
    "нашей строки нет в снимке", which is true and useless: the row was there under a
    transposed name, and the message sent the reader looking for a disqualification. A
    refusal that does not say what it looked for costs more than it saves.
    """
    b = snap["boards"][ep]
    rows = [r for r in b["data"] if len(r) > 1]
    for name in NAMES:
        for r in rows:
            if r[1] == name:
                return dict(zip(b["headers"], r))
    # Nothing matched: show what was tried, and anything that looks close enough to be a rename.
    near = sorted({r[1] for r in rows
                   if isinstance(r[1], str)
                   and set(r[1].lower().split()) & {w for n in NAMES for w in n.lower().split()}})
    msg = [f"{ep}: ни одно из известных имён не найдено -- сверять нечего.",
           f"    искали: {', '.join(NAMES)}",
           f"    строк на вкладке: {len(rows)}"]
    msg.append(f"    похожие имена на доске: {', '.join(near)}" if near else
               "    похожих имён на доске нет")
    msg.append("    Если команда переименована, добавьте новое имя ПЕРВЫМ в NAMES.")
    raise SystemExit("\n".join(msg))


def control(snap):
    """Is the board showing the submission these rules were written for?

    Two ways to get a meaningless verdict, and the original version guarded only the first.
    TOO EARLY: the recalibrated file is not up yet, so every number below would confirm the
    state BEFORE the change. TOO LATE: a further submission has replaced it, so the numbers
    describe an object these rules were never about -- item 313's step moved the same
    predictions again, and scoring it against item 308's R2 band would be a verdict on the
    wrong file. Returns a code rather than a bool so the caller can tell them apart.
    """
    now = our_row(snap, "partial_9")["Submitted"]
    print("КОНТРОЛЬ: та ли это подача, для которой писались правила\n")
    print(f"    до перекалибровки (пункт 308)   {PRE_SUBMITTED}")
    print(f"    оцениваемая этими правилами     {SCORED_SUBMITTED}")
    print(f"    на доске сейчас                 {now}")
    if now == PRE_SUBMITTED:
        print("\n    ОТКАЗ: на доске всё ещё ТА ЖЕ подача. Ниже стояли бы числа состояния ДО")
        print("    изменения, и любое 'подтверждение' было бы подтверждением того, что ничего")
        print("    не менялось. Это не вердикт. Загрузите файл и прогоните снова.")
        return "рано"
    if now != SCORED_SUBMITTED:
        print("\n    ОТКАЗ: на доске БОЛЕЕ ПОЗДНЯЯ подача. Правила 2 и 3 относятся к шагу")
        print("    пункта 308, а не к тому, что лежит на доске сейчас, и применять их к другому")
        print("    файлу -- это вердикт не о том объекте. Правило 1 (ранг) справедливо для")
        print("    любого монотонного преобразования и проверяемо отдельно.")
        print(f"    Чтобы воспроизвести исходный вердикт: --snap {rel_anchor()}")
        return "поздно"
    print("    -> это она; сверка осмысленна\n")
    return "ок"


def rel_anchor():
    """The committed snapshot on which item 308's verdict was actually rendered."""
    return "results/leaderboard_2026-09-17T1053Z.json"


def verdict(snap):
    m9 = our_row(snap, "partial_9")
    print("\nРАНГ --- резкая проверка: преобразование монотонно, ро двигаться НЕ МОЖЕТ\n")
    print("    %-8s %10s %10s %11s %s" % ("фермент", "ро до", "ро сейчас", "разность", "вердикт"))
    rank_ok = True
    for c in CYPS:
        got = our_row(snap, TABS[c])["Spearman's ρ"]
        d = got - PRE[c]["rho"]
        ok = abs(d) <= RANK_TOL
        rank_ok &= ok
        print("    %-8s %10.4f %10.4f %+11.4f %s"
              % (c, PRE[c]["rho"], got, d, "держится" if ok else "СДВИНУЛСЯ"))
    dm = m9["MA-Spearman's ρ"] - PRE_MACRO["rho"]
    macro_rank_ok = abs(dm) <= RANK_TOL
    rank_ok &= macro_rank_ok
    print("    %-8s %10.4f %10.4f %+11.4f %s"
          % ("макро", PRE_MACRO["rho"], m9["MA-Spearman's ρ"], dm,
             "держится" if macro_rank_ok else "СДВИНУЛСЯ"))
    if not rank_ok:
        print("\n    Сдвиг ранга означает, что изменилось НЕ ТОЛЬКО преобразование. Это")
        print("    опровергает подачу, а не только арифметику предсказания.")

    print("\n\nМАСШТАБ --- предсказано закрытой формой R2 = 2*r*k' - k'^2\n")
    print("    %-8s %10s %12s %10s %11s" % ("фермент", "R2 до", "R2 предск", "R2 факт", "промах"))
    for c in CYPS:
        got = our_row(snap, TABS[c])["R²"]
        print("    %-8s %+10.4f %12.4f %+10.4f %+11.4f"
              % (c, PRE[c]["R2"], PRED_R2[c], got, got - PRED_R2[c]))
    r2 = m9["MA-R²"]
    print("    %-8s %+10.4f %12.4f %+10.4f %+11.4f"
          % ("макро", PRE_MACRO["R2"], PRED_MACRO_R2, r2, r2 - PRED_MACRO_R2))

    rise = r2 - PRE_MACRO["R2"]
    dir_ok = rise >= MIN_R2_RISE
    mag_ok = R2_BAND[0] <= r2 <= R2_BAND[1]
    print(f"\n    правило 2, НАПРАВЛЕНИЕ: рост >= {MIN_R2_RISE:+.2f}; факт {rise:+.4f} -> "
          f"{'ВЫПОЛНИЛОСЬ' if dir_ok else 'ПРОВАЛИЛОСЬ'}")
    print(f"    правило 3, ВЕЛИЧИНА:    внутри [{R2_BAND[0]}, {R2_BAND[1]}]; факт {r2:.4f} -> "
          f"{'ВЫПОЛНИЛОСЬ' if mag_ok else 'ПРОВАЛИЛОСЬ'}")

    print("\n\nНАБЛЮДЕНО, НО НЕ ПРЕДСКАЗАНО (закрытой формы для ST-RAE без полос не было)\n")
    print("    %-8s %10s %11s %11s" % ("фермент", "ST-RAE до", "ST-RAE факт", "MAE факт"))
    for c in CYPS:
        r = our_row(snap, TABS[c])
        print("    %-8s %10.4f %11.4f %11.4f" % (c, PRE[c]["ST-RAE"], r["ST-RAE"], r["MAE"]))
    print("    %-8s %10.4f %11.4f %11.4f"
          % ("макро", PRE_MACRO["ST-RAE"], m9["MA-ST-RAE"], m9["MA-MAE"]))
    print(f"\n    место: было {PRE_MACRO['rank']} из {PRE_MACRO['of']}, стало "
          f"{m9['Rank']} из {len([x for x in snap['boards']['partial_9']['data'] if x])}")

    held = sum([rank_ok, dir_ok, mag_ok])
    print(f"\n\nИТОГ: выполнилось {held} из 3 предрегистрированных правил")
    print("    (1 ранг держится, 2 направление, 3 величина --- правила зафиксированы в этом")
    print("     файле до того, как результат существовал, и расширению не подлежат)")
    return held == 3


def main(argv=None):
    """Returns the exit code rather than calling sys.exit, so the guard below is testable.

    The anchor guard is the one branch that must never be exercised by actually hitting the
    network -- a test that pulls the live board to check that it refuses to overwrite is a test
    nobody runs twice. With the CLI behind a function, a harness can stub `pull` and prove both
    halves offline: that an existing target refuses BEFORE the pull, and that a free one reaches it.
    """
    ap = argparse.ArgumentParser()
    ap.add_argument("--snap", default="", help="снимок доски на диске")
    ap.add_argument("--pull", action="store_true", help="стянуть свежий снимок и сохранить")
    ap.add_argument("--out", default="", help="куда писать; по умолчанию с датой пулла")
    ap.add_argument("--force", action="store_true",
                    help="разрешить перезапись существующего снимка")
    a = ap.parse_args(argv)
    if not a.snap and not a.pull:
        raise SystemExit("нужен либо --snap <путь>, либо --pull")
    if a.pull:
        # Resolve and guard the destination BEFORE the network call, not after. A pull on 15
        # September lands on leaderboard_2026-09-15.json, which is the committed anchor every
        # comparison in item 308 is measured against; and a guard that can only be exercised by
        # first spending a network round trip is a guard nobody will exercise.
        stamp = utcstamp()
        # Date AND time in the default name. A date alone collides on the day that matters most:
        # the submission window opens at 22:01 UTC on 15 September, still the 15th in UTC, so a
        # dated default would resolve to leaderboard_2026-09-15.json -- the anchor -- and the
        # guard below would refuse at exactly the moment someone is trying to score the upload.
        # The guard stays, for an explicit --out that names an existing file; it should not be
        # something an ordinary run walks into.
        auto = RES + "leaderboard_" + stamp[:10] + "T" + stamp[11:13] + stamp[14:16] + "Z.json"
        out = a.out or auto
        if _pl.Path(out).exists() and not a.force:
            raise SystemExit(
                f"{out} уже существует, и перезаписан не будет.\n"
                f"Снимок доски --- свидетельство, а не черновик: затерев его, вы потеряете якорь,\n"
                f"к которому привязаны сравнения пункта 308. Укажите --out другим именем,\n"
                f"или --force, если перезапись действительно нужна.")
        snap = pull(stamp)
        json.dump(snap, open(out, "w"), ensure_ascii=False, indent=1)
        print(f"снимок сохранён: {out}, вкладок {len(snap['boards'])}, тянут {snap['pulled_utc']}\n")
    else:
        snap = json.load(open(a.snap))
        print(f"снимок: {a.snap}, тянут {snap.get('pulled_utc')}\n")
    c = control(snap)
    if c == "рано":
        return 2
    if c == "поздно":
        return 3
    return 0 if verdict(snap) else 1


if __name__ == "__main__":
    _sys.exit(main())
