"""The two NCGC caches the journal never reported, read out of their own tables. Reads only.

WHY A SECOND FILE. `verify/k129_ncgcread.py` is the measurement: it re-scores the cached
predictions and nests the arm over the composition, and it needs the challenge labels and bands to
do either. This one needs nothing but the caches, so the arithmetic that says there is something to
go and measure can be reproduced on any checkout -- which is the whole reason this finding sat
unread, and it would be absurd to answer that with a number nobody can re-derive without the data.

WHAT IS HERE. `src/ablncgc.py` wrote four caches into `results/preds/`. Item 153 reports
`oof_ncgc.json` and item 157 reports `oof_ncgc_w.json`. Neither `oof_ncgc_w123.json` -- the
`--seeds 1,2,3` replication item 157 called "queued and what settles it" -- nor `oof_ncgc_per.json`
appears anywhere in `verify/README.md`, `METHOD.md` or `docs/`. Item 187 is the only trace of the
second: it records that the per-enzyme arm "died instantly three separate times" on a parser defect,
"twice inside a queue whose wrapper reported success", and that the parser was replaced with a
regular expression. The repaired run's output is that file.

THE CURRENCY IS THE PROJECT'S OWN, which is the first thing this file establishes rather than
assumes. `src/ablncgc.py` scores per enzyme on that enzyme's masked vector -- `spearmanr` for rank
and `strae` after `fit_apply` for the pair, macro as the mean of four -- and that is
`verify/k106_probe.py:score()` exactly. So these numbers sit in the same units as item 324's +0.0138
and item 346's +0.0047 and need no conversion before being compared with them.

THE CONTROLS THAT MUST SUCCEED, and neither needs a label.
  1. THE REPORTED CACHES MUST REPRODUCE THE PUBLISHED LADDER. Item 157's four rows are hard-coded
     below and read back out of `oof_ncgc_w.json` and `oof_ncgc.json`. If they do not match, this
     checkout's caches are not the ones the journal was written from and nothing else here counts.
  2. EVERY PREDICTION VECTOR MUST HAVE ITS ENZYME'S ROW COUNT, 1412 / 1285 / 1493 / 2335. A cache
     whose vectors are the wrong length is misaligned against `rows.csv`, and that is checkable
     without the labels even though the scoring is not.

WHAT THIS FILE CANNOT DECIDE, stated here because the number is large and the temptation is to read
it as a result. The per-enzyme gain is measured over the per-enzyme MEMBER, not over the shipped
composition, and `docs/teammate_asks.md` records the ceiling that eats member-level gains: a delta
sinks by averaging (~delta/N) and again by redundancy, 3 to 22x on the fused enzymes. Items 176,
182, 191 and 194 are four consecutive standalone gains that arrived at the ensemble worth nothing,
and items 290 and 319 are a real channel dying the same way. So this is a candidate; `k129` is the
test.

    uv run python verify/k130_ncgctable.py        (stdlib only; plain python3 also runs it)

Reads results/preds/oof_ncgc{,_w,_w123,_per}.json. Writes results/logs/k130_ncgctable.json.
Exit: 0 ran, 2 a control failed, 3 a cache is missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))

import json

ROOT = _pl.Path(__file__).resolve().parents[1]
PRED = ROOT / "results" / "preds"
OUT = ROOT / "results" / "logs" / "k130_ncgctable.json"

CYPS = ("CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4")
ROWS = {"CYP1A2": 1412, "CYP2C9": 1285, "CYP2D6": 1493, "CYP3A4": 2335}
FLOOR = {"CYP1A2": 0.0061, "CYP2C9": 0.0071, "CYP2D6": 0.0049, "CYP3A4": 0.0033}
FLOOR_MACRO = 0.0036        # item 165, four seeds
FLOOR_STRAE = 0.007         # item 70, fixed seed

# Item 157's published ladder, seed 0, macro rank. Hard-coded so the caches must reproduce the
# journal rather than the journal being re-derived from whatever the caches now hold.
ITEM157 = {"база": 0.5462, "+NCGC активные w0.3": 0.5392,
           "+NCGC активные w0.1": 0.5510, "+NCGC всё, с цензурой w0.3": 0.5483}
ITEM153 = {"база": 0.5462, "+NCGC активные": 0.5286}

FAMILIES = [
    ("поферментно", ["oof_ncgc_per.json"], "поферментно",
     "поферментно+NCGC w0.1", "поферментно+NCGC w0.1, перемешанный", False),
    ("пул", ["oof_ncgc_w.json", "oof_ncgc_w123.json"], "база",
     "+NCGC активные w0.1", None, True),
]


def load(names):
    table, preds = [], {}
    for n in names:
        p = PRED / n
        if not p.exists():
            print("нет %s" % p, file=_sys.stderr)
            raise SystemExit(3)
        d = json.load(open(p))
        table += d["table"]
        for k, v in d["preds"].items():
            preds.setdefault(k, v)
    return table, preds


def row(table, seed, arm):
    return next((r for r in table if r["seed"] == seed and r["рука"] == arm), None)


def stats(v):
    n = len(v)
    m = sum(v) / n
    sd = (sum((x - m) ** 2 for x in v) / (n - 1)) ** 0.5 if n > 1 else float("nan")
    return m, sd, sum(1 for x in v if x > 0), n


def main():
    # ---- control 1: the reported caches reproduce the published ladder ----
    for names, want, which in ((["oof_ncgc_w.json"], ITEM157, "пункт 157"),
                               (["oof_ncgc.json"], ITEM153, "пункт 153")):
        t, _ = load(names)
        for arm, v in want.items():
            r = row(t, 0, arm)
            if r is None or abs(r["MACRO rho"] - v) > 1e-4:
                print("КОНТРОЛЬ 1 ПРОВАЛЕН: %s, рука %r: кэш %r против опубликованного %.4f"
                      % (which, arm, None if r is None else r["MACRO rho"], v), file=_sys.stderr)
                return 2
        print("контроль 1: %s воспроизведён из кэша, %d рук" % (which, len(want)))

    out = {"meta": {"floor_macro": FLOOR_MACRO, "floor_per": FLOOR, "floor_strae": FLOOR_STRAE},
           "families": {}}
    checked = 0

    for fam, files, base_arm, panel_arm, shuf_arm, reported in FAMILIES:
        table, preds = load(files)
        seeds = sorted({r["seed"] for r in table if row(table, r["seed"], panel_arm)})
        print("\n" + "=" * 90)
        print("%s  (%s)   %s" % (fam, ", ".join(files),
                                 "в журнале есть" if reported else "В ЖУРНАЛЕ НЕТ"))
        print("=" * 90)

        # ---- control 2: every vector has its enzyme's row count ----
        for s in seeds:
            for arm in (base_arm, panel_arm, shuf_arm):
                if arm is None:
                    continue
                for c in CYPS:
                    k = "%d|%s|%s" % (s, arm, c)
                    if k not in preds:
                        continue
                    if len(preds[k]) != ROWS[c]:
                        print("КОНТРОЛЬ 2 ПРОВАЛЕН: %r даёт %d значений, у %s строк %d"
                              % (k, len(preds[k]), c, ROWS[c]), file=_sys.stderr)
                        return 2
                    checked += 1

        acc = {"rank": {c: [] for c in CYPS}, "pair": {c: [] for c in CYPS},
               "rank_macro": [], "pair_macro": [], "rank_shuf": [], "pair_shuf": []}
        print("  %-5s %9s %9s %9s %9s %10s %10s" % ("сид", "1A2", "2C9", "2D6", "3A4",
                                                    "MACRO ранг", "MACRO пара"))
        for s in seeds:
            b, a = row(table, s, base_arm), row(table, s, panel_arm)
            if b is None or a is None:
                continue
            d = {c: a["%s rho" % c] - b["%s rho" % c] for c in CYPS}
            dm = a["MACRO rho"] - b["MACRO rho"]
            dp = a["MACRO пара"] - b["MACRO пара"]
            for c in CYPS:
                acc["rank"][c].append(d[c])
                acc["pair"][c].append(a["%s пара" % c] - b["%s пара" % c])
            acc["rank_macro"].append(dm)
            acc["pair_macro"].append(dp)
            print("  %-5d %+9.4f %+9.4f %+9.4f %+9.4f %+10.4f %+10.4f"
                  % (s, d["CYP1A2"], d["CYP2C9"], d["CYP2D6"], d["CYP3A4"], dm, dp))
            if shuf_arm is not None:
                sh = row(table, s, shuf_arm)
                if sh is not None:
                    acc["rank_shuf"].append(sh["MACRO rho"] - b["MACRO rho"])
                    acc["pair_shuf"].append(sh["MACRO пара"] - b["MACRO пара"])

        m, sd, pos, n = stats(acc["rank_macro"])
        pm, psd, ppos, _ = stats(acc["pair_macro"])
        print("  %-5s %+9.4f %+9.4f %+9.4f %+9.4f %+10.4f %+10.4f"
              % ("сред", *[stats(acc["rank"][c])[0] for c in CYPS], m, pm))
        print("  %-5s %9s %9s %9s %9s %8d/%d %8d/%d"
              % ("знак", *["%d/%d" % (stats(acc["rank"][c])[2], n) for c in CYPS],
                 pos, n, n - ppos, n))
        print("  %-5s %9.4f %9.4f %9.4f %9.4f %10.4f %10s"
              % ("пол", *[FLOOR[c] for c in CYPS], FLOOR_MACRO, "%.3f" % FLOOR_STRAE))
        print("  sd макро: ранг %.4f, пара %.4f" % (sd, psd))
        if acc["rank_shuf"]:
            sm = stats(acc["rank_shuf"])
            sp = stats(acc["pair_shuf"])
            print("  КОНТРОЛЬ перемешанной панели: ранг %+.4f (знак %d/%d), пара %+.4f"
                  % (sm[0], sm[2], sm[3], sp[0]))
            if sm[0] >= m or sp[0] <= pm:
                print("КОНТРОЛЬ 3 ПРОВАЛЕН: перемешанная панель не хуже настоящей",
                      file=_sys.stderr)
                return 2
        out["families"][fam] = {"files": files, "reported_in_journal": reported,
                                "seeds": seeds, "rank_macro": acc["rank_macro"],
                                "pair_macro": acc["pair_macro"],
                                "rank_per": {c: acc["rank"][c] for c in CYPS},
                                "pair_per": {c: acc["pair"][c] for c in CYPS},
                                "rank_shuf": acc["rank_shuf"], "pair_shuf": acc["pair_shuf"],
                                "rank_mean": m, "rank_sd": sd, "rank_sign": [pos, n],
                                "pair_mean": pm, "pair_sd": psd}

    print("\nконтроль 2: %d векторов, у каждого длина своего фермента" % checked)
    print("\nЭТО ЧЛЕН-УРОВНЕВЫЕ ЧИСЛА. Решает вложение в состав --- verify/k129_ncgcread.py,")
    print("которому нужны метки и полосы; см. docs/transfer_axis.md.")
    OUT.parent.mkdir(parents=True, exist_ok=True)
    json.dump(out, open(OUT, "w"), ensure_ascii=False, indent=1)
    print("\nзапись: %s" % OUT.relative_to(ROOT))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
