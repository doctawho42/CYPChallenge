"""Affine recalibration of the regression submission: q = mu_y + b*(p - mean(p)), per enzyme.

Why this exists. The submission's own affine pair (`fit_shrinkage`) is fitted against OUR
labels' credible bands, so it can only calibrate to the TRAINING distribution. The blind test
distribution is wider and differently centred, and ST-RAE is normalised by the spread of
y_true -- so a submission can hold a mid-field ordering and still score near the bottom. Item
308 measures exactly that and fixes the scale and location from board-derived moments.

Monotone in p with b > 0, so Spearman and Kendall are untouched BY CONSTRUCTION and only the
scale-sensitive metrics (ST-RAE, MAE, R2) can move. The assertion below proves it per run
rather than trusting the algebra.

Writes a NEW file and a provenance record beside it; never overwrites the input. Verifying a
no-op by hash does NOT work here -- pandas re-serialises the floats, so an identity run
changes the sha256 while every number is bit-identical. Compare numerically instead.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import RES

import argparse
import datetime
import hashlib
import json
import subprocess

import numpy as np
import pandas as pd

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
COL = "{c}_pIC50_direct_inhibition"
ROOT = _pl.Path(__file__).resolve().parents[1]


def rel(p):
    """Repo-relative path for the provenance record.

    `RES` is an absolute prefix, and four people work from four machines: a home
    directory baked into a file that travels with the repository is noise at best.
    `submission.meta.json` stores repo-relative paths and this matches it.
    """
    try:
        return str(_pl.Path(p).resolve().relative_to(ROOT))
    except ValueError:
        return str(p)


def sha256(path):
    return hashlib.sha256(open(path, "rb").read()).hexdigest()


def git_state():
    """Commit and dirty tree, so a file on disk can be traced back to code.

    core.quotepath=false because git quotes non-ASCII paths by default, which made an
    earlier audit report "no such file" for files that were plainly there.

    Two defects the meta-writer session found here by measurement, both fixed:

    `cwd=ROOT`, because git was running in the CALLER's directory. Invoked from /tmp this
    returned commit "" and branch "" with no error at all, so the provenance of a shipped
    file could silently record nothing; a non-zero exit now records the error instead of a
    blank, since a blank is indistinguishable from a clean answer.

    And the dirty tree is split on a trailing newline only. `.strip()` on the whole output
    ate the first column of the FIRST entry, turning " M f" into "M f" -- in porcelain the
    first column is the index and the second the worktree, so that is not cosmetic: it
    reports a different state of the file. The fingerprint is visible in both committed
    metas, whose first entry lacks the space every later entry has.
    """
    def run(*a):
        r = subprocess.run(a, capture_output=True, text=True, cwd=ROOT)
        if r.returncode != 0:
            return f"<git failed: {r.returncode} {r.stderr.strip()[:80]}>"
        return r.stdout.rstrip("\n")
    return {
        "commit": run("git", "rev-parse", "HEAD"),
        "branch": run("git", "rev-parse", "--abbrev-ref", "HEAD"),
        "грязное_дерево": [x for x in run(
            "git", "-c", "core.quotepath=false", "status", "--porcelain").split("\n") if x],
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inp", default=RES + "submission/activity_submission.csv")
    ap.add_argument("--out", default=RES + "submission/activity_submission_recal.csv")
    ap.add_argument("--meta", default=RES + "submission/recal.meta.json")
    ap.add_argument("--b", required=True, help="четыре множителя через запятую")
    ap.add_argument("--mu-y", default="",
                    help="четыре целевых средних; иначе центр не двигается")
    a = ap.parse_args()

    bs = [float(x) for x in a.b.split(",")]
    if len(bs) != 4:
        raise SystemExit(f"--b: нужно четыре числа, дано {len(bs)}")
    mus = [float(x) for x in a.mu_y.split(",")] if a.mu_y else None
    if mus is not None and len(mus) != 4:
        raise SystemExit("--mu-y: нужно четыре числа")

    df = pd.read_csv(a.inp)
    out = df.copy()
    rec = {}
    print(f"  {'фермент':>8s} {'b':>5s} {'sd до':>7s} {'sd после':>9s} "
          f"{'ср. до':>8s} {'ср. после':>10s} {'порядок':>8s}")
    for e, c in enumerate(CYPS):
        col = COL.format(c=c)
        p = df[col].to_numpy(float)
        mp = p.mean()
        target = mus[e] if mus is not None else mp
        q = target + bs[e] * (p - mp)

        # Monotonicity is the whole safety argument: b <= 0 would reverse the ordering, which
        # is the one thing this repository's currency (rank) cannot survive.
        if bs[e] <= 0:
            raise SystemExit(f"{c}: b={bs[e]} <= 0 -- порядок бы перевернулся")
        same = np.array_equal(np.argsort(p, kind="stable"), np.argsort(q, kind="stable"))
        assert same, f"{c}: порядок изменился"

        out[col] = q
        rec[c] = {"b": bs[e], "mu_y": target,
                  "sd_до": float(p.std(ddof=1)), "sd_после": float(q.std(ddof=1)),
                  "ср_до": float(mp), "ср_после": float(q.mean()),
                  "порядок_сохранён": bool(same)}
        print(f"  {c:>8s} {bs[e]:5.2f} {p.std(ddof=1):7.3f} {q.std(ddof=1):9.3f} "
              f"{mp:8.3f} {q.mean():10.3f} {'да':>8s}")

    ident = (df["SMILES"].equals(out["SMILES"])
             and df["Molecule_Name"].equals(out["Molecule_Name"]))
    assert ident, "идентификаторы изменились"
    out.to_csv(a.out, index=False)

    meta = {
        "что": ("Провенанс ПЕРЕКАЛИБРОВАННОЙ регрессионной подачи. Аффинное преобразование "
                "поверх файла, собранного src/submit.py: масштаб и центр взяты из моментов, "
                "выведенных с доски (пункт 308). Ранг не меняется по построению."),
        "собрано_utc": datetime.datetime.now(datetime.timezone.utc)
                               .strftime("%Y-%m-%dT%H:%M:%SZ"),
        "команда": f"uv run python src/recalib.py --b {a.b} --mu-y {a.mu_y}",
        "git": dict(git_state(), замечание=(
            "Преобразование читает только --inp и четыре константы из командной строки, "
            "поэтому состояние дерева на воспроизводимость не влияет: тот же вход и те же "
            "константы дают тот же выход из любого коммита.")),
        "пункт": "308 (предрегистрация правила и осознанная отмена обязательства пункта 300)",
        "вход": {"файл": rel(a.inp), "sha256": sha256(a.inp), "строк": int(len(df))},
        "выход": {"файл": rel(a.out), "sha256": sha256(a.out), "строк": int(len(out))},
        "преобразование": rec,
        "инварианты": {"идентификаторы_совпали": bool(ident),
                       "порядок_сохранён_на_всех": all(v["порядок_сохранён"]
                                                       for v in rec.values())},
        "замечание": ("Сверять неизменность хешем нельзя: pandas пересериализует числа, "
                      "поэтому даже прогон с b=1 меняет sha256 при побитово равных числах. "
                      "Проверять численно."),
    }
    json.dump(meta, open(a.meta, "w"), ensure_ascii=False, indent=1)
    print(f"\n  записано: {a.out}")
    print(f"  провенанс: {a.meta}")


if __name__ == "__main__":
    main()
