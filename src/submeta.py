"""Provenance record for the submitted files, assembled by the run that writes them.

Item 309 found that results/submission/submission.meta.json had no writer anywhere in this
repository: it says of itself that it is assembled programmatically, and nothing could
assemble it. A rerun therefore produced the expensive artefact and no account of it.

`build` is a pure function of what the run knows, so the record can be tested without a
161-to-230-minute rerun; `write` is the thin wrapper src/submit.py calls once the
organisers' validators have accepted both files.

Three keys of the hand-written record are not reproduced, because no run can know them: the
git note that the dirty files do not touch submit.py's imports (a judgement), the gate's
historical remark about commit fde6b22, and the path of the previous submission's archive
(a manual step). tests/test_submission_meta.py names them.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))

import datetime
import importlib.metadata
import json
import platform
import subprocess

import pandas as pd

from recalib import rel, sha256

ROOT = _pl.Path(__file__).resolve().parents[1]


def git_state(root=ROOT):
    """Commit, branch and dirty tree of the repository at `root`, not of the working directory.

    recalib.git_state asks git in the caller's working directory, and outside a checkout it
    returns an empty commit and a clean tree without complaint. Here git is pointed at the
    repository, and a failure is recorded as a failure rather than as an empty answer. The
    status output is split without stripping it first: stripping the whole output eats the
    leading status column of the first entry only, which is visible in the hand-made record.
    """
    def run(*args):
        return subprocess.run(["git", "-C", str(root), "-c", "core.quotepath=false", *args],
                              capture_output=True, text=True)
    head = run("rev-parse", "HEAD")
    if head.returncode:
        return {"commit": None, "branch": None, "грязное_дерево": None,
                "ошибка": head.stderr.strip() or f"git вернул код {head.returncode}"}
    return {"commit": head.stdout.strip(),
            "branch": run("rev-parse", "--abbrev-ref", "HEAD").stdout.strip(),
            "грязное_дерево": [x for x in run("status", "--porcelain").stdout.split("\n") if x]}


def versions():
    """Versions of the interpreter doing the writing, read from installed metadata so that
    recording torch does not import it."""
    def installed(dist):
        try:
            return importlib.metadata.version(dist)
        except importlib.metadata.PackageNotFoundError:
            return None
    return {"sklearn": installed("scikit-learn"), "numpy": installed("numpy"),
            "torch": installed("torch"), "python": platform.python_version()}


WHAT = ("Провенанс ПОДАВАЕМЫХ файлов: каким кодом, в каком составе и на каком разбиении собран "
        "лежащий на диске сабмит. Пишется src/submit.py (через src/submeta.py) после того, как "
        "валидаторы организаторов приняли оба файла, а не переписывается руками.")


def _composition(mode, solo, cyps, deadzone):
    # Item references are true of the arms that ran, so they are assembled from the flags
    # rather than copied from the record of one particular run.
    items = ["282-285 (состав)"]
    if mode == "ансамбль5":
        items.append("120 (пятый член)")
    if deadzone:
        items.append("164/213 (мёртвая зона)")
    full = "все пять членов" if mode == "ансамбль5" else f"все члены режима {mode}"
    out = {"mode": mode,
           "SOLO": {c: list(members) for c, members in solo.items()},
           **{c: f"{full} (в SOLO не входит)" for c in cyps if c not in solo},
           "пункты": ", ".join(items)}
    if not deadzone:
        # Without the dead zone, submit.py's test path sends only the trunk through _keep,
        # while oof_predictions fits the affine pair on the SOLO-filtered average.
        reached = "только ствол" if mode == "ансамбль5" else "ни один член"
        out["SOLO_на_тесте"] = (f"{reached}: без мёртвой зоны остальные члены теста SOLO не "
                                f"фильтрует, хотя усадка подогнана на отфильтрованном составе")
    return out


def _transform(delta, lams, grid, cyps):
    lo, hi = grid
    if lams is None:
        return {"delta": delta, "lambda_по_ферментам": None, "сдвиг_предсказаний": None,
                "GRID": f"{lo}..{hi}", "пункт": "усадка не применялась (--no-shrink)"}
    # fit_shrinkage returns (lambda, centre, shift); the centre is a fitting internal.
    edge = [c for c, (lam, _, _) in zip(cyps, lams) if min(abs(lam - lo), abs(lam - hi)) < 1e-9]
    where = ("все четыре оптимума внутри, край не задет" if not edge
             else "оптимум на краю сетки: " + ", ".join(edge))
    return {"delta": delta,
            "lambda_по_ферментам": {c: lam for c, (lam, _, _) in zip(cyps, lams)},
            "сдвиг_предсказаний": {c: shift for c, (_, _, shift) in zip(cyps, lams)},
            "GRID": f"{lo}..{hi}",
            "пункт": f"277 (сетка расширена выше 1.0; {where})"}


def _classification(bundle, positives, total_rows):
    if bundle:
        arm = "связка: ворота(ансамбль+мёртвая зона) * сдвиг, калибровка, plug-in"
        items = "244 (предрегистрация), 250 (принято)"
    else:
        arm = "прямой классификатор + Платт, plug-in (--no-bundle)"
        items = "250 (связка отключена: поведение до этого пункта)"
    return {"плечо": arm, "пункты": items, "положительных": {**positives, "из": total_rows}}


def build(*, paths, argv, mode, solo, cyps, deadzone, bundle, delta, lams, grid, positives,
          total_rows, digest, clusters, golden, gate, git, versions, elapsed_s):
    return {
        "что": WHAT,
        "собрано_utc": (datetime.datetime.now(datetime.timezone.utc)
                        .replace(microsecond=0).isoformat()),
        # A documented runtime is what tells the next reader whether a rerun is affordable.
        "время_прогона": {"секунд": int(round(elapsed_s)), "минут": round(elapsed_s / 60, 1),
                          "от": "старта main() до записи этого файла"},
        "команда": "uv run python src/submit.py" + (" " + " ".join(argv) if argv
                                                    else "   (все флаги по умолчанию)"),
        "git": dict(git),
        "файлы": {rel(p): {"sha256": sha256(p), "строк": int(len(pd.read_csv(p)))}
                  for p in paths},
        "состав": _composition(mode, solo, cyps, deadzone),
        "преобразование": _transform(delta, lams, grid, cyps),
        "классификация": _classification(bundle, positives, total_rows),
        "разбиение": {"digest": digest, "кластеров": clusters, "золотое": golden,
                      "совпало": digest == golden},
        "версии": dict(versions),
        "гейт_валидатора": {name: "принято" if passed else "ОТКЛОНЕНО"
                            for name, passed in gate.items()},
    }


def write(path, **context):
    """Same serialisation as the committed record, which round-trips through it byte for
    byte: Cyrillic kept as written, one-space indent, no trailing newline."""
    meta = build(**context)
    _pl.Path(path).write_text(json.dumps(meta, ensure_ascii=False, indent=1), encoding="utf-8")
    return meta
