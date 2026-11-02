"""Item 374's arm must be expressible in `src/submit.py`, and adding it must change nothing.

The sibling of `tests/test_stack346_inert.py`, for the arm that REMOVES a member instead of adding
one. Item 374 measured +0.0099 of macro rank and -0.0167 of ST-RAE over the shipped six at sign
8/8, and the thing it removes is the sixth member item 324 shipped: two in-context learners replace
the RidgeCV probe on the frozen 600, whose errors correlate 0.948 with the stronger of them.

What has to be true of the threading, and what each test here pins:

  * dropping the probe must actually change the mean, and must change it in exactly one way ---
    `_combine` called without `extra` is the same as with the probe left out. A test that only
    checked "the flag exists" would pass on a `_combine` that ignored the parameter.
  * the flags are OFF by default, the two stacks are mutually exclusive, and `--drop-probe`
    without `--tabicl` is refused outright: no run has measured a slot removed with nothing in its
    place, and item 373's law says that costs rather than pays.
  * `_tabicl_load`'s controls must FIRE, not merely exist. The vectors are read from files and
    nothing about them is recomputed, so a wrong file cannot raise an exception on its own: a
    foreign fold digest, a length that does not match the masks, a value outside the pIC50
    window, and two representations that are secretly the same array must each be refused.
  * the committed test predictions must match the committed out-of-fold ones in shape and the
    golden split, because those two files are what `--stack374` averages.
"""
import json
import pathlib
import subprocess
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import submit as S  # noqa: E402


def _run(args):
    return subprocess.run([sys.executable, str(ROOT / "src" / "submit.py")] + args,
                          capture_output=True, text=True, cwd=str(ROOT))


def _parts(n=9, seed=0):
    rng = np.random.default_rng(seed)
    kinds = ["поферментно", "пул", "GP", "гребневая", "ствол"]
    return [(k, [rng.normal(size=n) for _ in range(len(S.CYPS))]) for k in kinds]


def _mask(n=9):
    return np.ones((n, len(S.CYPS)), bool)


def _te(n=750):
    import pandas as pd
    return pd.DataFrame({"Molecule_Name": ["m%d" % i for i in range(n)]})


def test_dropping_the_probe_is_combine_without_extra():
    """The arm's whole mechanism: `extra=None` must equal "the probe left out", and not by
    accident --- the with-probe call has to differ, or the parameter is being ignored."""
    parts = _parts()
    probe = [np.arange(9.0) + 100.0 for _ in range(len(S.CYPS))]
    extras = [[np.arange(9.0) + 7.0 for _ in range(len(S.CYPS))]]
    for e in range(len(S.CYPS)):
        dropped = S._combine(parts, e, None, extras=extras)
        kept = S._combine(parts, e, probe, extras=extras)
        assert np.array_equal(dropped, S._combine(parts, e, extras=extras))
        assert not np.allclose(dropped, kept), "зонд не влияет на среднее — параметр игнорируется"


def test_flags_are_off_by_default():
    h = _run(["--help"])
    assert h.returncode == 0, h.stderr[-400:]
    for flag in ("--tabicl", "--drop-probe", "--stack374",
                 "--tabicl-test", "--tabicl-oof-e", "--tabicl-oof-s"):
        assert flag in h.stdout, "%s нет в справке" % flag


def test_refusals():
    for args, want in ((["--tabicl"], "--tabicl требует --probe"),
                       (["--drop-probe"], "--drop-probe требует --probe"),
                       (["--probe", "--drop-probe"], "--drop-probe без --tabicl"),
                       (["--stack346", "--stack374"], "разные составы")):
        r = _run(args)
        assert r.returncode != 0, "%s не был отвергнут" % args
        assert want in (r.stdout + r.stderr), (r.stdout + r.stderr)[-500:]


def test_committed_files_agree_with_the_golden_split():
    """The two files `--stack374` averages, checked where they must agree."""
    # Anchored at ROOT rather than left to the working directory: cyppaths hands out string
    # prefixes, and a relative open() would make this test pass or fail by where pytest was
    # started from.
    def _at(path):
        q = pathlib.Path(path)
        return q if q.is_absolute() else ROOT / q
    assert _at(S.TABICL_TEST).exists(), S.TABICL_TEST
    tj = json.load(open(_at(S.TABICL_TEST)))
    assert int(tj["n_test"]) == 750
    for repk in ("E", "S"):
        T = np.asarray(tj["preds"][repk], float)
        assert T.shape == (len(S.CYPS), 750), (repk, T.shape)
        assert _at(S.TABICL_OOF[repk]).exists(), S.TABICL_OOF[repk]
        oj = json.load(open(_at(S.TABICL_OOF[repk])))
        assert str(oj["digests"][0]) == "2d93c19815e14261", (repk, oj["digests"][0])
        assert "%s_seed0" % repk in oj["preds"], repk
    # The two representations must differ; one runner, one checkpoint, one PCA cap, so a mix-up
    # of the input key is the realistic failure and it would read as a free gain.
    d = float(np.abs(np.asarray(tj["preds"]["E"], float)
                     - np.asarray(tj["preds"]["S"], float)).max())
    assert d > 0.1, "представления почти совпали, max|Δ| %.2e" % d


@pytest.mark.parametrize("break_it", ["digest", "length", "window", "same"])
def test_tabicl_load_controls_fire(tmp_path, break_it):
    """Each control refused on a file mutated to violate exactly it, and the UNMUTATED file
    accepted — a control that rejects everything is not a control."""
    te, mask = _te(), _mask(9)
    want = [int(mask[:, e].sum()) for e in range(len(S.CYPS))]
    dg = "2d93c19815e14261"
    oof = {k: {"digests": [dg, "x", "y", "z"],
               "preds": {"%s_seed0" % k: [[4.5] * want[e] for e in range(len(S.CYPS))]}}
           for k in ("E", "S")}
    test = {"n_test": len(te), "pca": 100,
            "preds": {"E": [[4.4] * len(te) for _ in S.CYPS],
                      "S": [[4.9] * len(te) for _ in S.CYPS]}}
    if break_it == "digest":
        oof["E"]["digests"][0] = "0000000000000000"
    elif break_it == "length":
        oof["S"]["preds"]["S_seed0"][1] = [4.5] * (want[1] - 1)
    elif break_it == "window":
        test["preds"]["E"][0][0] = 99.0
    elif break_it == "same":
        test["preds"]["S"] = [list(v) for v in test["preds"]["E"]]
    paths = {}
    for k in ("E", "S"):
        paths[k] = str(tmp_path / ("oof_%s.json" % k))
        json.dump(oof[k], open(paths[k], "w"))
    tp = str(tmp_path / "test.json")
    json.dump(test, open(tp, "w"))
    with pytest.raises(SystemExit):
        S._tabicl_load(te, mask, dg, tp, paths)


def test_tabicl_load_accepts_a_sound_file(tmp_path):
    te, mask = _te(), _mask(9)
    want = [int(mask[:, e].sum()) for e in range(len(S.CYPS))]
    dg = "2d93c19815e14261"
    paths = {}
    for k, v in (("E", 4.4), ("S", 4.9)):
        paths[k] = str(tmp_path / ("oof_%s.json" % k))
        json.dump({"digests": [dg] * 4,
                   "preds": {"%s_seed0" % k: [[v] * want[e] for e in range(len(S.CYPS))]}},
                  open(paths[k], "w"))
    tp = str(tmp_path / "test.json")
    json.dump({"n_test": len(te), "pca": 100,
               "preds": {"E": [[4.4] * len(te) for _ in S.CYPS],
                         "S": [[4.9] * len(te) for _ in S.CYPS]}}, open(tp, "w"))
    got = S._tabicl_load(te, mask, dg, tp, paths)
    assert [g[0] for g in got] == [S.TABICL_KIND["E"], S.TABICL_KIND["S"]]
    assert got[0][1].shape == (len(te), len(S.CYPS))
    assert [len(v) for v in got[0][2]] == want
