"""Guards on item 326: the dead-zone targets are clipped from members already computed.

`main()` used to run `oof_members` TWICE under default flags -- once inside the shrinkage
pass and once more to build the dead-zone targets, on folds drawn by a second
`butina_folds` call over the same SMILES. Those folds are equal element for element, and
every member is deterministic, so the second pass reproduced the first bit for bit at a
cost of 943 s, about a ninth of the run.

Two things are pinned here, and the second exists because the first is not enough. The
TRANSPORT: that `raw_out` hands back the members as they stood BEFORE `dz_pass`, with the
aliasing behaviour the real `dz_pass` actually has rather than a convenient one. And the
PAYLOAD: that `main()` in fact reuses them. Without the second, reverting the whole change
leaves every other test in this repository passing -- which was true of the first version
of this file and is the defect it was rewritten to close.

The arithmetic is pinned by the run instead: `src/submit.py --probe` must still reproduce
`results/submission/activity_submission.csv` element for element through its
`activity_submission_noprobe.csv` twin. Everything in this file is stubbed or static, so
it costs milliseconds and needs no data.
"""
import ast
import pathlib
import sys

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

SUBMIT = ROOT / "src" / "submit.py"
MEMBERS = ("поферментно", "пул", "GP", "гребневая", "ствол")


@pytest.fixture()
def S():
    import submit
    saved = {n: getattr(submit, n) for n in ("oof_members", "dz_pass")}
    yield submit
    for n, f in saved.items():
        setattr(submit, n, f)


def _members():
    """Five members shaped like the real ones. Every array is DISTINCT, so an assertion
    about order or identity is able to fail; the first version of this file gave every
    member the same values and could not."""
    return [(k, [np.full(3, 10.0 * i + j) for j in range(4)])
            for i, k in enumerate(MEMBERS)]


def _dz_stub(parts, *a, **k):
    """dz_pass's own behaviour, including the part that is easy to stub away: it allocates
    for every member it passes through, and appends "ствол" BY REFERENCE. A stub that
    allocated for the trunk too would give the code under test a property it does not
    have, and the aliasing test below could not fail."""
    out = []
    for kind, P in parts:
        if kind == "ствол":
            out.append((kind, P))
            continue
        out.append((kind, [p + 100.0 for p in P]))
    return out


def _install(S, counter):
    S.oof_members = lambda *a, **k: (counter.append(1), _members())[1]
    S.dz_pass = _dz_stub


# --------------------------------------------------------------------- transport

def test_raw_out_receives_the_pre_dead_zone_members(S):
    calls = []
    _install(S, calls)
    raw = []
    out = S.oof_parts(None, None, None, None, "ансамбль5", bands=("lo", "hi"), raw_out=raw)

    assert len(raw) == 1, "ожидался ровно один проход вне фолда"
    assert calls == [1], "зеркалирование членов не должно стоить лишнего прохода"
    assert [k for k, _ in raw[0]] == list(MEMBERS)
    for (kind, A), (_, B) in zip(raw[0], out):
        for a_, b_ in zip(A, B):
            if kind == "ствол":
                continue
            assert not np.array_equal(a_, b_)
            assert a_.max() < 100.0, "raw_out отдал предсказания ПОСЛЕ прохода мёртвой зоны"


def test_non_trunk_members_are_not_aliased(S):
    """Writing through the returned parts must not reach what main() will clip."""
    calls = []
    _install(S, calls)
    raw = []
    out = S.oof_parts(None, None, None, None, "ансамбль5", bands=("lo", "hi"), raw_out=raw)
    before = {k: [a_.copy() for a_ in A] for k, A in raw[0] if k != "ствол"}
    for kind, P in out:
        if kind == "ствол":
            continue
        for p in P:
            p += 1000.0
    after = {k: A for k, A in raw[0] if k != "ствол"}
    for k in before:
        for b_, a_ in zip(before[k], after[k]):
            assert np.array_equal(b_, a_), f"член {k} записался сквозь raw_out"


def test_trunk_member_IS_shared_and_that_is_why_main_skips_it(S):
    """The honest half. dz_pass appends the trunk by reference, so raw_out and the return
    value share those arrays; the docstring of oof_parts says so. It is harmless only
    because main()'s target loop skips "ствол" -- pinned separately, statically, below."""
    calls = []
    _install(S, calls)
    raw = []
    out = S.oof_parts(None, None, None, None, "ансамбль5", bands=("lo", "hi"), raw_out=raw)
    r = dict(raw[0])["ствол"]
    o = dict(out)["ствол"]
    assert all(x is y for x, y in zip(r, o)), (
        "заглушка dz_pass перестала повторять настоящую: ствол должен идти по ССЫЛКЕ")


def test_oof_predictions_forwards_raw_out(S):
    """The branch the verifying --probe run does not exercise: plain default goes here."""
    calls = []
    _install(S, calls)
    raw = []
    S.oof_predictions(None, None, None, None, "ансамбль5", bands=("lo", "hi"), raw_out=raw)
    assert len(raw) == 1, "oof_predictions не передала raw_out дальше"
    assert calls == [1]


def test_default_call_is_unchanged(S):
    """No raw_out: the old behaviour, with bands and without."""
    calls = []
    _install(S, calls)
    out = S.oof_parts(None, None, None, None, "ансамбль5", bands=("lo", "hi"))
    assert calls == [1]
    assert any(p.max() >= 100.0 for k, P in out if k != "ствол" for p in P), \
        "проход мёртвой зоны не отработал"

    calls.clear()
    plain = S.oof_parts(None, None, None, None, "ансамбль5")
    assert calls == [1]
    assert all(p.max() < 100.0 for _, P in plain for p in P), \
        "проход мёртвой зоны отработал там, где полос нет"


# ----------------------------------------------------------------------- payload
#
# main() is one long function with no seam a unit test can reach, and the reuse it performs
# is exactly the thing a well-meaning revert would undo. These read the source instead.

def _main_def():
    tree = ast.parse(SUBMIT.read_text(encoding="utf-8"))
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name == "main":
            return node
    raise AssertionError("в src/submit.py нет функции main --- тест ослеп, а не прошёл")


def _calls_to(node, name):
    return [n for n in ast.walk(node)
            if isinstance(n, ast.Call) and isinstance(n.func, ast.Name) and n.func.id == name]


def _fallback_guarded_calls(main):
    """Direct oof_members calls that sit in the `else` of a test on `raw_members`."""
    ids = set()
    for node in ast.walk(main):
        if not isinstance(node, ast.If):
            continue
        if "raw_members" not in {n.id for n in ast.walk(node.test)
                                 if isinstance(n, ast.Name)}:
            continue
        for b in node.orelse:
            ids.update(id(c) for c in _calls_to(b, "oof_members"))
    return ids


def test_main_never_computes_the_members_on_the_main_path():
    """The payload of item 326, and the only assertion in this repository that a revert
    cannot satisfy.

    Counting the calls is NOT enough and the first version of this test made exactly that
    mistake: before the change main() also contained a single DIRECT oof_members call --
    the second pass -- because the first pass went through oof_parts. The count was one on
    both sides of the change. What distinguishes them is WHERE the surviving call sits: it
    must be reachable only when raw_members is empty."""
    main = _main_def()
    direct = _calls_to(main, "oof_members")
    guarded = _fallback_guarded_calls(main)
    assert guarded, ("в main() нет вызова oof_members в запасной ветке по raw_members --- "
                     "переиспользование пропало или ветвление переписано")
    unguarded = [c for c in direct if id(c) not in guarded]
    assert not unguarded, (
        f"{len(unguarded)} из {len(direct)} вызовов oof_members в main() стоят на основном "
        f"пути; после пункта 326 там не должно быть ни одного --- члены берутся из raw_out")
    assert len(direct) == 1, (
        f"main() зовёт oof_members {len(direct)} раз, ожидался один запасной проход")


def test_the_trunk_is_still_skipped_when_targets_are_built():
    """dz_pass leaves the trunk un-reprojected (item 164), and the target loop has always
    skipped it. The reuse hands main() an array the return value SHARES, so clipping it
    here would write through -- this is the assertion that keeps that unreachable."""
    main = _main_def()
    skips = [n for n in ast.walk(main)
             if isinstance(n, ast.If)
             and any(isinstance(c, ast.Constant) and c.value == "ствол"
                     for c in ast.walk(n.test))
             and any(isinstance(b, ast.Continue) for b in n.body)]
    assert skips, 'в main() пропал `if kind == "ствол": continue` при сборке мишеней'
