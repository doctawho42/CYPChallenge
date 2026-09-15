"""Guards on the provenance record for the submitted files.

Item 309 found that results/submission/submission.meta.json had no writer anywhere in the
repository: the file says of itself that it is assembled programmatically, and nothing could
assemble it, so a rerun produced the expensive artefact and no account of it. `src/submeta.py`
is that writer, and these checks hold it to recording what the run did.

The committed record is the specification. `build` must produce every key it has except the
three no run can derive, which are named in UNDERIVABLE so that dropping any other key fails
here rather than quietly shrinking the record.
"""
import json
import pathlib
import sys

import pytest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / "src"))

import submeta  # noqa: E402

COMMITTED = (pathlib.Path(__file__).resolve().parents[1]
             / "results" / "submission" / "submission.meta.json")

# Two prose judgements about one tree and one commit, and a manual archiving step. A run can
# report the dirty tree and the verdicts; it cannot decide that the dirty files do not matter.
UNDERIVABLE = {"top": ["архив_предыдущего"],
               "git": ["замечание"],
               "гейт_валидатора": ["замечание"]}

# sha256 of FIXTURE_CSV, computed with shasum rather than by the code under test.
FIXTURE_SHA = "92abeff4cc56d66835cb80f2c437a57984fdc6fb8f189f2d0a56fd2754dd2309"
FIXTURE_CSV = "SMILES,Molecule_Name\nC,a\nCC,b\nCCC,c\n"

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
SOLO = {"CYP1A2": ("поферментно", "GP", "ствол"),
        "CYP2C9": ("GP", "ствол"),
        "CYP3A4": ("GP", "ствол")}


@pytest.fixture
def context(tmp_path):
    """What the submitting run knows at the moment the gate accepts its files."""
    act = tmp_path / "activity_submission.csv"
    act.write_text(FIXTURE_CSV, encoding="utf-8")
    return dict(
        paths=[act],
        argv=[],
        mode="ансамбль5",
        solo=SOLO,
        cyps=CYPS,
        deadzone=True,
        bundle=True,
        delta="0.1,0.3,0,0.7",
        # (lambda, centre, shift) per enzyme, exactly as fit_shrinkage returns them. The
        # centre is a fitting internal and must not reach the record.
        lams=[(1.08, 4.21, 0.13), (1.22, 4.80, 0.17), (1.1, 3.07, 0.01), (1.04, 4.81, 0.24)],
        grid=(0.2, 2.0),
        positives={"CYP3A4": 360, "CYP2D6": 285},
        total_rows=750,
        digest="2d93c19815e14261",
        clusters=4703,
        golden="2d93c19815e14261",
        gate={"регрессия": True, "классификация": True},
        git={"commit": "0" * 40, "branch": "main", "грязное_дерево": []},
        versions={"sklearn": "1.3.2", "numpy": "1.26.4", "torch": "2.13.0", "python": "3.12.5"},
        elapsed_s=9725.4,
    )


def test_it_records_every_section_the_committed_record_has(context):
    """The committed record is the specification, and a section lost in a refactor is a
    provenance hole rather than a test failure anyone would otherwise notice."""
    meta = submeta.build(**context)
    committed = json.loads(COMMITTED.read_text(encoding="utf-8"))
    assert [k for k in committed if k not in meta] == UNDERIVABLE["top"]
    for section, expected in [("git", UNDERIVABLE["git"]), ("состав", []),
                              ("преобразование", []), ("классификация", []),
                              ("разбиение", []), ("версии", []),
                              ("гейт_валидатора", UNDERIVABLE["гейт_валидатора"])]:
        assert [k for k in committed[section] if k not in meta[section]] == expected, section
    assert meta["собрано_utc"].endswith("+00:00")


def test_it_hashes_and_counts_the_files_it_names(context):
    meta = submeta.build(**context)
    (entry,) = [v for k, v in meta["файлы"].items() if k.endswith("activity_submission.csv")]
    assert entry == {"sha256": FIXTURE_SHA, "строк": 3}


def test_the_affine_pair_is_recorded_per_enzyme_without_the_fitting_centre(context):
    """`fit_shrinkage` returns (lambda, centre, shift). The centre is the mean of the
    out-of-fold predictions, not a property of the submission, and the committed record has
    no place for it: a writer that dumps the tuple would put it there."""
    moved = submeta.build(**context)["преобразование"]
    assert moved["lambda_по_ферментам"] == {"CYP1A2": 1.08, "CYP2C9": 1.22,
                                            "CYP2D6": 1.1, "CYP3A4": 1.04}
    assert moved["сдвиг_предсказаний"] == {"CYP1A2": 0.13, "CYP2C9": 0.17,
                                           "CYP2D6": 0.01, "CYP3A4": 0.24}
    assert "4.21" not in json.dumps(moved, ensure_ascii=False)
    assert moved["GRID"] == "0.2..2.0"


def test_the_record_carries_how_long_the_run_took(context):
    """A documented runtime is what tells the next reader whether they dare rerun something,
    and src/submit.py takes 161 to 230 minutes. 9725.4 seconds is 162.09 minutes."""
    took = submeta.build(**{**context, "elapsed_s": 9725.4})["время_прогона"]
    assert took["секунд"] == 9725
    assert took["минут"] == 162.1


def test_an_enzyme_outside_solo_is_named_rather_than_left_out(context):
    """CYP2D6 has no SOLO entry because it ships every member, which is a fact about the
    submission and reads as an omission if the record simply lacks the key."""
    assert submeta.build(**context)["состав"]["CYP2D6"] == "все пять членов (в SOLO не входит)"


def test_the_positive_calls_are_recorded_with_the_total_they_are_out_of(context):
    assert submeta.build(**context)["классификация"]["положительных"] == {
        "CYP3A4": 360, "CYP2D6": 285, "из": 750}


def test_a_moved_split_is_recorded_as_a_mismatch(context):
    """The golden digest is the repository's central agreement. A run on a moved split must
    say so in the record rather than carry a digest nobody compares."""
    meta = submeta.build(**{**context, "digest": "0000000000000000"})
    assert meta["разбиение"]["digest"] == "0000000000000000"
    assert meta["разбиение"]["совпало"] is False
    assert submeta.build(**context)["разбиение"]["совпало"] is True


def test_the_gate_verdicts_are_the_ones_the_validators_gave(context):
    meta = submeta.build(**{**context, "gate": {"регрессия": True, "классификация": False}})
    assert meta["гейт_валидатора"] == {"регрессия": "принято", "классификация": "ОТКЛОНЕНО"}


def test_the_command_records_the_flags_the_run_was_given(context):
    """Provenance that cannot distinguish a default run from a flagged one cannot answer the
    question it exists for: which configuration produced the file on disk."""
    assert submeta.build(**context)["команда"] == ("uv run python src/submit.py"
                                                   "   (все флаги по умолчанию)")
    flagged = submeta.build(**{**context, "argv": ["--mode", "ансамбль", "--no-deadzone"]})
    assert flagged["команда"] == "uv run python src/submit.py --mode ансамбль --no-deadzone"


def test_the_file_keeps_the_committed_record_s_formatting(context, tmp_path):
    """The committed record round-trips through json.dumps(ensure_ascii=False, indent=1)
    byte for byte, and a writer that escapes the Cyrillic keys makes every future diff of
    this file unreadable."""
    out = tmp_path / "written.meta.json"
    submeta.write(out, **context)
    text = out.read_text(encoding="utf-8")
    assert '"что"' in text
    assert text == json.dumps(json.loads(text), ensure_ascii=False, indent=1)


def test_a_run_without_shrinkage_records_that_no_pair_was_applied(context):
    """`--no-shrink` leaves `lams` as None in src/submit.py. The record has to say so,
    not crash after both files are already on disk and not invent a pair."""
    moved = submeta.build(**{**context, "lams": None})["преобразование"]
    assert moved["lambda_по_ферментам"] is None
    assert moved["сдвиг_предсказаний"] is None
    assert "--no-shrink" in moved["пункт"]


def test_a_lambda_on_the_grid_edge_is_named_rather_than_reported_as_interior(context):
    """Item 277 widened the grid because optima sat on its edge. The committed record's
    "край не задет" is a claim about one run, so the writer has to check it on each run."""
    interior = submeta.build(**context)["преобразование"]["пункт"]
    assert interior == "277 (сетка расширена выше 1.0; все четыре оптимума внутри, край не задет)"
    lams = [(1.08, 4.21, 0.13), (1.22, 4.80, 0.17), (2.0, 3.07, 0.01), (1.04, 4.81, 0.24)]
    edge = submeta.build(**{**context, "lams": lams})["преобразование"]["пункт"]
    assert "CYP2D6" in edge and "внутри" not in edge


def test_a_non_default_run_describes_the_arms_it_actually_ran(context):
    """The item numbers and arm descriptions in the committed record are true of the default
    run only. Copied as constants, they would describe a --no-bundle run as the bundle."""
    meta = submeta.build(**{**context, "mode": "ансамбль", "deadzone": False, "bundle": False})
    assert "пять" not in meta["состав"]["CYP2D6"]
    assert "120" not in meta["состав"]["пункты"] and "164" not in meta["состав"]["пункты"]
    assert "Платт" in meta["классификация"]["плечо"]
    assert "244" not in meta["классификация"]["пункты"]


def test_solo_without_the_dead_zone_is_recorded_as_reaching_only_the_trunk(context):
    """In src/submit.py the dead-zone loop sends every member through `_keep`. Without it the
    test path filters only the trunk, while `oof_predictions` still fits the affine pair on
    the SOLO-filtered average. Recording SOLO as such a run's composition would be false."""
    assert "SOLO_на_тесте" not in submeta.build(**context)["состав"]
    flagged = submeta.build(**{**context, "deadzone": False})["состав"]
    assert "ствол" in flagged["SOLO_на_тесте"]


def test_git_state_describes_the_repository_whatever_the_working_directory(monkeypatch,
                                                                           tmp_path):
    """recalib.git_state runs git in the caller's working directory, and outside a checkout
    it returns an empty commit and a clean tree without a word (measured, not supposed). A
    record that says "no commit, nothing modified" when it could not look is worse than none."""
    monkeypatch.setenv("GIT_CEILING_DIRECTORIES", str(tmp_path.parent))
    monkeypatch.chdir(tmp_path)
    state = submeta.git_state()
    assert len(state["commit"]) == 40 and state["branch"]
    blind = submeta.git_state(tmp_path)
    assert blind["commit"] is None and blind["ошибка"]


def test_git_state_keeps_the_status_column_of_the_first_dirty_entry(monkeypatch, tmp_path):
    """Stripping git's whole status output eats the leading space of the first entry only,
    so " M" (modified, not staged) reads as "M" (staged). The hand-made record shows exactly
    that: its first entry is "M verify/README.md" while its second keeps " M"."""
    import subprocess
    repo = tmp_path / "repo"
    repo.mkdir()

    def git(*args):
        subprocess.run(["git", "-C", str(repo), "-c", "user.name=t", "-c", "user.email=t@t",
                        "-c", "commit.gpgsign=false", "-c", "core.hooksPath=/dev/null", *args],
                       check=True, capture_output=True)
    git("init", "-q")
    (repo / "a.txt").write_text("1", encoding="utf-8")
    git("add", "a.txt")
    git("commit", "-q", "-m", "one")
    (repo / "a.txt").write_text("2", encoding="utf-8")
    monkeypatch.setenv("GIT_CEILING_DIRECTORIES", str(tmp_path))
    assert submeta.git_state(repo)["грязное_дерево"] == [" M a.txt"]


def test_versions_are_read_from_the_running_interpreter():
    import platform

    import numpy
    found = submeta.versions()
    assert set(found) == {"sklearn", "numpy", "torch", "python"}
    assert found["python"] == platform.python_version()
    assert found["numpy"] == numpy.__version__


def _import_submit():
    """src/submit.py puts the organisers' repository on sys.path when it is imported, which
    raises where the submodule is not checked out. Say so rather than error."""
    try:
        import submit
    except FileNotFoundError as err:
        pytest.skip(f"src/submit.py needs the organisers' repository: {err}")
    return submit


def test_submit_hands_the_writer_what_the_run_actually_did(context):
    """The call sits at the end of a 161-to-230-minute run, after both files are written and
    accepted, so a slip in it would surface three hours in. This exercises it in milliseconds."""
    import argparse

    import numpy as np
    import pandas as pd
    S = _import_submit()
    args = argparse.Namespace(mode="ансамбль5", deadzone=True, bundle=True,
                              delta="0.1,0.3,0,0.7")
    cls = pd.DataFrame({"CYP3A4_is_TDI": [True, False, True],
                        "CYP2D6_is_TDI": [False, False, True]})
    ctx = S._provenance(args, [], context["lams"], np.array([0, 1, 2, 0, 1]), 3, cls,
                        {"регрессия": True, "классификация": True}, context["paths"],
                        elapsed_s=5400.0)
    meta = submeta.build(**ctx)
    assert meta["время_прогона"]["минут"] == 90.0
    assert meta["классификация"]["положительных"] == {"CYP3A4": 2, "CYP2D6": 1, "из": 3}
    assert meta["преобразование"]["GRID"] == "0.2..2.0"
    assert set(meta["состав"]["SOLO"]) == set(S.SOLO)
    assert meta["разбиение"]["кластеров"] == 3 and meta["разбиение"]["совпало"] is False
    assert len(meta["git"]["commit"]) == 40
