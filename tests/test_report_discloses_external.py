"""The outward-facing documents must disclose the third-party encoder the submission uses.

Why this test exists. On 26 September a ridge head on a third-party pretrained checkpoint shipped as
the sixth member of the regression ensemble (items 324, 325, 327). In the same commit,
`docs/method_report.md` kept saying "External data: the organisers' own files only" and `METHOD.md`
kept saying the pretrained representations were "all closed negative in the journal". Both stood for
three weeks on a public repository whose report field the leaderboard links.

Nothing could have caught it. `tests/test_journal_citations.py` guards the journal's citations and
`tests/test_split.py` guards the fold digest, but no test read the documents, so a disclosure could
drift away from the code indefinitely and silently. This closes that, loudly: if the shipped path
still carries a third-party encoder, both documents must name where it came from.

The premise is asserted rather than assumed. If `src/submit.py` ever stops carrying the encoder the
test's own reason disappears, and it says so instead of passing vacuously.
"""

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
# The third-party repository the checkpoint is committed in. Not a URL: the string that must appear.
SOURCE = "jeremycheminf/openadmet_scripts"
OUTWARD = ["METHOD.md", "docs/method_report.md"]


def _read(rel):
    p = ROOT / rel
    assert p.exists(), f"{rel} отсутствует, а тест написан про него"
    return p.read_text(encoding="utf-8")


def test_the_shipped_path_still_carries_a_third_party_encoder():
    """The premise. Without it this test would pass for the wrong reason."""
    src = _read("src/submit.py")
    assert "PROBE_EMB" in src, (
        "src/submit.py больше не содержит PROBE_EMB: сторонний энкодер, видимо, убран из "
        "отгружаемого пути. Тогда этот тест потерял смысл и его надо переписать или удалить, "
        "а не чинить документы."
    )
    assert "emb_chemprop_medium" in src, (
        "ожидался тот самый чекпойнт в отгружаемом пути; если он сменился, обновите SOURCE "
        "и документы вместе с ним"
    )


def test_both_outward_documents_name_the_encoder_source():
    missing = [f for f in OUTWARD if SOURCE not in _read(f)]
    assert not missing, (
        "в %s не назван источник стороннего чекпойнта (%s). Отгружаемая модель использует внешние "
        "данные, и публичные документы обязаны это раскрывать --- именно это и разошлось с кодом "
        "на три недели." % (", ".join(missing), SOURCE)
    )


def test_neither_document_still_claims_only_the_organisers_files():
    """The two exact sentences that were false, so they cannot come back unnoticed.

    Matched outside quotation marks only: both documents keep the old wording as a dated quote of
    their own error, and a test that could not tell a claim from a quotation of it would fire on the
    correction itself.
    """
    bad = {
        "docs/method_report.md": r"External data:\s*the organisers' own\s*\n?files only",
        "METHOD.md": r"all of them are closed negative in the\s*\n?journal",
    }
    for f, pat in bad.items():
        text = _read(f)
        unquoted = re.sub(r'"[^"]*"', "", text)      # drop quoted spans, errors are quoted there
        hit = re.search(pat, unquoted)
        assert hit is None, (
            "%s снова утверждает опровергнутое: %r. Если состав подачи изменился и утверждение "
            "снова верно, правьте тест вместе с документом." % (f, pat)
        )


def test_the_search_can_fail():
    """A control: the matcher must be capable of finding something in these files.

    Without it, a typo in a path or an encoding mishap would make every assertion above pass by
    finding nothing anywhere, which is the defect this repository calls a query that could not have
    succeeded.
    """
    for f in OUTWARD:
        text = _read(f)
        assert "ST-RAE" in text, f"{f} прочитан, но не содержит 'ST-RAE' --- читается не тот файл"
        assert "zzz-no-such-string-" + "sentinel" not in text, "контроль-пустышка сработал"
