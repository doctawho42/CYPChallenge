"""Every journal item number cited anywhere in the repository must exist as a heading.

WHY THIS EXISTS. Item 355 was committed on 17 October with its script and its log, and its finding
was left in the commit message: the heading `**355.` was never written into `verify/README.md`. For
two days nothing failed, while `verify/k134_probealpha.py` and `verify/k135_probehead.py` each
carried a control that must "reproduce item 355's +0.0148", three passages of the journal cited item
355, and `docs/board_gap.md` built an action list on it.

`CLAUDE.md` already says to cite a symbol rather than a line number, because "a symbol either exists
or the grep for it fails loudly". An item NUMBER is the same hazard one level up: it is prose, it
resolves in the reader's head, and nothing checks it. This test is the loud failure.

THE CONTROL THIS TEST NEEDS ON ITSELF. A matcher that silently matches nothing would pass forever,
so `test_matcher_finds_real_citations` asserts the premise: the pattern must find known-good
citations in a file that certainly contains them. Verified to fail on the real defect --- run
against `verify/README.md` as of `6f0c80f` it reports four files citing the absent item 355, and
against the repaired journal it reports none.
"""

import pathlib
import re

ROOT = pathlib.Path(__file__).resolve().parents[1]
JOURNAL = ROOT / "verify" / "README.md"

# "item 324", "items 324/325", "пункт 346", "пункты 110 и 111", "пп. 244/250", "п. 77"
CITATION = re.compile(
    r"(?:items?|пункт(?:ы|а|е|ов)?|пп?\.)\s*(\d+(?:\s*(?:[/,]|and|и)\s*\d+)*)",
    re.IGNORECASE,
)
# Numbering starts at 22 and 1--21 were never used (CLAUDE.md), so anything below 22 is prose
# ("item 3 of the list") rather than a journal reference. The upper guard keeps a four-digit
# figure that happens to follow the word "item" from being read as a citation.
LOW, HIGH = 22, 400

GLOBS = ("verify/*.py", "src/*.py", "tests/*.py", "docs/*.md",
         "verify/README.md", "README.md", "CONTRIBUTING.md", "METHOD.md")


def _headings():
    text = JOURNAL.read_text(encoding="utf-8")
    return {int(m) for m in re.findall(r"(?m)^\*\*(\d+)\.[ \t]", text)}


def _cited(path):
    text = path.read_text(encoding="utf-8", errors="replace")
    out = set()
    for m in CITATION.finditer(text):
        out |= {int(x) for x in re.findall(r"\d+", m.group(1))}
    return {c for c in out if LOW <= c <= HIGH}


def _files():
    seen, out = set(), []
    for g in GLOBS:
        for p in sorted(ROOT.glob(g)):
            if p not in seen:
                seen.add(p)
                out.append(p)
    return out


def test_journal_has_headings():
    """Premise: the heading pattern must actually find the journal's items."""
    h = _headings()
    assert len(h) > 300, f"нашлось только {len(h)} заголовков — сломан шаблон, а не журнал"
    assert min(h) == LOW, f"нумерация начинается с {min(h)}, а не с {LOW}"
    # The range is not contiguous by design, but a hole is what this test is about: report it.
    holes = sorted(set(range(min(h), max(h) + 1)) - h)
    assert not holes, f"дыры в нумерации журнала: {holes}"


def test_matcher_finds_real_citations():
    """A matcher that matches nothing would pass forever. It must find known-good citations."""
    h = _headings()
    probe = ROOT / "verify" / "k132_sixthsub.py"
    assert probe.exists(), "контрольный файл исчез — проверьте имя, а не вывод"
    found = _cited(probe)
    assert found & h, f"шаблон не нашёл ни одной существующей ссылки в {probe.name}"


def test_every_cited_item_exists():
    h = _headings()
    missing = {}
    for p in _files():
        bad = sorted(c for c in _cited(p) if c not in h)
        if bad:
            missing[str(p.relative_to(ROOT))] = bad
    assert not missing, (
        "ссылки на пункты, которых нет в verify/README.md:\n"
        + "\n".join(f"  {f}: {b}" for f, b in sorted(missing.items()))
        + "\nЛибо пункт не записан, либо номер неверен. Второе чаще."
    )
