"""Line citations in verify/README.md must still hold what their sentence names.

Item 308 found citations into src/submit.py that had stopped resolving without anyone
noticing. The journal cited line numbers, submit.py grew to three times its size, and every
diff that moved those lines was a diff to submit.py, so nobody re-read the journal. Its rule
is to cite a symbol, not a line. This guard makes the line citations that remain fail loudly.

Read as a line citation, in prose only:

    `src/x.py:123`, `src/x.py:123-130`       the colon form
    `src/x.py` line 123, lines 375 and 598   the file is the one after "of" or "in" if
                                             there is one, else the nearest .py named
                                             before it in the same item
    (`src/x.py:fn`, line 542)                the line must lie inside `fn`

A citation holds when every cited range contains at least one symbol its sentence names in
backticks, as text or by lying inside that definition. A sentence that names no symbol can
only be held to a file and a line that exist and are not blank: nobody means to cite a blank
line, and item 245 cited `src/submit.py:400` with no symbol until `SOLO` moved away and left
it blank. That is the weak guarantee, and the reason to name the symbol.

Deliberately not read:

- indented blocks and blockquotes, which are the journal's tables and quotations. Item
  308's audit table quotes the broken citations, and a scanner that read it would fail on
  its own case study;
- anything in double quotes: a quoted citation is a mention, not a use. Item 308 quotes two
  of item 307's retired self-citations that way, and a sentence that records a stale citation
  has to do the same, since a stale citation is precisely one whose line no longer holds. The
  quotes go outside the backticks, and they count only when they pair up across the
  paragraph: one stray quote character turns every quotation in that paragraph back into a
  citation, which fails loudly on the stale one it quoted rather than hiding a live one;
- a sentence that names a commit. Its line numbers belong to that commit's file, not to the
  working tree: in 61182bd the two calls item 202 cited as 593 and 598 stood at 600 and 605,
  and a record of where a line WAS is not a claim about where it is;
- a symbol right after "no", "not", "without" or "never", which is a claim of absence and
  not a needle. Item 257's pooled member fits with no `sample_weight`, and zero occurrences
  is the claim;
- backticked expressions, which are not symbols a line holds. Item 164 cited line 230 in the
  sentence that gives `(4 * четыре + ствол) / 5`; line 230 then held np.mean and none of the
  formula's names, so reading them as needles would have failed that citation on the day it
  was correct.

Not parsed: "at 553" without the word "line", which cannot be told apart from the numbers
this journal is full of.
"""
import ast
import bisect
import functools
import keyword
import os
import pathlib
import re
from dataclasses import dataclass

import pytest

ROOT = pathlib.Path(__file__).resolve().parents[1]
JOURNAL = ROOT / "verify" / "README.md"

ITEM = re.compile(r"\*\*(\d+)\.\s")          # whitespace after the dot: "**0.393" is a number
HEADING = re.compile(r"#+\s+(.*)")
LIST_ENTRY = re.compile(r"\s*(?:\d+\.|[-*])\s")
TICKS = re.compile(r"`[^`]*`")
QUOTES = re.compile(r'"[^"]*"|“[^”]*”|«[^»]*»')
PY = re.compile(r"(?<![\w/.-])((?:[\w.-]+/)*[\w-]+\.py)(?!\w)")
NUM = r"\d+(?!\.\d)(?:\s*[-–]\s*\d+(?!\.\d))?"
COLON = re.compile(rf":(?:({NUM})|([A-Za-z_]\w*))")
LINES = re.compile(rf"\b[Ll]ines?\s+({NUM}(?:\s*(?:,\s*and|,|and)\s+{NUM})*)")
OF_FILE = re.compile(r"\s+(?:of|in)\s+`?")
OWNER_GAP = re.compile(r"[`\s,(]*")          # between `x.py:fn` and the ", line N" it owns
SENTENCE_END = re.compile(r"[.!?][\"”»)*]*\s+(?=[A-ZА-ЯЁ*(\"“«\0])")
ABBREVIATION = re.compile(r"\b(?:e\.g|i\.e|cf|vs|etc)\.$")
COMMIT = re.compile(r"(?<![\w/.-])(?=[0-9a-f]*[0-9])(?=[0-9a-f]*[a-f])[0-9a-f]{7,40}(?![\w-])")
NAME = re.compile(r"[A-Za-z_]\w*(?:\.[A-Za-z_]\w*)*")
CALL_OR_ASSIGNMENT = re.compile(r"\s*(?:\(.*\))?\s*(?:=(?!=).*)?")
FILE_NAME = re.compile(r"\.(?:py|json|csv|md|npz|npy|txt|log|pdf|tex|toml|ya?ml|sh|ipynb)$")
ABSENCE = re.compile(r"\b(?:no|not|without|never)\s*$", re.IGNORECASE)


@dataclass(frozen=True)
class Citation:
    item: str                            # "item 202", or the section heading before items
    line: int                            # where in the journal, for the failure message
    path: str                            # as written: "src/submit.py", or a bare "submit.py"
    ranges: tuple[tuple[int, int], ...]
    symbol: str | None = None            # from the (`x.py:fn`, line N) form
    needles: tuple[str, ...] = ()


def scan(text):
    """Every line citation in the journal's prose, in order."""
    found, seen, current = [], [], None
    for item, para in _paragraphs(text):
        if item != current:
            seen, current = [], item      # .py files named so far in this item
        found += _citations(item, para, seen)
    return found


def problems(cite, locate):
    """Why a citation no longer holds; empty when it does. `locate` returns every file the
    cited path could mean, so that finding none and finding several are both reported."""
    hits = locate(cite.path)
    if len(hits) != 1:
        return [f"{len(hits)} files end with {cite.path}" if hits else
                "no such file (is the CYP-Challenge-Tutorial submodule checked out?)"]
    source = hits[0].read_text(encoding="utf-8")
    lines, spans = source.splitlines(), _definitions(source)

    def inside(name, a, b):
        return any(lo <= a and b <= hi for lo, hi in spans.get(name, ()))

    out = []
    for a, b in cite.ranges:
        where = f"line {a}" if a == b else f"lines {a}-{b}"
        if not 1 <= a <= b <= len(lines):
            out.append(f"{where} is past the end of a {len(lines)}-line file")
        elif cite.symbol:
            if not inside(cite.symbol, a, b):
                out.append(f"{where} is not inside `{cite.symbol}`")
        elif cite.needles:
            held = "\n".join(lines[a - 1:b])
            if not any(inside(s, a, b) or re.search(rf"(?<!\w){re.escape(s)}(?!\w)", held)
                       for s in cite.needles):
                out.append(f"{where} holds none of " + ", ".join(f"`{s}`" for s in cite.needles))
        elif not any(line.strip() for line in lines[a - 1:b]):
            out.append(f"{where} is blank, and the sentence names no symbol to check it by")
    return out


def _paragraphs(text):
    """(item, [(journal line, text), ...]) for each prose paragraph.

    An indented block after a blank line is a table or a log and a '>' line is a quotation;
    neither is prose, so neither is read."""
    item, para, code, after_blank = "preamble", [], False, True
    for n, raw in enumerate(text.splitlines(), 1):
        if not raw.strip():
            if para:
                yield item, para
                para = []
            after_blank = True
            continue
        if raw.startswith(("    ", "\t")) and (after_blank or code):
            code, after_blank = True, False
            continue
        code, after_blank = False, False
        heading, header = HEADING.match(raw), ITEM.match(raw)
        quotation = raw.lstrip().startswith(">")
        if heading or header or quotation or LIST_ENTRY.match(raw):
            if para:
                yield item, para
                para = []
            if heading:
                item = heading.group(1).strip()
            if header:
                item = f"item {header.group(1)}"
            if heading or quotation:
                continue
        para.append((n, raw.strip()))
    if para:
        yield item, para


def _citations(item, para, seen):
    starts, text = [], ""
    for _, line in para:
        starts.append(len(text))
        text += line + " "
    code = TICKS.sub(lambda m: "\0" * len(m.group()), text)    # quotes inside backticks are code
    balanced = (code.count('"') % 2 == 0 and code.count("“") == code.count("”")
                and code.count("«") == code.count("»"))
    quoted = [m.span() for m in QUOTES.finditer(code)] if balanced else []
    bounds = [0, *(m.end() for m in SENTENCE_END.finditer(code)
                   if not ABBREVIATION.search(code[:m.start() + 1])), len(text)]

    def mentioned(pos):
        return any(a <= pos < b for a, b in quoted)

    def cite(pos, path, spec, symbol=None):
        i = bisect.bisect_right(bounds, pos) - 1
        a, b = bounds[i], bounds[i + 1]
        if COMMIT.search(text, a, b):
            return                        # a record of that commit's file, not of this one
        ranges = tuple((int(x), int(y or x))
                       for x, y in re.findall(r"(\d+)(?:\s*[-–]\s*(\d+))?", spec))
        out.append((pos, Citation(item, para[bisect.bisect_right(starts, pos) - 1][0], path,
                                  ranges, symbol, _needles(text, a, b))))

    files = [m for m in PY.finditer(text) if not mentioned(m.start())]
    out, owners = [], []
    for f in files:
        c = COLON.match(text, f.end())
        if c and c.group(1):
            cite(f.start(), f.group(1), c.group(1))
        elif c:
            owners.append((c.end(), f.group(1), c.group(2)))
    for m in LINES.finditer(text):
        if mentioned(m.start()):
            continue
        of = OF_FILE.match(text, m.end())
        after = of and PY.match(text, of.end())
        owner = next((o for o in owners
                      if o[0] <= m.start() and OWNER_GAP.fullmatch(text, o[0], m.start())), None)
        before = [f.group(1) for f in files if f.end() <= m.start()] or seen
        if after:
            cite(m.start(), after.group(1), m.group(1))
        elif owner:
            cite(m.start(), owner[1], m.group(1), owner[2])
        elif before:
            cite(m.start(), before[-1], m.group(1))
    seen += [f.group(1) for f in files]
    return [c for _, c in sorted(out, key=lambda pc: pc[0])]


def _needles(text, a, b):
    """The symbols a sentence names in backticks, less the ones it says are absent."""
    names = []
    for m in TICKS.finditer(text, a, b):
        body = m.group()[1:-1].strip()
        name = NAME.match(body)
        if not name or not CALL_OR_ASSIGNMENT.fullmatch(body, name.end()):
            continue                      # an expression or a path, not a symbol
        sym = name.group()
        if (len(sym) < 3 or keyword.iskeyword(sym) or FILE_NAME.search(sym)
                or ABSENCE.search(text[a:m.start()])):
            continue
        names.append(sym)
    return tuple(dict.fromkeys(names))


def _definitions(source):
    """name -> [(first line, last line)] for every def, class and plain assignment."""
    spans = {}
    for node in ast.walk(ast.parse(source)):
        if isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names = [node.name]
        elif isinstance(node, ast.Assign):
            names = [t.id for t in node.targets if isinstance(t, ast.Name)]
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names = [node.target.id]
        else:
            continue
        for name in names:
            spans.setdefault(name, []).append((node.lineno, node.end_lineno))
    return spans


@functools.lru_cache(maxsize=None)
def _python_files():
    found = []
    for folder, subdirs, names in os.walk(ROOT):
        subdirs[:] = [d for d in subdirs if not d.startswith(".") and d != "__pycache__"]
        found += [pathlib.Path(folder, n) for n in names if n.endswith(".py")]
    return tuple(found)


def _locate_in_repo(cited):
    direct = ROOT / cited
    if direct.is_file():
        return [direct]
    return [p for p in _python_files() if p.as_posix().endswith("/" + cited)]


# A stand-in source file with hand-counted line numbers: 3-5 is `pooled`, 6 is blank,
# 7-8 is `combine`, and there is no line past 8.
SRC = """\
import numpy as np

def pooled(X, y):
    model = fit(X, y)
    return model

def combine(parts):
    return np.mean(parts, axis=0)
"""


@pytest.fixture
def check(tmp_path):
    src = tmp_path / "m.py"
    src.write_text(SRC, encoding="utf-8")

    def locate(cited):
        return [src] if cited == "m.py" else []

    def run(journal):
        cites = scan(journal)
        return cites, [(c.item, c.ranges) for c in cites if problems(c, locate)]
    return run


def test_a_moved_citation_is_reported_and_a_correct_one_is_not(check):
    cites, bad = check(
        "**22. Right.** `m.py` line 8 combines the members with `np.mean`.\n\n"
        "**23. Moved.** `m.py` line 5 combines the members with `np.mean`.\n")
    assert len(cites) == 2
    assert bad == [("item 23", ((5, 5),))]


def test_a_citation_wrapped_across_lines_is_still_read(check):
    """grep reads one physical line, and the journal wraps at about a hundred characters."""
    cites, bad = check(
        "**22. Wrapped.** The members are combined in `m.py` at\n"
        "line 8 by `np.mean`, and at\n"
        "lines 3 and 4 by `pooled`.\n")
    assert [(c.path, c.ranges) for c in cites] == [("m.py", ((8, 8),)),
                                                   ("m.py", ((3, 3), (4, 4)))]
    assert bad == []


def test_a_bare_line_number_binds_to_a_file_named_earlier_in_the_same_item_only(check):
    """Item 202's defect 2 names src/submit.py only in the item's heading, two paragraphs up."""
    cites, _ = check(
        "**22. Two defects in `m.py`.** Found by reading.\n\n"
        "**Defect 1.** The member pools (`pooled(X, y)`, line 3).\n\n"
        "**23. Next item.** No file is named here (`pooled(X, y)`, line 3).\n")
    assert [(c.item, c.path, c.ranges) for c in cites] == [("item 22", "m.py", ((3, 3),))]


def test_a_range_of_a_named_file_binds_to_that_file_not_to_the_nearest_one(check):
    """Item 308's form, "per lines 74-77 of" the organisers' validator, in a sentence that has
    already named a different file."""
    cites, bad = check(
        "**22. Of.** `other.py` is unrelated, and lines 7-8 of `m.py` hold `combine`.\n")
    assert [(c.path, c.ranges) for c in cites] == [("m.py", ((7, 8),))]
    assert bad == []


def test_a_symbol_qualified_citation_must_lie_inside_the_symbol(check):
    """Item 257's form: the line is cited as a line OF `_oof_one`, not as a line naming it."""
    cites, bad = check(
        "**22. Inside.** It fits (`m.py:pooled`, line 4).\n\n"
        "**23. Outside.** It fits (`m.py:pooled`, line 8).\n")
    assert [c.symbol for c in cites] == ["pooled", "pooled"]
    assert bad == [("item 23", ((8, 8),))]


def test_a_claim_of_absence_is_not_a_needle(check):
    """Trap (a): searching for the absent symbol reports a sound citation as rot, which is
    exactly what item 308's first audit pass did to item 257."""
    _, bad = check(
        "**22. Absent.** `m.py` line 4 fits with no `sample_weight`.\n\n"
        "**23. Present.** `m.py` line 4 fits with `sample_weight`.\n")
    assert bad == [("item 23", ((4, 4),))]


def test_quoted_and_tabulated_citations_are_mentions_not_uses(check):
    """Trap (b): item 308 documents broken citations by quoting them, in a table and in
    double quotes, and a scanner that reads those fails on its own case study."""
    cites, bad = check(
        "**22. Audit.** Item 21 wrote \"`m.py` line 99\", and the table quotes it again:\n\n"
        "    ссылка        символ    держится\n"
        "    `m.py:99`     pooled    НЕТ\n\n"
        "> `m.py` line 99 is quoted here as well.\n\n"
        "But this one is live: `pooled` is at `m.py:99`.\n")
    assert [(c.item, c.ranges) for c in cites] == [("item 22", ((99, 99),))]
    assert bad == [("item 22", ((99, 99),))]


def test_a_sentence_naming_a_commit_describes_that_commit_not_the_working_tree(check):
    """Trap (b) again, in the form item 308 uses to say where item 202's calls stood when
    item 202 was written."""
    cites, bad = check(
        "**22. Then.** In `abc1234`, `m.py` line 99 held `pooled`.\n\n"
        "**23. Now.** `m.py` line 99 holds `pooled`.\n")
    assert [(c.item, c.ranges) for c in cites] == [("item 23", ((99, 99),))]
    assert bad == [("item 23", ((99, 99),))]


def test_a_backticked_expression_is_not_a_needle(check):
    _, bad = check(
        "**22. Formula.** `m.py` (line 8) averages, so it is `(4 * members + trunk) / 5`.\n\n"
        "**23. Symbol.** `m.py` line 4 is where `combine` lives.\n")
    assert bad == [("item 23", ((4, 4),))]


def test_a_citation_naming_no_symbol_must_still_point_at_a_real_line(check):
    """Item 245 cited `src/submit.py:400` with no symbol, and line 400 became blank. Item 26
    is how a sentence records such a citation without being read as making it."""
    cites, bad = check(
        "**22. Blank.** The constant is in `m.py:6`.\n\n"
        "**23. Real.** The constant is in `m.py:4`.\n\n"
        "**24. Past the end.** The constant is in `m.py:40`.\n\n"
        "**25. Gone.** The constant is in `gone.py:4`.\n\n"
        "**26. Recorded.** The audit missed \"`m.py:6`\", which item 22 cites.\n")
    assert len(cites) == 4
    assert bad == [("item 22", ((6, 6),)), ("item 24", ((40, 40),)), ("item 25", ((4, 4),))]


def test_every_line_citation_in_the_journal_resolves():
    bad = [f"{c.item} (journal line {c.line}): `{c.path}` -- {why}"
           for c in scan(JOURNAL.read_text(encoding="utf-8"))
           for why in problems(c, _locate_in_repo)]
    assert not bad, ("line citations in verify/README.md that no longer hold what their "
                     "sentence names. Cite the symbol rather than renumbering (item 308); "
                     "if the sentence only quotes an old citation, put it in double quotes:\n  "
                     + "\n  ".join(bad))
