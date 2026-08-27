"""Repository paths. Scripts import this instead of hard-coding directories.

DATA     - challenge data, the data/ directory (not stored in git, see data/fetch.sh)
RESULTS  - run logs and saved predictions, the results/ directory
TUTORIAL - the organisers' repository holding the official metric implementation.
           Resolved in this order: $CYP_TUTORIAL, the pinned git submodule at the
           repository root, then a sibling clone next to this repository.

Importing this module creates the output directories it names. Scripts write straight
into them and none of them calls mkdir, so without this they would fail on savefig or
np.save in a fresh clone.
"""
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT / "data"
RESULTS = ROOT / "results"
PREDS = RESULTS / "preds"
LOGS = RESULTS / "logs"
FIG = RESULTS / "fig"
DOCS = ROOT / "docs"
TEXFIG = DOCS / "tex" / "fig"

for _d in (DATA, RESULTS, PREDS, LOGS, FIG, TEXFIG):
    _d.mkdir(parents=True, exist_ok=True)

# Most scripts were written against a string prefix, so keep that form too.
D = str(DATA) + os.sep
RES = str(RESULTS) + os.sep


def tutorial() -> Path:
    """Return the CYP-Challenge-Tutorial directory, having put it on sys.path.

    The submodule inside this repository is preferred over a sibling clone: it is
    pinned to a commit, so every machine sees the same metric implementation.
    """
    env = os.environ.get("CYP_TUTORIAL")
    candidates = [Path(env)] if env else []
    candidates += [ROOT / "CYP-Challenge-Tutorial", ROOT.parent / "CYP-Challenge-Tutorial"]
    for c in candidates:
        if (c / "evaluation" / "custom_scoring_functions.py").exists():
            if str(c) not in sys.path:
                sys.path.insert(0, str(c))
            return c
    raise FileNotFoundError(
        "CYP-Challenge-Tutorial not found. It ships as a submodule of this repository:\n"
        "  git submodule update --init\n"
        "or point CYP_TUTORIAL at your own clone."
    )
