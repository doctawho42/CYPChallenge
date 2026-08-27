# Working on this together

Four people, four machines, one set of numbers that has to keep meaning the same thing.
This file is about the second part.

The repository was written by one person for one machine, and several of its habits only
worked under that assumption. They have been fixed, but the reasoning is worth knowing,
because most of the ways to break this repository are quiet ones — nothing errors, the
numbers just stop being comparable.

## Setup, once per machine

```bash
git clone git@github.com:doctawho42/CYPChallenge.git
cd CYPChallenge
git submodule update --init     # the organisers' repo, pinned to a commit
uv sync                         # exact versions from uv.lock
bash data/fetch.sh              # challenge data, ~7 MB
uv run python src/feats.py      # ~45 s, builds data/feats.npz and data/rows.csv
uv run pytest                   # should be 6 passed
```

Run everything through `uv run`, never a bare `python`. That is what guarantees all four
of us are on the same interpreter and the same library versions.

## The three things that break comparability

**1. The environment.** `uv.lock` is the contract. `scikit-learn` is capped below 1.9 on
purpose: 1.3.2 through 1.8.0 regenerate `results/preds/oof.json` bit for bit, and 1.9.0
does not — it shifts every prediction by about 0.135 pIC50 and doubles the apparent effect
of the mechanistic block. If you need a new package, `uv add` it, commit `uv.lock` in the
same commit, and tell everyone to `uv sync`. Never edit versions by hand.

**2. The split.** It lives in one place now, `cypsplit.py`. `tests/test_split.py` pins the
resulting fold vector to a golden digest and CI runs it on every push. If that test goes
red, the split moved and every table in the document is describing different folds than
the code now produces. Fix the code; do not edit the digest to make CI green. If you are
deliberately changing the split, update the digests in the same commit, say so in the
message, and re-run the pipeline.

Note that `mech2.py`, `mech3.py` and `mech4.py` cluster at 1024 bits rather than 2048, so
their folds have never matched the ablation table. That is preserved deliberately and
pinned by its own test. Those three scripts compare models against each other on their own
folds, which is sound; just do not put their numbers in the same table as `ablate.py`'s.

**3. The expensive artefacts.** `results/preds/oof.json` is 750 KB on a single line and
takes an hour to regenerate. Git cannot merge it and neither can you. It is marked binary
in `.gitattributes`, so a concurrent change is a clean "both modified" conflict rather than
a corrupted file.

The rule: **one person regenerates `oof.json` at a time, on `main`, and says so.** If you
hit a conflict on it, do not try to resolve it — take one side (`git checkout --theirs` or
`--ours`), then re-run `src/ablate.py` yourself if you need it to reflect your change.
`results/preds/oof.meta.json` records which split, which library versions and which commit
produced the current file, so you can always tell what you are holding.

## Branches

`main` is what the document quotes. Work on a branch, open a pull request, let CI run.

Keep a branch scoped to one question — "does X help?" — because that is the unit the
paired bootstrap answers, and it is the unit a reviewer can check.

## Reporting a result

The repository's own rule, from `README.md`, applies to conversation too: a difference
without an interval is not a result. When you report a number, say which split seed, which
feature set, and give the bootstrap interval. If it is a decision about what to submit, it
needs at least four seeds — a seed change alone moves macro ST-RAE by 0.016, which is more
than most effects we are measuring.

Negative results go in the repository, not just in chat. There are four already, and they
are some of the most useful content here — they are what stops the next person spending a
week on the same idea.

## Splitting the work

The pipeline forks cleanly after `feats.py`, so the natural division is by track rather
than by file:

- the regression track and decision layer (`ablate`, `score`, `decision*`);
- CYP2D6 and the range-restriction question (`range`, `rangectl`, `mech*`);
- the TDI classification track (`tdibase`, `tdiprob`, `f9`, `f10`);
- the document and the verification sweep (`docs/`, `verify/`).

Whoever touches `cypsplit.py`, `cyppaths.py` or `pyproject.toml` is changing something all
four tracks stand on — that is a pull request others read, not a quick push.

## Things that are still rough

- `docs/tex/figs.py` regenerates every figure in the document. It now reads and writes the
  right paths, but it needs `results/nn_seed0.npy` from `verify/f12_cvhard.py` (~10 min).
  `docs/build.sh` checks for it and tells you.
- The PDF is committed and is a build artefact, so it conflicts like one. Rebuild it from
  `docs/tex/` rather than merging it.
- Three numbers in the PDF are known to be wrong: the mechanistic block has 30 features
  and not 28, and there are fourteen verification scripts and not twenty-three. See the
  bottom of `verify/README.md`.
