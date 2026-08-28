# CYP Challenge

Our entry to the OpenADMET CYP Inhibition Blind Challenge: predicting inhibition of
CYP1A2, CYP2C9, CYP2D6 and CYP3A4 from molecular structure. Submission closes
3 November 2026.

The primary artefact is the document [`docs/CYP — модель и данные.pdf`](docs/) (in
Russian, and staying that way — it is the team's working document). Thirty-five pages: the
biochemistry of the task, what is in the data and what we measured from it, the model
with a justification for each part, the verification protocol, and a list of what has
already been tried and does not work. Everything else in this repository is the code
those numbers came from.

## Getting started

```bash
uv sync                             # environment, exact versions from uv.lock
git submodule update --init         # organisers' repo, holds the official metric
bash data/fetch.sh                  # challenge data (~7 MB, not stored in git)
uv run python src/feats.py          # features: fingerprint, descriptors, mechanistic block
uv run python src/score.py          # scores the saved predictions, macro ST-RAE 0.7729 / 0.7673
```

If you do not have [uv](https://docs.astral.sh/uv/) yet:
`curl -LsSf https://astral.sh/uv/install.sh | sh`. It fetches the pinned Python 3.12
itself, so nothing else needs installing. `make help` lists the pipeline targets.

Working on this with others? Read [CONTRIBUTING.md](CONTRIBUTING.md) first — it covers
who regenerates the expensive artefacts and how not to invalidate each other's numbers.

**Never used git or GitHub?** Start with [ONBOARDING.md](ONBOARDING.md) instead — it is
written in Russian, assumes nothing, and covers the traps specific to this repository.

## Layout

| Directory | Contents |
|---|---|
| `docs/` | the document and its sources (`tex/`, built with XeLaTeX), archive of earlier versions |
| `src/` | the pipeline: features, baselines, ablations, TDI, decision layer |
| `eda/` | data exploration and the figures for the document |
| `verify/` | twenty-two verification scripts, see `verify/README.md` |
| `results/` | run logs and saved out-of-fold predictions |
| `data/` | challenge data, not stored in git, see `data/README.md` |
| `tests/` | the golden-value guard on the cross-validation split |

Two modules at the root are shared by everything: `cyppaths.py` (paths, and locating the
organisers' repository) and `cypsplit.py` (the cross-validation split).

## Pipeline

Each step depends on the one before it:

1. `src/feats.py` — builds `data/feats.npz`: Morgan counts at r=2 (2048), 217 RDKit
   descriptors, 30 mechanistic features. The mechanistic block uses `src/pka.py`, a rule
   for estimating the pKa of the most basic centre.
2. `src/base.py` — baseline on fingerprint plus descriptors, five-fold cross-validation
   over Butina clusters at threshold 0.35.
3. `src/ablate.py` — six feature sets on the same folds, saves out-of-fold predictions to
   `results/preds/oof.json`.
4. `src/score.py` — metrics and a paired bootstrap over the saved predictions.
5. `src/range.py`, `src/rangectl.py` — testing the range-restriction explanation for CYP2D6.
6. `src/mech2.py`–`src/mech4.py` — three ways to bring in the primary screen, all three negative.
7. `src/tdibase.py`, `src/tdiprob.py` — baselines for the classification track.
8. `src/decision.py`, `src/decision2.py` — the decision layer for ST-RAE and its oracle ceiling.

Row 5 of the ablation grid — the joint likelihood — has its own chain, because its
predictions live in their own files and take about half an hour per arm to make:

9. `src/trunk.py` — shared trunk, two heads, `lambda_scr` as the switch on the screening
   loss term. `--mode twohead` gives the screen a free head; `--mode calibrated` routes it
   through the fixed instrument calibration instead, with no new parameters. One output
   file per mode, `results/preds/trunk_{mode}.json`; `--noise` adds an eta suffix. Both are
   committed, so everything below runs in minutes.
10. `src/trunkscore.py` — the lambda response with a paired bootstrap against lambda = 0.
11. `src/trunkdose.py` — the dose curve in lambda: whether the control is exact, what the
    damage is made of per enzyme, and how the arms compare with the boosting under a tilted
    label marginal.
12. `src/trunknoise.py` — the dose curve in *noise*. Degrades the screening channel in
    known steps so the four-enzyme ordering becomes four separate curves; needs the
    `--noise` runs of `src/trunk.py` first.

Two more that stand outside the grid:

13. `src/submit.py` — builds both submission files and refuses to write anything the
    organisers' validators reject.
14. `src/reweight.py` — scores under a test-like label marginal, since the test's own
    geometry cannot be reproduced by any re-split of the training data.

## Key numbers

Five-fold cross-validation over Butina clusters, threshold 0.35, seed 0.

| | 1A2 | 2C9 | 2D6 | 3A4 | Macro |
|---|---|---|---|---|---|
| ST-RAE, fingerprint and descriptors | 0.880 | 0.702 | 0.993 | 0.518 | 0.773 |
| ST-RAE, plus mechanistic block | 0.879 | 0.691 | 0.980 | 0.519 | 0.767 |
| Spearman, plus mechanistic block | 0.496 | 0.597 | 0.403 | 0.765 | 0.565 |
| sd of ST-RAE across four split seeds | 0.0109 | 0.0078 | 0.0052 | 0.0057 | 0.0019 |

One on the ST-RAE scale is what you get by predicting a single number for every compound.
CYP2D6 sits almost exactly there, and section three of the document is about why.

## How we work

All of this comes from measurement rather than from first principles; the reasoning is in
the document.

- Compare variants only on the same split and only with a paired bootstrap over
  compounds. A difference without an interval is not a result.
- One split is not enough: changing the seed moves macro ST-RAE by 0.016, which is more
  than the typical effect being measured. Any decision about what to submit is taken over
  at least four seeds.
- Our cross-validation is 0.153 harder than the real test set, by median similarity to the
  nearest training compound. Local numbers are a lower bound, but they also select for a
  different skill than the leaderboard pays for.
- Negative results are recorded and stay in the repository. There are four so far, and
  they are some of the most useful content here.
- The environment is pinned for a reason: scikit-learn 1.3.2 through 1.8.0 reproduce
  `results/preds/oof.json` bit for bit, 1.9.0 does not. See the comment in `pyproject.toml`.
