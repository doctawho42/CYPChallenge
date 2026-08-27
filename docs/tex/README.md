# Document sources

The document itself is in Russian and stays that way: it is the team's working artefact,
and all four of us read Russian. Only the surrounding engineering documentation is in
English.

`main.tex` pulls in `preamble.tex` and sections `s01`…`s14` in order. Build with
`../build.sh`, or three runs of `xelatex main.tex` once the figures exist.

`figs.py` writes the figures into `fig/` as both pdf and png, at the page's natural size:
the 6.1-inch text block is filled by figures exactly 6.10 wide, so labels are not rescaled
on insertion and sit at the same point size as the body text. Colour is fixed per enzyme
across every figure: 1A2 blue, 2C9 brick, 2D6 green, 3A4 ochre. The palette was checked
for colour-blind discriminability and for contrast against the background.

`figs.py` needs, beyond the challenge data: `results/preds/oof.json` and
`results/preds/tdi_probs.json` (both committed), `results/seeds_all.csv` (committed), and
`results/nn_seed0.npy`, which is **not** committed — run `verify/f12_cvhard.py` (~10 min)
to produce it. `build.sh` checks for all of these before starting.

| Section | Topic |
|---|---|
| `s01` | the task and both metrics |
| `s015` | biochemistry: P450, IC50, TDI, how the four enzymes differ |
| `s02` | the data and the structure of the test set |
| `s03` | the baseline and three versions of the CYP2D6 failure |
| `s04` | generative model, likelihoods, instrument calibration |
| `s05` | identifiability, slope and depth of the curve |
| `s06` | missingness mechanism and factoring of the joint distribution |
| `s07` | architecture, coregionalisation, the mechanistic block |
| `s08` | loss function and training |
| `s09` | uncertainty |
| `s10` | decision layer for both metrics |
| `s11` | transductive layer |
| `s12` | evaluation protocol, robustness to the split |
| `s13` | what has been measured and does not work |
| `s14` | expected failures, what is untested, summary of numbers |
