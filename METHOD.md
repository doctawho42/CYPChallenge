# Method

Our entry to the OpenADMET CYP Inhibition Blind Challenge: predicting inhibition of CYP1A2, CYP2C9,
CYP2D6 and CYP3A4 from SMILES, on two tracks — a regression on pIC50 for four isoforms, and a binary
classification of time-dependent inhibition for two.

The models are ordinary. Histogram gradient boosting, an exact Gaussian process, a ridge, and a
small two-headed network, averaged without weights. Nothing here is an architecture we would ask
anyone to adopt.

What is not ordinary is what surrounds them. This document is organised around four things that are:

**A loss derived from the competition's own metric.** ST-RAE charges zero for a prediction inside a
compound's published confidence band. That makes the point estimate the wrong training target, and
the right one — the prediction projected onto the band — turns out to be reachable by a stock
gradient booster under absolute error, fitting the metric's own gradient rather than a surrogate.
It is worth +0.0197 of rank against a floor of 0.0036, with the sign holding in all sixteen
per-enzyme cells.

**A label read correctly, late.** The classification target is not measured; it is computed from two
measured arms, and its rule folds into a conjunction. Folding it explains the whole difference
between the two scored endpoints with one number, and shows that our own classifier had been
solving a potency problem rather than a time-dependence one.

**An apparatus for telling an effect from a floor.** Noise floors measured per metric rather than
assumed; size-matched null blocks; deployment decisions pre-registered and committed before the
seeds finish; and an audit that found several of our own headline numbers were maxima and therefore
biased upward by construction. That apparatus closed the majority of our own ideas. The negative
results are the deliverable.

**A number written down before the reveal.** Three independent routes established that the test
regime cannot be checked from inside. A falsifiable prediction with a date is the only way we found
to turn that into a result — including when the first version of it had to be withdrawn.

Everything below carries the sample it was measured on, the number of split seeds, and the noise
floor it is judged against. Where a conclusion of ours was later overturned by our own measurement,
it is said so; that history is in the journal deliberately and is quoted here for the same reason.

The code is at the repository root; `verify/README.md` is the journal, numbered to item 303, and is
the primary record. This document is a reading of it, not a substitute.

---

## The loss the metric implies

ST-RAE scores a prediction inside the compound's published confidence band `[lo, hi]` at exactly zero and one outside it at the distance to the nearest edge; the denominator puts a constant predictor at `mean(y)` through the same rule (`evaluation/custom_scoring_functions.py`). The point estimate is the usual training target, under squared error, with the band appearing only at scoring time. We trained against the band.

**Shape.** The loss the metric charges per compound is

L(p) = max(0, lo − p, p − hi) = |p − clip(p, lo, hi)|

— piecewise linear, slope −1 below `lo`, zero on `[lo, hi]`, +1 above `hi`: a flat bottom of width `hi − lo`, no curvature where it is differentiable.

**Gradient identity.** dL/dp = sign(p − clip(p, lo, hi)), exactly the gradient of |p − t| at t = clip(p, lo, hi). At the kinks it is containment rather than equality: the dead zone's subdifferential is [−1, 0] at p = lo and [0, 1] at p = hi, both inside the surrogate's [−1, 1], and all three contain zero, so the stationary sets agree.

**Majorisation.** For any t in [lo, hi], |p − t| ≥ L(p), with equality at t = clip(p, lo, hi) — three cases: above `hi` both sides are p − hi at t = hi and the left grows as t falls; below `lo` symmetric; inside the band the right side is zero. Absolute error against the clipped target therefore majorises the metric's loss and is tight at the current prediction, so "reproject, refit" is majorise–minimise:

Σ L(p⁽ᵏ⁺¹⁾) ≤ Σ |p⁽ᵏ⁺¹⁾ − t⁽ᵏ⁾| ≤ Σ |p⁽ᵏ⁾ − t⁽ᵏ⁾| = Σ L(p⁽ᵏ⁾)

— majorisation, a refit no worse than the incumbent on its own problem, tightness. The loss cannot increase, subject to the middle step assuming the fit attains its minimum, which boosting does only approximately.

Nothing new has to be built: a stock booster under `loss="absolute_error"` fitted to `clip(p, lo, hi)` descends the metric's own gradient rather than a surrogate for it, and `DZ_KW` in `src/submit.py` is the baseline's pins plus `absolute_error`. The reduction is exact only for the absolute-loss members; the Gaussian process and the ridge are refitted on the same target with a square outside the band, an approximation measured separately (item 149).

**The trap.** The target cannot come from a model's predictions on its own training rows: an overfitted model puts nearly every training row inside its band, the target equals the prediction, the gradient is zero, and the pass becomes the identity — while the run finishes cleanly and prints a plausible table. Projection therefore uses out-of-fold predictions (`_dz_oof`), which doubles the cost of building a submission. The remaining leak, that a training row's target saw models trained on other folds, can only make the target easier; rank against the truth is what is scored.

**Loss, or target?** The pass changes both, and `src/abldead.py` separates them — per-enzyme booster, four seeds, macro after the affine pair:

| arm | ST-RAE | rank | rank vs squared |
|---|---|---|---|
| squared on the point label | 0.7156 | 0.5630 | — |
| absolute on the point label | 0.7136 | 0.5646 | +0.0016, sign 3/4 seeds |
| absolute on the clipped target | 0.6887 | 0.5969 | +0.0339, sign 4/4 seeds |

Against a fixed-seed floor of 0.007 the loss change alone is null; the reprojection is the whole effect. Contrast item 77: L1 on the label, once the largest raw gain here at −0.0455 of ST-RAE, left 0.0019 after the affine pair. The dead zone survives that pair because its target depends on per-compound `lo` and `hi`, so it reorders compounds, and one increasing function of p cannot manufacture a reordering. Only the first pass pays: at seed 0, passes two and three give 0.5904 and 0.5885 against 0.5923.

**Measurement.** `verify/k58_dzsubmit.py`, seeds 0–3, the submission's own code, five members (item 213):

| arm | ST-RAE | rank | 1A2 | 2C9 | 2D6 | 3A4 |
|---|---|---|---|---|---|---|
| no pass | 0.6658 | 0.6145 | 0.5433 | 0.6586 | 0.4648 | 0.7911 |
| pass in four of five | 0.6483 | 0.6297 | 0.5568 | 0.6820 | 0.4744 | 0.8056 |
| pass in all five | 0.6459 | 0.6342 | 0.5615 | 0.6855 | 0.4820 | 0.8078 |

+0.0197 of rank and −0.0199 of ST-RAE against a four-seed macro floor of 0.0036, sign holding four seeds of four, all sixteen per-enzyme cells positive (+0.0181, +0.0270, +0.0172, +0.0167). Rank and metric improve together, which is rare here.

**Every member, or it does not transfer.** Four members of five gives +0.0152, all five +0.0197; the fifth member's +0.0045 holds its sign on four seeds of four but sits at the 0.0036 floor, so the sign carries it, not the size. Applied to the boosters alone in the four-member ensemble the pass scored 0.6001 against 0.6029 without it — inside the floor and on the wrong side (item 149). Two earlier attempts on a single member arrived at nothing (standalone +0.031 and +0.029 reaching the ensemble as −0.0008 and +0.0014, items 176 and 182); here +0.0320 standalone delivers +0.0045, because the objective changed under every member at once rather than one member being replaced by a better one.

**What we got wrong.** `abldead.py` pre-registered that the gain must be largest where bands are widest — 29.4 % of CYP3A4 rows exceed 1.0 pIC50 of band width against 9.9 % on CYP2D6. The order inverted, and on four seeds the last cell inverts further: CYP2D6 +0.0523, CYP2C9 +0.0407, CYP1A2 +0.0402, CYP3A4 **+0.0025** — inside CYP3A4's own floor of 0.0033 and holding its sign on only three seeds of four. (Item 148 prints +0.0552 / +0.0431 / +0.0407 / +0.0056, which are seeds 1–3; the four-seed figures are recomputed from the saved predictions and strengthen the point rather than weakening it — on the enzyme with the widest bands the dead zone buys nothing measurable at all.) On CYP3A4, 73.9 % of predictions on wide-banded rows already lie inside the band, where clip(p) = p and there is nothing to win. The dead zone does not stop effort being wasted on wide bands; it stops us paying for the last tenth of a pIC50 on narrow ones.

## The classification label is a conjunction, and that determines what is achievable

The TDI endpoint is not measured. `is_TDI` is computed from two measured arms — the direct-inhibition pIC50 and the pre-incubation arm — by a piecewise rule: above pIC50 4 in the direct arm, positive means the shift exceeds two-fold; below it, positive means the pre-incubation arm exceeds 4.301. The second branch yields what the organisers call "inferred positives", because a true direct value of at most 4 makes a pre-incubation arm above 4 + log10 2 sufficient to guarantee a greater-than-two-fold shift.

The two branches fold into one conjunction:

    is_TDI  <=>  (Delta > log10 2)  AND  (pi_TDI > 4 + log10 2),    Delta = pi_TDI - pi_dir

The proof is two lines. If `pi_dir > 4`, then `Delta > log10 2` already implies `pi_TDI > 4.301`, so the second conjunct is slack. If `pi_dir <= 4`, then `pi_TDI > 4.301` already implies `Delta > 0.301`, so the first is slack. There are no branches — one conjunction seen from two sides. Both forms reproduce the published label exactly, 2334/2334 rows on CYP3A4 and 1493/1493 on CYP2D6, and `src/tdi.py` re-checks the equivalence elementwise on every run (item 234).

The folded form makes visible what the piecewise form hides: the shift is necessary for every positive, and the second conjunct is a gate that can only strike rows out. Measured over the same 3827 rows, positives with no shift: 0. Positives with the gate closed: 0.

| | gate closed | positive rate, gate open | oracle MCC, gate alone | oracle MCC, shift alone | gate recovered |
|---|---|---|---|---|---|
| CYP3A4 | 39.8 % | 0.543 | 0.5667 | 0.6875 | 53.7 % |
| CYP2D6 | 5.8 % | 0.230 | 0.1302 | 0.9010 | 1.2 % |

The first four columns are properties of the label, computed on the true arms, with range 0.0000 over four split seeds; the recovery column is four seeds (item 238).

That table is the whole difference between the two endpoints. On CYP3A4 the gate strikes out 39.8 % of rows and recognising them is a pure potency question. On CYP2D6 it strikes out 5.8 %, so there is almost no potency sub-problem to win and the label is very nearly the shift alone, oracle MCC 0.9010. AUC 0.745 against 0.588 follows from that and needs no other explanation.

The consequence for our own model is uncomfortable. On CYP3A4 the predicted gate alone — a bare HistGB regression on the pre-incubation arm, cut at 4.301 — scores MCC 0.3043 (four seeds, range 0.0203) against 0.315 for the whole structural classifier with its plug-in threshold. Adding the predicted shift brings it to 0.3317, so the shift contributes +0.027. Both gaps are inside the CYP3A4 MCC floor of 0.0281 (item 235, four seeds). The comparison is exact at a fixed threshold rule: 0.315 is the structural classifier under the same plug-in cut. To within that floor, that classifier is a potency threshold: it finds compounds too weak to clear 4.301 and calls them negative, which is automatically correct on 39.8 % of the set. It is not modelling time-dependence. CYP2D6 is the mirror image — gate 0.0015, shift 0.0910, both together 0.0928, against the classifier's 0.118.

In conjunction coordinates the two error sources price separately. A predicted gate with the true shift scores 0.7563 on CYP3A4 and 0.7722 on CYP2D6; a true gate with the predicted shift scores 0.5285 and 0.1334. Error in the gate costs 0.244 and 0.228 of MCC from a perfect 1.0000; error in the shift costs 0.471 and 0.867. The shift is the larger loss on both, but on CYP3A4 a quarter of the loss sits in the tractable half.

The view also predicts a result we would otherwise mis-read. Fitting the gate's cut out of fold on CYP2D6 (optimum 4.69 to 4.98, not 4.301) improves the gate, MCC against the true gate rising from 0.078 to 0.182, and makes the label worse, 0.0015 to −0.0330. The gate's only job is to subtract, so a better-centred but still noisy gate strikes out more true positives than true negatives.

The history belongs with the result. Item 229 read the mechanism off the piecewise form and attributed CYP3A4's score to "the rule's second branch". Item 232 corrected the arithmetic — per-branch AUC 0.664 and 0.773, branch one supplying 680 of CYP3A4's 764 positives — but kept the frame, which was the error: conditioning on `pi_dir` cuts the population in a way the rule has no counterpart for. Item 234 folded the rule and superseded both. Worse, item 23 had stated the conjunction in prose 211 entries earlier — "the label is 'potent AND shifted'" — with a knife-edge argument neither successor recovered: the median Delta over CYP3A4 actives is +0.294 against a cutoff of 0.301, so half the actives sit within a hundredth of the line, and a predictor of Delta can be good while a predictor of the label looks worthless. Two wrong accounts, settled by one boolean comparison.

This diagnosis is what the submitted model was changed in response to. The classifier that now ships does not learn the label directly at all: it builds both conjuncts — a gate probability from the pre-incubation arm predicted by a four-member ensemble with the dead-zone pass, times a shift probability from a classifier on Δ > log₁₀2 — and thresholds the calibrated product. Acceptance was fixed in writing before the seeds finished and the conditions were met: macro MCC +0.0127 against a floor of 0.0076, sign holding in six cells of eight, with the shuffled-shift control collapsing to −0.1453 at zero of eight (items 244, 250). Per enzyme it clears nothing — +0.0083 on CYP3A4 against 0.0281 — so the claim is macro and is written that way. It moves 17 to 20 per cent of the submitted labels for a gain of 1.7 floors, which is recorded beside it.

## What we made of the released data

The release gives us dose-response curves with confidence bands, a single-concentration screen, and a time-dependent-inhibition arm. Only the first is the scored target; most of what follows came out of the other two. Every external source we tried is closed.

**The screen, joined through a shared latent curve.** The single-concentration file is 17504 rows: 4376 molecules on all four enzymes at 4.95·10⁻⁵ M. 11509 of those readings have no dose-response curve, 1.8 times the 6525 labelled pairs (item 83). They cannot serve as labels, because curves were run on what the screen flagged: on CYP2D6 every curve compound is screen-active against 5.7 per cent of the rest.

That selection is what makes a joint likelihood legitimate rather than a hazard. The missingness is not at random, but the selection variable is observed: written through a latent activity, the factor p(R | y_screen) carries no parameters and leaves the gradient, provided the screen enters as an observation and not as an input. It is a Heckman correction whose first stage is measured rather than modelled. Simulation with known truth (40 repeats, 6000 compounds, curves on the top 28 per cent by screen) sizes it: fitting the selected labels alone recovers the coefficients at a calibration slope of 0.731, the joint fit at 1.012. Select on the latent value instead, hiding the screen for the unselected, and the joint scheme fails as badly as the naive one, 0.681.

In the model this is a shared trunk with two heads — pIC50 from the curves, log2fc from the screen — and no screening column at inference, since the 750 test molecules have none. It multiplies the trunk's supervision by 2.76. Against its own control, the identical network with the screening head receiving no gradient, the channel is worth −0.0264 of ST-RAE after the affine pair on four seeds (p = 0.002) and +0.0350 of rank; the architecture alone gives −0.0008 (items 79, 80, 121). It ships as a fifth ensemble member, +0.0054 of rank on four seeds against a macro floor of 0.0036.

**The confidence bounds as a training target**, not only as a scoring device — the dead-zone projection described above. Band width is an isotonic function of the label to within 3 to 7 per cent of its variance (item 114), so the bounds carry potency information the point estimate does not present directly. In all five members the projection is worth +0.0197 of rank on four seeds, sign holding in all sixteen per-enzyme cells (item 213).

**The pre-incubation arm as a second regression target.** It has its own confidence bounds, so the dead zone applies. Trained that way, the four-member ensemble plus the dead zone lifts CYP3A4's AUC against the true TDI gate from 0.8714 to 0.8928 (sign 4/4, four seeds) but moves MCC by only +0.0161 against a floor of 0.0281: real in the ordering, not demonstrable in the scored metric, so nothing is deployed (item 243).

**The 1238 rows outside the model.** The TDI table holds 1238 molecules absent from the modelled row set. Their `is_TDI` is a placeholder — all negative, while 1048 of them (84.6 per cent) clear the rule's second branch, pi_TDI > 4.301 (item 231). Their pi_TDI is a real measurement and 53 per cent more rows for CYP3A4's arm regressor: +0.0051 on the regression route, −0.0055 on the classifier route, four seeds, both inside the 0.0281 floor (item 240).

**The mechanistic block, and why it works.** Thirty derived columns — nitrogen and acid counts, protonation state under a pKa rule, distances from a basic nitrogen to an aromatic ring — worth +0.0163 of rank on four seeds, +0.0313 with rows reweighted toward the test set's regime (item 119). Two measurements from opposite directions say what carries it. A full unit of Gaussian noise in the pKa estimate, flipping the basicity indicator on 4.2 to 4.9 per cent of the set, moves macro rank −0.0020 on four seeds, below the 0.0036 floor — so a trained pKa predictor is bounded at nothing (item 230; on seed 0 alone that run read as a gain and was withdrawn). And permuting the ten salt-bridge columns costs +0.1545 of rank on CYP2D6 above a size-matched random ten-column block, while the 20 charge-derived descriptor columns fall below their own size-matched null on all four enzymes (item 239). The block delivers geometry of the basic centre, not the protonation call and not the charge.

**Tried from outside, and closed.**

| Source | Measurement | Item |
|---|---|---|
| Pretrained representations, three checkpoints | chemprop `rdkit2d` concatenated: −0.0075 rank; MoLFormer-XL via in-fold PCA-64: 0.5637 against 0.5651 | 61, 117, 154 |
| ChEMBL, 8068 compounds | a raw gain of 0.0128 becomes +0.0068, harmful, after the affine pair; seed 0 | 61, 77 |
| NCGC panel, 54177 rows | +0.0048 macro rank at the best of three weights tried, below the 0.007 floor | 144, 157 |
| Quantum descriptors, 10 columns | null against potency and against the TDI shift, where the shuffled block (+0.0041) beats the real one (−0.0018) | 186, 242 |

The binding constraint is the 1285 to 2335 labelled rows per enzyme, and the readings that relieve it were already in the release.

## Uncertainty, in the units the metric pays in

Conventional uncertainty quantification puts an interval around the prediction. Under ST-RAE that
is the wrong object. A prediction anywhere inside the compound's published confidence band scores
exactly zero, so distance within the band is not paid for, and an interval that straddles the band
edge says nothing about what it will cost. The quantity the metric actually cares about is the
per-compound probability of not paying at all:

```
P(hit_i) = P( lo_i <= ŷ_i <= hi_i )
```

This is the natural partner to the dead zone. One trains predictions *into* the band; the other
says whether they landed. Both are derived from the definition of the metric rather than imported
from elsewhere.

The problem is well posed because the organisers' own bands make it so. Out of fold on the submitted
model, **31.3 % of predictions already score exactly zero** (per enzyme: 0.247, 0.439, 0.198, 0.369).
And the rate is strongly structured by a quantity the model itself produces. Band width is a
function of the label to within three per cent of its variance, and the relation is *decreasing* —
weakly active compounds have wide bands. On CYP3A4, by quintile of predicted potency:

| quintile | mean prediction | mean band width | hit rate | mean loss |
|---|---|---|---|---|
| 1 | 3.08 | 2.000 | 0.788 | 0.101 |
| 5 | 5.14 | 0.209 | 0.161 | 0.309 |

A five-fold spread in the probability of scoring zero, visible entirely from the model's own output.

### What it delivers, and on how much of the problem

| enzyme | base rate | AUC | Brier | Brier of a constant at the base rate |
|---|---|---|---|---|
| CYP1A2 | 0.247 | 0.534 | 0.1928 | **0.1861** |
| CYP2C9 | 0.439 | 0.597 | 0.2490 | **0.2463** |
| CYP2D6 | 0.198 | 0.486 | 0.1635 | **0.1585** |
| CYP3A4 | 0.369 | **0.750** | **0.1814** | 0.2328 |

**On three enzymes of four the Brier score is worse than quoting the base rate.** The estimator adds
noise there and we report it as such. On CYP3A4 it is real: AUC 0.750 and a Brier score that beats
the constant by 0.051.

The mechanism is a property of the assay rather than of the model, and it predicts which enzyme
works before the model is fitted. CYP3A4's bands span a nineteen-fold range between the tenth and
ninetieth percentiles; the other three span five to eight. Its predictions also span twice the
spread of CYP2D6's. Where the band barely varies there is nothing for a hit probability to
discriminate, and no model will change that.

The better form of the estimator is a single constructed feature rather than a learned one — the
signal-to-noise ratio `z = ŵ(ŷ) / (2·σ̂)`, calibrated to a hit probability by isotonic regression.
On CYP3A4 it beats the learned classifier on both metrics and reduces the whole object to one
interpretable number: the correlation of `z` with hitting is **+0.435**.

### Two failures, reported because the apparatus is the point

The signal-to-noise feature first scored AUC 0.485 on CYP3A4 — worse than the classifier it was
meant to improve. The cause was that `IsotonicRegression` defaults to requiring a non-decreasing
fit, while the width-potency relation is decreasing; the fit collapsed to a constant of range
exactly 0.000, and `z` degenerated into an inverse transform of the residual scale alone. It was
caught by a control rather than by reading the code: the correlation between predicted band width
and hitting came out at ±0.01 on all four enzymes, which cannot happen for a monotone function of an
informative feature.

A second head predicts the score itself, and reports macro ST-RAE of 0.6594 against a true 0.6600.
That reads as a strong result and is not one. Predicting the **mean loss** out of fold does equally
well — 0.7657 / 0.5790 / 0.8523 / 0.4426 against true 0.7656 / 0.5793 / 0.8520 / 0.4430 — because
the folds are exchangeable, not because the model knows anything about individual compounds. The
control belonged in the first version and was absent; it is the same defect this project had
identified in a different measurement two days earlier.

What the section claims is therefore narrow: on one of four enzymes, a calibrated per-compound
probability of scoring zero under the competition's own metric, at AUC 0.750. It is not a macro
result, it is not deployed, and it is not a prediction of the leaderboard score.

## The apparatus

What this project offers is methodological rather than architectural: an instrument for telling an effect from a floor, pointed at its own ideas, which then closed most of them. The journal, `verify/README.md`, runs to item 303; negative results are the majority.

**Noise floors are measured per metric, not assumed.** Item 70 measured the pipeline's chaotic sensitivity at a fixed seed and fixed folds: two arms whose tilt targets coincide on CYP2D6 to within 0.01 still differ there by 0.0071 in ST-RAE. That 0.007 is the fixed-seed macro floor. Changing only the split seed moves macro ST-RAE by 0.016 (`verify/f3_seeds.py`), more than most effects measured here.

Item 165 found the macro floor had been used throughout for per-enzyme claims and cannot be: a macro average over four enzymes is quieter than its parts. Spread across four seeds at a fixed arm:

| | rank, sd across 4 seeds |
|---|---|
| CYP1A2 | 0.0061 |
| CYP2C9 | 0.0071 |
| CYP2D6 | 0.0049 |
| CYP3A4 | 0.0033 |
| macro | 0.0036 |

The classification track had no floor until item 235, where a four-seed sweep of the deployed arm gave one as the range across seeds: 0.0281 on CYP3A4, 0.0419 on CYP2D6, 0.0076 macro. The macro MCC floor landed within a thousandth of the macro rank floor of 0.007 — a different metric on a different track. That is a coincidence worth noticing and not worth theorising about.

**Only rank survives.** The submission fits an affine pair — a shrink and a shift — to the metric per fold before writing predictions out; being monotone, it undoes any intervention that changed only scale or location. Item 77 refitted that pair out of fold on five arms measured in one day: raw gains of 0.0084 to 0.0455 of ST-RAE became at most 0.0033, all under the 0.007 floor, and three of the five changed sign to harmful. On four seeds the largest holds its raw sign 4/4; after the pair it is 0.0019 and the sign breaks on seed 3. Item 117 re-read the sharpest of those five by rank: a concatenated chemprop embedding read as a raw gain of 0.0084 was a loss of 0.0113 after the pair, −0.0075 by rank, and −0.0182 under test-regime weights.

**Named feature blocks are scored against a size-matched random block** — same width, same source, ten draws — because bigger blocks cost more for being bigger. Item 239 answered an outside proposal for DFT-quality ESP charges. Permuting the 20 charge-derived DESC columns costs 0.0228 of rank on CYP3A4 — a channel, until twenty random DESC columns cost 0.0391. Against its own null it is negative on all four enzymes: −0.0058, −0.0076, −0.0013, −0.0163. The control is not merely a killer — the 10 salt-bridge features sit +0.1545 above their null on CYP2D6.

**Deployment decisions are pre-registered before the seeds finish.** Item 233 fixed three conditions for putting Platt calibration into the TDI path, with seed 0 in hand and seeds 1 to 3 still running: a positive mean macro-MCC gain, the sign holding in at least 6 of 8 cells, and a positive Platt slope in all 40 folds. Item 235 checked them in order — +0.0133, 7 of 8, minimum slope +0.0836 — and it shipped. It also bounded the claim: +0.0235 on CYP3A4 against a floor of 0.0281 and +0.0031 on CYP2D6 against 0.0419, so there is a macro result and no per-enzyme one. Item 244 pre-registered a bundle of sub-floor gains the same way, requiring it to clear the 0.0076 macro floor none of its parts cleared. It failed its own audit before the run reported.

**That audit is item 245.** A sign count is evidence only for an arm whose null expectation is zero; an arm that is itself a maximum has a positive null expectation by construction. Three of the project's own numbers were maxima:

- item 235's threshold oracle — best of 91 thresholds;
- item 241's "best MCC" column — best of 3 threshold rules;
- item 218's solo-GP constant — best of 4 enzymes, then best of 31 subsets, and the only one that ships.

Two of item 244's three bundled arms were that argmax column: +0.0187 (6/8) and +0.0162 (7/8) at the maximum, but −0.0031 (4/8) and −0.0123 (3/8) under the plug-in rule the submission uses. The inflation from best-of-three, +0.0121 and +0.0152, is larger than either apparent gain. Item 218's +0.0098 of rank on CYP3A4 becomes +0.006 to +0.008 once corrected, about +0.002 in macro at quarter weight; it still ships, but no longer as a leaderboard gain.

**Permutation and shuffle controls are standard, and they fire.** Item 242 measured a ten-column quantum-chemistry block against the TDI shift, four seeds, eight cells: the block raw scores −0.0018 at 3/8, its out-of-fold residual +0.0035 at 5/8, ten Gaussian-noise columns −0.0033 at 3/8 — and the shuffled block +0.0041 at 4/8, higher than either real arm. The block closed on its own control.

**The log overturns itself on the record.** Item 61 said the representation was never the bottleneck, item 68 said 61's headline was wrong, item 117 reinstated 61. A CYP2D6 gain of +0.0174 at seed 0, three and a half times that enzyme's floor, became +0.0006 at 1/4 over four seeds (items 227, 230). Items 229 and 232 both read the TDI label's mechanism wrongly, and item 234 killed both with one boolean comparison over 3827 rows. Items 68 and 118 stand marked withdrawn and retracted; items 135 and 161 retract earlier work (item 134's causal half, and a pairwise dead-zone run). The negative results are the deliverable; the apparatus is what makes them worth reading.

## What we tried and closed

The project's journal (`verify/README.md`) runs to item 303, and most of them record something that did not work. Two conventions govern the catalogue. A gain counts only after the affine pair — a per-fold shrink-and-shift fitted to the metric — has run: item 77 found five separately measured raw improvements collapsing into a band of 0.013 once it did, three of them reversing sign. And every claim is quoted against a noise floor: 0.007 macro ST-RAE at fixed seed (item 70), 0.0036 macro rank and 0.0033 to 0.0071 per enzyme over four seeds (item 165), 0.0033 to 0.0071 per enzyme (item 165), 0.0281 on CYP3A4 and 0.0419 on CYP2D6 for MCC (item 235), normally over four split seeds.

| Idea | What was measured | Why it failed | Item |
|---|---|---|---|
| **External data sources** | | | |
| ChEMBL set, 8068 compounds | −0.013 macro raw, sign held on two seeds | +0.0068 after the affine pair: the pair already removes the offset and variance the extra rows corrected | 61, 77 |
| NCGC panel, 18711 fitted-AC50 rows | −0.0176 macro rank; +0.0048 at weight 0.1, floor 0.007 | at full weight the panel is 74 % of the training table, so a foreign protocol offset by +0.44 to +0.87 dominates the objective | 153, 157 |
| The organisers' own second CYP release | 14 of 1340 molecules overlap our training set, 4 carry a CYP3A4 label | a protocol offset is estimable only on shared molecules; overlap is a property of the library, not of its size | 201 |
| 1238 unread CYP3A4 measurements | manufactured label keeps 98.6 % of rank; −0.0108 against a floor of 0.0033, 4/4 seeds | good but not good enough; and 1048 of them are flagged `False` while their own pre-incubation arm clears the rule's cutoff — placeholders, not measurements | 216, 231 |
| **Representation learning** | | | |
| Pretrained encoders: chemprop, ChemBERTa, MoLFormer-XL | concatenated, raw −0.0084 but +0.0113 after the pair, rank −0.0075; MoLFormer alone 0.4370 against 0.5651 | three checkpoints add nothing to the 2295-column matrix; the reported win was raw ST-RAE and stood against the document for a month before item 117 withdrew it | 61, 68, 117, 154 |
| Free-Wilson on the fifty informative bits | −0.15 to −0.13 rank in every similarity stratum, four seeds | those bits are the shared coordinating pharmacophore, what analogues have in common; additivity needs what distinguishes them | 191 |
| Two-view low-rank completion, decoder from the screen | −0.016 macro at rank 3, −0.064 at rank 2 | costs everywhere except CYP2C9 (+0.0227, floor 0.0071), which has the fewest labels (1285) — the one place pre-registered to gain | 193 |
| **Ensemble composition** | | | |
| Neural trunk as a fifth member | +0.031 rank standalone, −0.0008 in the ensemble | it improved by fitting the target the other four already fit; a mean pays for disagreement | 174, 176 |
| Screening-table arm as a fifth member | +0.029 standalone, +0.0058 plain and +0.0014 in the best composition | its errors correlate with the ensemble's at 0.935 to 0.969 — wrong on the same compounds | 177, 182 |
| Negative correlation learning | member error correlation 0.939 to −0.302; ensemble +0.0002, −0.0067, −0.1813 | diversity is manufacturable to any degree and does not pay | 194 |
| **Post-processing and thresholds** | | | |
| Any per-compound correction on the predictions | monotone budget above the affine pair: 0.0076 macro, 0.0002 on CYP3A4 | the oracle correction is the model's own residual at r = 0.87–0.95, so it is information that belonged in the model | 128 |
| Best-of-N threshold rules | argmax over three rules inflates by +0.0121 and +0.0152 | a maximum has a positive null expectation, so a consistent sign is a tautology; under the deployed rule the arms were 4/8 and 3/8 | 244, 245 |
| **Feature blocks** | | | |
| Quantum block (HOMO, LUMO, gap, Fukui f−) | null against potency and against the TDI shift, four seeds each; the shuffled block scores higher | passes the derivability gate on paper and fails in measurement — the gate is necessary, not sufficient | 186, 242 |
| DFT-quality ESP charges | the 20 charge columns already present fall below a size-matched random block on all four enzymes, −0.0058 to −0.0163 | wider blocks cost more for being wider, so the control must be width-matched rather than permuted; this is refinement, not addition | 239 |
| Analogue-series shrinkage layer | within-series spread is 0.29 to 0.58 of what chemistry allows | already over-smoothed inside series, so shrinking compounds the error; settled from existing predictions in half an hour, with no labels | 133 |
| **TDI classification** | | | |
| Mechanism-based structural alerts | the strongest, a tertiary-amine indicator, reconstructs from DESC+MECH at R² 0.998 out of fold | already in the matrix; and the open-gate subsets that flatter alerts condition on a quantity that does not exist at prediction time | 23, 236 |

Three of the ensemble rows are one finding. Items 176, 182 and 191 are consecutive standalone gains — +0.031, +0.029, and a stratum profile deliberately opposite to the Gaussian process's — arriving at the ensemble worth −0.0008, +0.0014 and −0.0025, with member-to-ensemble error correlations of 0.90 to 0.97. Item 194 then drove the members to anti-correlation and the ensemble did not move. Four independent measurements say the same thing: the members are limited by the 1285 to 2335 labelled rows they all share, not by a common inductive bias. Model diversity does not cure a shortage of data.

The classification track has the mirror-image constraint. Twice we found a real improvement to the ordering: training on `Delta > log10 2` and ranking `is_TDI` with it, +0.0196 AUC on CYP2D6 at sign 4/4 (item 241); and the regression ensemble plus dead zone applied to the pre-incubation arm, +0.0214 AUC on CYP3A4 at 4/4 (item 243). Neither is submittable — the MCC gains are +0.0183 against a floor of 0.0419 and +0.0110 against 0.0281. Those two runs give the exchange rate: 0.021 of AUC buys 0.011 of label MCC, so clearing CYP3A4's floor demands about +0.05 of AUC, more than twice the largest ordering gain anything here has produced. The floor binds, not the modelling.

Two families were closed more than once. The quantum block was measured null against potency (186) and again against the shift that carries the whole remaining classification signal (242). The kinetic scheme arrived three times — reaction phenotyping, bond dissociation energies, a Markov model — and item 211 closes it on identifiability rather than on the tally of eight refutations: two pre-incubation times identify one parameter beyond Ki, and that parameter has now been measured three ways.

## A falsifiable prediction, and what we cannot check

We wrote our score down before it was revealed, and then had to withdraw the first version and
recompute it — which is itself the point of writing it down.

The organisers publish two external samples: a live leaderboard on half the test set, split so that
compounds from one parent stay together, and an interim reveal on the full 750 at the halfway mark.
Resampling the submitted model's out-of-fold predictions 3000 times under a common draw of a
test-sized molecule sample gives a band for each.

| Sample | Predicted macro ST-RAE | sd | 95 per cent band |
|---|---|---|---|
| interim reveal, full test, n = 750 | 0.6537 | 0.0228 | 0.6114 – 0.6997 |
| live leaderboard, n = 375 | 0.6562 | 0.0333 | 0.5945 – 0.7265 |

**Restated 13 September (item 294); the measurement is item 287's.** This table used to carry
0.6666 / 0.6691 with bands 0.6263 – 0.7090 and 0.6126 – 0.7322, which are item 256's figures on the
composition that shipped *then*. Item 287 recomputed the band on 11 September over the regenerated
`oof_submitted.json`, on the composition that ships *now* (items 282–285, plus item 277's lambda
grid); the numbers above are that recomputation, read from `results/preds/band.json`. The shift is
0.013 of macro — the same size as the composition change itself — so the old row was not a rounding
difference but a band belonging to a different arm. What item 294 did was not measure it again but
notice that this table had gone on quoting the superseded one for two days. Item 256's own table is
left standing where it is, as the record of what was true then. One caveat travels with the band,
from item 287: its centre is the seed-0 realisation and seed 0 is inside the set that generated the
composition hypothesis, so the effect size is clean on fresh seeds while this centre is partly
selected.

The band carries two named assumptions, not one. The first is that the test's label spread and band
widths resemble the training set's — that is what makes it a sampling band. The second is newer and
sharper: the submitted predictions are transformed by an affine pair fitted under an *assumed* shift
of the test distribution, and measured against our own labels that bet **costs 0.0131 of macro
ST-RAE** on the arm that ships. It is a bet on a quantity nobody can observe before the reveal.
A score below the band means the bet paid; a score above it is where to look first.

Both halves of that trade were restated on 13 September, and only one of them has a number. The
cost is measured, and it is now 0.0131 rather than the 0.0223 published here before: it is
`band.json`'s own `macro_oof` minus `macro_plain` (0.6506 against 0.6375) on the shipped
composition, and per enzyme it is carried almost entirely by CYP3A4 (+0.0430) and CYP2C9 (+0.0154),
with CYP1A2 (−0.0038) and CYP2D6 (−0.0021) very slightly *against* the tilt — CYP2D6 because item
256 set its delta to zero, so there is nothing left to pay there. The **gain** side is a different
matter: the +0.0473 once quoted beside it comes from `src/shrinkchoice.py`'s posterior for a delta
vector that has since been retired, and it has never been recomputed on this ensemble or this
composition. So it must not be set against the 0.0131 as though the two were commensurable. The
honest statement before the reveal is that the bet's price is 0.0131 and measured, and its expected
return is unquantified on the arm that ships.

The first version of this band was computed on the wrong model. Three configurations were in play:
a four-member ensemble at macro 0.6599, a five-member one at 0.6459, and the one that actually
ships — five members with the dead-zone pass, but the Gaussian process alone on CYP3A4, where a
search over all thirty-one member subsets chose it. That third configuration scored **0.6416**, and
nobody had ever computed it: not the journal, not the saved predictions. The published band was
0.0186 too high.

What settles the arm's identity now, and was missing before, is an anchor: CYP3A4 comes out at
0.4092 against 0.4129 from the independent subset search for the Gaussian process alone on that
enzyme. The cell lands where the selection decision puts it rather than merely looking plausible.
The failure mode is worth naming — a saved arm carried the right *name* and the wrong *composition*,
and the check that catches that is reconciling an arm's own score against the row for the
configuration being claimed, not against any row showing the same figure.

This is a sampling band and not a prediction interval: it covers the draw, not the distribution shift, because it assumes the test's label spread and confidence-band widths resemble the training set's. A score outside it falsifies that named assumption, not the model; a score inside it confirms nothing.

The band was tightened after it was first written. Item 147 built the original by averaging four per-enzyme percentile bounds as though the enzymes' sampling errors moved together. They do not — the cross-enzyme correlation of per-draw ST-RAE under a common draw is 0.02 — and averaging four nearly independent errors halves the spread. The half-width fell from 0.082 to 0.042, so the corrected test is stricter than the one it replaces. A band widened after the fact is worthless; one narrowed is defensible only if the narrowing is announced with its reason before the measurement, which is what item 246 does.

We pre-register because three independent routes each found the test regime uncheckable from inside. Importance reweighting is capped at an effective sample of 26.1 per cent by a chi-square divergence of 2.838 between test and out-of-fold similarity (item 123). Building a test-like split is blocked by composition: 93.6 per cent of Butina clusters are singletons, and the training set is a 4375-compound diversity screen glued to a 530-compound CYP3A4-only campaign (item 129). And a rigorous chi-square bound on the score is vacuous: ST-RAE's denominator, a constant predictor at the mean through the same soft threshold, is zero for a large share of compounds, so its standard deviation exceeds its mean and the bound admits a reweighting that drives it to zero (item 147).

One number calibrates the rest. The denominator belongs to the test set, not to the submission, so it cancels in the difference between two submissions scored on the same test set; that paired difference is the only quantity a leaderboard measures cleanly. Run on the paired difference rather than a single score, the same bootstrap gives sd 0.0144 at n = 750 and 0.0203 at n = 375 for our per-enzyme arm against our pooled arm, against 0.0220–0.0231 and 0.0306–0.0322 for a single score. Our whole trajectory from the baseline is 0.0583 of post-pair ST-RAE: about four sigma paired — taking the per-enzyme against pooled pair as the proxy for paired spread, since the paired sd of the baseline-against-submitted comparison was never measured directly — and under three unpaired.

Four things stay unchecked. The shift itself; nothing above measures it. The similarity geometry: the test sits at median nearest-neighbour similarity 0.587 to the training set, against 0.435 for our cluster split and 0.450 for leave-one-out, the most generous re-slicing available, so none reaches the test (item 22). In the one place a test-like split is constructible — the CYP3A4 analogue campaign, held-out similarity 0.61 to 0.62 — the mechanistic block's gain fell from +0.0056 of rank to +0.0018, both inside the floor (item 155). Pooling, worth +0.0141 of rank over four seeds against a macro floor of 0.0036, has a mechanism only by elimination: neighbour-borrowing, shared-function transfer, sample size and enzyme level were each refuted by measurement (items 110, 111, 125, 132), leaving contrast — and whether that enzyme-conditioning structure survives the shift is exactly what nothing internal can test. And the classification track, a third of the leaderboard, has a split-seed floor of 0.0076 macro MCC (item 235) and no sampling band at all: the bootstrap above was written for ST-RAE and has never been run for MCC. Internal agreement settles none of it: a head predicting the score itself matched the submitted arm's macro to within 0.0006, which item 247 records as a tautology.

---

## Reproducing this

Everything runs through `uv run`, which pins the interpreter and the library versions.

```
make setup                          # once per machine
uv run python src/feats.py          # data/feats.npz, data/rows.csv
uv run python src/ablate.py         # ~1 h, rewrites results/preds/oof.json
uv run python src/score.py          # metrics and paired bootstrap
uv run python src/submit.py         # the two submission files
uv run pytest                       # golden-value guard on the split
```

The cross-validation split is defined once, in `cypsplit.py`, and pinned to a golden digest that CI
checks on every push: Butina clustering over Morgan fingerprints at Tanimoto 0.35, each cluster
assigned whole to one of five folds, 4703 clusters over 4905 molecules. Numbers from different
scripts are comparable only while their folds agree element for element, so the digest is the
guarantee that the tables in this document describe the code in the repository.

The environment pins are measured rather than cautious: scikit-learn 1.3.2 through 1.8.0 regenerate
`results/preds/oof.json` bit for bit, and 1.9.0 does not — it moves every prediction by about 0.135
pIC50 and doubles the measured effect of the mechanistic block. numpy is deliberately unpinned;
1.26.4 and 2.5.2 both reproduce exactly.

Run logs are results here and are committed under `results/logs/`. Saved out-of-fold predictions are
under `results/preds/`, and every table in this document can be recomputed from them without
refitting a model.

No proprietary data is used. The external sources we tried — a public NCGC panel, ChEMBL, and three
pretrained molecular representations — are all public, and all of them are closed negative in the
journal.
