# Method report — OpenADMET CYP Inhibition Blind Challenge

*Outward-facing document, linked from the leaderboard's Model Report field. English by intent: the
rest of `docs/` is the team's internal Russian-language artefact, this one is written for the
organisers and for anyone who wants to reuse the negative results.*

Code: this repository, public. Proprietary training data: none. External data: the organisers' own
files only — the four dose-response tables and the single-concentration screen.

---

## 1. What we submit

Two files, one per track.

**Direct inhibition.** A five-member ensemble over a shared feature matrix (2048 Morgan counts at
r=2, 217 RDKit descriptors, 30 mechanistic columns), with the member set chosen **per enzyme** by a
search over all thirty-one subsets, and a dead-zone refit applied inside every member. Members:
per-enzyme gradient boosting, a contrast-pooled model, a Gaussian process, a neural trunk, and a
re-projection of that trunk. The per-enzyme compositions that ship are

    CYP1A2   per-enzyme + GP + trunk
    CYP2C9   GP + trunk
    CYP2D6   all five members
    CYP3A4   GP + trunk

**Time-dependent inhibition.** Not a classifier on the label. The shipped rule builds the label's
two conjuncts separately — a gate probability from the pre-incubation arm predicted by a four-member
ensemble, times a shift probability from a classifier on Δ > log₁₀2 — and thresholds the calibrated
product.

Out-of-fold performance, Butina-clustered five-fold split, four seeds:

    трек            метрика                 значение
    регрессия       macro ST-RAE            0.6459
    регрессия       macro Spearman          0.6342
    TDI             MCC CYP3A4              0.3510
    TDI             MCC CYP2D6              0.1235
    TDI             macro MCC               0.2373

From a conventional starting point (same features, per-enzyme boosting, no ensemble) at 0.7150 pair
and 0.5651 rank, the trajectory is **+0.0579 of rank and −0.0583 of ST-RAE**.

We will say the uncomfortable part first: **the architecture is the least novel thing in this
report.** Gradient boosting, a Gaussian process and a small trunk are 2015 technology. What we think
is worth reading is where the gains came from, and — more of the document — where they did not.

## 2. The finding we would most want cited

**Physics enters this problem through the measurement model and returns zero through the feature
matrix.**

That sentence is not a slogan; it is what the log says once the interventions are sorted into two
piles. Everything we tried that modelled *how the assay produces a number* paid. Everything we tried
that added a *descriptor of the molecule* did not.

    вошло в подачу (модель измерения)                 прирост ранга   знак
    мёртвая зона во всех членах                            +0.0197    4/4
    механистический блок                                   +0.0163    4/4
    пулирование контрастом                                 +0.0141    4/4
    ствол пятым членом                                     +0.0054    4/4
    уравнение прибора (скрининговая голова)                +0.0308    4/4

    отвергнуто (матрица признаков)              над ансамблем   пол   как член   пункт
    3D-форма: PBF, NPR, асферичность, PMI (16)        +0.0004  0.0049   +0.0087     281
    рукотворные блоки активного центра, 2D6           −0.0094  0.0049         —     179
    перекрытие с со-кристальной позой (контраст)      +0.0017  0.0061         —     296
    внешние CYP-панели, вспомогательная голова        −0.0104      —    +0.0085     290
    мультитаск-деревья                                −0.0062      —    +0.0059     292
    SMARTCyp, сайты метаболизма (16)                  +0.0009      —         —     179
    квантовый блок: HOMO/LUMO/щель/Фукуи (10)             ~0       —         —     186
    CYP2C19 пятой изоформой как донор                 −0.0082      —         —     153

The right-hand columns are the point. Three of these arms carry a **real** signal as a standalone
member — the shape block +0.0087 on CYP2D6, the external-panel channel +0.0085 at sign 4/4, the
multi-task channel +0.0059 at sign 4/4 — and every one of them arrives at or below its floor once it
has to share an ensemble with five members that already correlate 0.91–0.99 with each other.

The dead zone is the largest single effect in the project and the clearest example of the law. The
organisers' metric forgives any prediction inside the credible interval, so the loss a model should
minimise is not the distance to the point estimate — it is the distance to the nearest bound. We
refit every member by majorise–minimise against `clip(p, lo, hi)`. That is not a feature and not an
architecture; it is reading the metric as a statement about the instrument.

The instrument equation is the same move applied to the screening data. The single-concentration
screen and the dose-response tables measure the same enzymes under different conditions; writing
down the relation between them, rather than concatenating one onto the other, gave +0.0308 of macro
rank at 4/4 seeds, t = +18.61, p = 0.0003.

Contrast pooling is the third. Pooling enzymes helps — but we spent five items on *why*, and the
answer is neither borrowed neighbours nor shared function nor sample size. Zeroing the enzyme
indicator drops the model 0.076 of rank **below** training each enzyme separately, and handing the
model each enzyme's mean for free recovers none of the gain. It works by **contrast**, and the
contrast is purchased with residual magnitude: switching the pooled member to absolute error
destroys 95 per cent of it, and the effect is **entirely CYP2D6** — the enzyme whose pooled member
does the work is also the only one that loses it.

## 3. Only rank survives, and this is a warning to the field

Before writing predictions out we fit an affine pair — a per-fold shrink and shift — to the metric.
It is strictly increasing, so it preserves Spearman exactly and rewrites everything else.

We recommend anyone reading this check their own gains against that transform, because we did it
late. Five separately measured improvements in our own log collapse into a band of 0.013 once the
pair runs — none clears the 0.007 noise floor and **three of the five change sign to harmful.** The
largest finding of that day was worth 0.046 of raw ST-RAE and 0.002 after.

The sharpest single case was a concatenated Chemprop embedding. Raw macro ST-RAE improved by 0.0084
(0.7673 → 0.7589). After the affine pair it **lost** 0.0113 (0.7150 → 0.7263). By rank it lost
0.0075, and under the test set's regime it lost 0.0182 — more than twice as much. The raw number and
the post-processed number carry opposite signs, and we published the raw one first. A raw
improvement in this challenge is not evidence of anything until it has been through the
post-processing that ships.

The mechanism is worth stating because it generalises: the affine pair is a shrinkage fitted
directly to the metric, per fold, so anything that merely reduces prediction variance or removes a
systematic offset is something the pair already does. It does not add to such improvements — it
substitutes for them.

The practical consequence for the leaderboard's secondary metrics: the organisers publish Spearman
ρ and Kendall τ with bootstrap intervals, and those are the columns an affine transform cannot
touch. We optimised them on purpose.

## 4. Uncertainty, over conclusions rather than over predictions

We do not ship predictive intervals, and we say so rather than dress up what we have. What we do
have is a quantified answer to "is this number real", used as a gate on every claim.

    поле шума (четыре сида)                        значение
    макро по рангу                                   0.0036
    CYP1A2 / CYP2C9 / CYP2D6 / CYP3A4     0.0061 / 0.0071 / 0.0049 / 0.0033
    парный (разность двух подач)                     0.0052
    MCC: CYP3A4 / CYP2D6 / макро          0.0281 / 0.0419 / 0.0076

Nothing below its floor is reported as a result. The split seed alone moves macro ST-RAE by 0.016 —
more than a typical effect — which is why every claim is quoted over four seeds and why the sign
count (how many of four seeds agree) is printed beside every mean.

We also pre-registered what the reveal should show, before it happens. On the arm that ships:
macro ST-RAE inside **[0.6114, 0.6997]** at n = 750 (half-width 0.0442) and **[0.5945, 0.7265]** at
n = 375; macro Spearman half-width **0.0157** at n = 750 and **0.0270** at n = 375. The n = 375 rows
are stated as *lower bounds on the spread*, because the live half is split by chemical series and a
clustered split has strictly higher variance than the random resampling those intervals assume.

One consequence we want on the record because it cuts against us: our absolute score is
unpredictable to about ±0.04, while our position relative to a *similar* submission is pinned to
about ±0.02. Anyone comparing two entries a few thousandths apart on this leaderboard is reading
noise, ours included.

## 5. How claims were tested

Every intervention in this project went through the same harness, and the harness is the part we
would defend as the contribution.

- **Pre-registration.** The acceptance threshold and a written prediction are committed **before**
  the run. When the prediction is wrong, the item says so — including the cases where we were right
  about the outcome and wrong about the mechanism.
- **Wrong-isoform controls, not permutation.** For any quantity derived from a (ligand, cavity) pair,
  the control is the same quantity computed against a *different* enzyme's cavity. A permutation
  control cannot separate real pocket information from a disguised volume descriptor; the wrong-cavity
  control can. On CYP3A4 our cheap structural proxy correlates +0.283 with heavy-atom count even
  after centring, and a permutation test never noticed.
- **Harness controls that must read exactly zero.** Where an arm rebuilds a component that is kept on
  only some enzymes, the untouched enzymes are required to print `0.000000` to every digit. Twice
  this caught nothing and confirmed the rig; that is what it is for.
- **Fresh seeds for confirmation.** A composition chosen on seeds 0–3 is confirmed on seeds 4–7, and
  the item states which seeds generated the hypothesis.
- **Knockout, not addition.** Named contributions sum to +0.0555 when each is measured by *adding* it
  to its own base, against a measured total of +0.0579. Measured by *removing* each from the shipped
  configuration they sum to **+0.0315** — and contrast pooling, worth +0.0141 by addition, is worth
  **+0.0008 at sign 0/4** as a marginal member. The table of contributions is not additive and we
  stopped presenting it as though it were.

The journal is public and runs to item 302. Most of it records something that did not work, and we
think that is the useful part — a correct idea takes one item, a wrong one takes three (proposal,
refutation, correction to the refutation).

Item numbers in this report are checkable: the affine-pair collapse is 77 and its by-rank re-reading
is 117; pooling's mechanism is 131, 132 and 264; the shape block is 281; the estimator-family screen
is 289; transfer and multi-task are 290 and 292; the structural proxy is 296; and the pre-registered
reveal bands are 294.

## 6. What we refused, and what that closed

Five conclusions we would hand to anyone continuing this problem.

1. **Descriptors of the enzyme are closed by arithmetic, not by measurement.** Four enzymes, all four
   seen in training, no zero-shot requirement — so a one-hot indicator is a sufficient statistic for
   enzyme identity, and any vector of protein descriptors is a four-row table, i.e. a change of
   basis. Measured: +0.0007. A fifth isoform as donor: −0.0082. What works instead is the
   *interaction* form: ligand columns conditioned on a known site.
2. **Ensemble redundancy is the binding constraint, not accuracy.** Member error correlations run
   0.91–0.99 across every pair (item 280). Four separate arms, measured in four separate items, each
   produced a genuine single-model gain that died on insertion: the largest was the shape block's
   **+0.0087 on CYP2D6 diluting to +0.0004** at sign 2/4 against that enzyme's floor of 0.0049
   (item 281), and CYP2D6's per-enzyme member effect of **+0.0149** at 4/4 — three times its floor —
   reached **+0.0011** over the ensemble. Item 269 sharpened the mechanism and is the one to read:
   cross-enzyme stacking arrived at **−0.0007** over the ensemble and died of **redundancy, not
   dilution** — the other members already carried the information. Escaping dilution is therefore
   necessary and not sufficient, which is the sentence we would put on the wall.
3. **No estimator family is both accurate and decorrelated** (item 289). Random forest and
   extra-trees decorrelate to only 0.95–0.96; SVR-RBF reaches 0.898 and kernel ridge 0.661, and both
   fail on accuracy instead. Sweeping the kernel width does not find a point that is both: kernel
   ridge peaks at macro rank 0.5417 and never reaches the weakest shipped member at 0.5776.
4. **The (ligand, cavity) pair-function class is not empty, but its magnitude is small.** Our targeted
   structural proxy beat its wrong-isoform control on CYP1A2 at sign 4/4 — the discrimination a
   permutation test cannot supply — and still came in at +0.0017 against that enzyme's floor of
   0.0061, from a univariate precondition of +0.173. A roughly hundred-fold collapse between a raw
   correlation with the label and a rank gain over a three-member mean.
5. **Transfer and multi-task both fail the same way.** An auxiliary-head trunk trained on external CYP
   panels and multi-task trees each showed a real channel over the ensemble (+0.0085 and +0.0059, 4/4)
   and each lost by rank once inserted (−0.0104 and −0.0062). The channel is real; the vehicle is
   redundant.

Several of these were reached independently by other entrants whose code is public, which we take as
corroboration rather than coincidence: 3D shape descriptors net-negative for an ensemble, tabular
models on ECFP4 correlating r ≈ 0.9 with a gradient-boosting baseline, and pseudo-labelling from the
single-concentration screen degrading accuracy monotonically with pseudo-weight.

## 7. Reproducibility

- Environment pinned by `uv.lock`; the `scikit-learn` ceiling is measured, not cautious — 1.3.2
  through 1.8.0 regenerate the saved out-of-fold predictions bit for bit, while 1.9.0 moves every
  prediction by ~0.135 pIC50.
- The split lives in one module and is pinned to a golden digest (`2d93c19815e14261`, 4703 clusters)
  by a test that CI runs on every push and every pull request. If the digest moves, the published
  tables no longer describe the code.
- CI additionally re-scores the committed predictions against documented values and imports every
  script, so a file that stops running is caught rather than discovered later.
- Structural inputs that are too large to commit are rebuilt by committed scripts, both default-safe:
  `src/prep_receptors.py` (chain A plus its own heme, parsed by column) and `src/prep_lig3d.py`
  (ETKDGv3, fixed seed, MMFF cleanup, names carrying the row index).
- Four people work from four machines; the protocol for the unmergeable prediction files, and the
  tripwire that verifies every commit on `main` arrived through a pull request with a green split
  guard, are in `CONTRIBUTING.md`.

## 8. Limitations

- **The architecture is conventional.** If the award weighs architectural novelty, this entry does not
  compete on that axis and we are not going to pretend otherwise.
- **No predictive uncertainty.** Our uncertainty work is about our own conclusions, not about
  individual predictions. We consider that an honest reframing, not a substitute.
- **One enzyme carries the rank.** CYP3A4 reaches Spearman 0.822 out of fold while CYP2D6 reaches
  0.480. Part of that gap is metric structure rather than model quality — CYP2D6 has the narrowest
  credible intervals of the four, so the soft threshold forgives least there, and CYP3A4 has two
  fifths of its labels below pIC50 4, where the metric downweights.
- **CYP2D6's TDI cell is weak** (MCC 0.1235 against a floor of 0.0419) and we have not fixed it. Public
  entries suggest it is hard for everyone, which is an explanation and not an excuse.
- **A campaign is still running.** Docking 5652 ligands into four cavities was pre-registered blind,
  with a written prediction that nothing passes. Whatever it returns will be in the journal.

## 9. A note on metric robustness

Offered as an observation about the scoring function, not as a complaint about anyone.

ST-RAE clips each error to the credible interval and divides by the same quantity for a
mean-predictor baseline. Both halves therefore depend on where the test labels sit and how wide their
intervals are, and neither depends on ordering. Across our four enzymes, out-of-fold ST-RAE runs
almost inverse to median interval width. The practical consequence is that **placement** — the mean
and spread a submission is centred on — moves the primary metric a long way independently of whether
the predictions are correctly ordered, and the leaderboard displays rows where ST-RAE is good while
Spearman is approximately zero.

We treated the test distribution as unknown and estimated it the only way that touches no leaderboard
feedback: nearest-training-neighbour label transfer for the 750 test compounds, with the estimator's
bias calibrated leave-one-out on the training set, where the true moments are known. The gate passed
(mean shifts within 0.062, sd ratios 0.948–1.021), and the estimate puts the blind set at or above
training potency on all four enzymes — which is what an analog expansion of the top hits implies by
construction. We pre-registered, before running it, that a disagreement with any other estimate would
be reported as *ambiguous* rather than as a refutation, because the live half is series-split and
need not represent the full 750. It disagreed, and we are holding to that.

Nothing in our submission is calibrated onto any estimate of the blind label distribution. If the
organisers find the placement sensitivity worth addressing in future challenges, a secondary metric
that is invariant to a monotone transform — Spearman, which they already publish — does the job.
