# Verification

Twenty-five scripts in four groups. `f*` was a sweep over everything that had been computed
and written by that point; `g*` answers four questions raised against the document; `h*`
tests two claims the document made about geometry and about reactivity; `k*` is about
post-hoc rescaling of the predictions and about how far the test set sits from the training
distribution. All of them expect the data in `data/` and the organisers' metric code
(`git submodule update --init`).

The only dependency between scripts is that `k6_shift1d.py` reads `data/test_pred.npz`,
which `k5_shift.py` writes. Everything else runs in any order.

| Script | What it checks | Time |
|---|---|---|
| `f1_formulas.py` | 44 numerical checks of every formula in the document against an independent implementation | ~40 s |
| `f2_splitnorm.py` | accuracy of the rule `sigma=(bound-mu)/1.96` for a split normal | ~30 s |
| `f3_seeds.py` | ST-RAE at split seeds 1, 2, 3 (seed 0 is already computed) | ~60 min |
| `f4_macro.py` | paired bootstrap at the macro level, 3000 resamples | ~3 min |
| `f5_rho.py` | paired bootstrap of Spearman per enzyme | ~4 min |
| `f6_data.py` | 21 numbers from the document reconciled against the source tables | ~20 s |
| `f7_testset.py` | structure of the test set and the percentiles of its anchors | ~5 min |
| `f8_mnar.py` | simulation of the MNAR to MAR transition, 40 repeats | ~20 s |
| `f9_mccthr.py` | cost of a plug-in threshold against a fitted one | ~30 s |
| `f10_calib.py` | whether calibration repairs the plug-in threshold | ~30 s |
| `f11_pka.py` | prevalence of the motif the pKa rule gets wrong | ~2 min |
| `f12_cvhard.py` | difficulty of our cross-validation against the real test set | ~10 min |
| `g1_calib.py` | instrument calibration: fitting E and h per enzyme, plate drift | ~1 min |
| `g2_factor.py` | whether the joint distribution factors into three marginals | ~1 min |
| `h1_geometry.py` | where the test set sits on the random-split / cluster-split axis | ~3 min |
| `h2_tdi_alerts.py` | structural reactivity alerts against the TDI label | ~20 s |
| `h3_alerts_delta.py` | the same alerts against the shift Delta rather than the label | ~20 s |
| `k1_shrink.py` | shrinkage to the mean: real effect or an artefact of the metric | ~2 min |
| `k3_center.py` | what the shrinkage centre actually is and where its gain comes from | ~4 min |
| `k4_enrich.py` | is the test set activity-enriched — nearest-neighbour-label proxy | ~2 min |
| `k5_shift.py` | covariate shift measured in the model's own prediction space | ~5 min |
| `k6_shift1d.py` | the same reweighting done along one axis, where it does not degenerate | ~1 min |
| `k7_2d6shift.py` | the CYP2D6 shift against basic-amine composition, and the raw shifts | ~2 min |
| `k8_kernel.py` | delta by inverting the kernel E[yhat\|y] — no instrument, no scalar propagation | ~5 min |
| `k9_shape.py` | the test is shifted *and* widened, and what that does to ST-RAE | ~2 min |

## What it found

Eight discrepancies between the document and the data, all corrected:

1. "75 series of 10" in the test set — clustering gives 194 groups with a median size of two.
2. Anchor percentiles 92.5 / 98.1 / 97.2 / 68.9 do not reproduce; the correct ones are 93.0 / 97.5 / 98.2 / 62.8.
3. 0.236 was called the "typical error" — that is the root mean square; the median is 0.069.
4. "92 % of the ranking ability of the full model" — the full model yields 80 %.
5. The 21.3 % positive TDI rate refers to a different population than the one MCC was computed on.
6. The plug-in threshold by expected MCC requires a calibration that was not mentioned.
7. Changing the split seed moves macro ST-RAE by more than the entire effect being measured.
8. The screening table was assembled from three runs with different baselines.

Items 1 and 2 have since been superseded by item 12 below and are kept for the record
rather than as current findings: the "194 groups with a median size of two" holds at a
Butina threshold of about 0.48 and not at the 0.35 this repository works at, and of the
four replacement percentiles the CYP3A4 one, 98.2, does not reproduce at any threshold.

Four claims were established for the first time: the MNAR to MAR transition was confirmed
by simulation; the Spearman gain from the mechanistic block was confirmed by bootstrap;
the Bayes optimum under ST-RAE was verified numerically; and our cross-validation was
measured to be 0.153 harder than the test set.

No leakage was found in the pipeline. There is no contamination between training and test:
the overlap by name, by canonical SMILES and by SMILES with stereochemistry stripped is
zero across all three training files.

## Thirteen more, found while setting up the environment

Mostly documentation errors — the code was right and the prose was wrong. Where a script
was at fault it is said so explicitly.

9. The mechanistic block has **30** features, not 28. `data/mech_names.csv`, written by
   `feats.py`, has always listed 30 and matches the copy committed in `results/`.
10. There are **fourteen** verification scripts, not twenty-three, in both this file and
    the top-level README. (Eight more were added later, in the `h*` and `k*` groups; the
    current count is twenty-five and both files say so.)
11. `f12_cvhard.py` could never have run to completion as committed: line 34 referenced an
    undefined `fRES`. Fixed, and the script now reproduces the documented -0.153.
12. The test-set clustering figures in §2 (194 groups, median size two, 55 groups of five
    or more) come from a Butina threshold of about 0.48. The repository works at 0.35,
    where the same computation gives 477 groups, median size one and 21 groups of five or
    more. `f7_testset.py` now sweeps the threshold so the dependence is visible: the
    anchor percentiles turn out to be stable across it (1A2 92–94, 2C9 97.5–99, 2D6 58–66)
    and only the 3A4 figure moves, 82–95 against the 98.2 the document quoted. The number
    of anchors carrying a label is nine to thirty-seven per enzyme, which the text did not
    say and which is what actually limits the claim.
13. §4's plate numbers had no source in the repository. `g1_calib.py` printed only the raw
    spread; the residual analysis the section quotes was computed nowhere. It now is, and
    the answer differs: plate explains 0.9–8 % of residual variance, not 1–3 %, with
    CYP2C9 at 8.1 % — nine times CYP3A4's 0.9 %. Subtracting with a free monotone
    calibration instead of the Hill form gives 0.9–6.7 %, so this is a property of the
    data rather than of the calibration. The raw spread was also understated: 0.293 and
    0.922 against the quoted 0.24 and 0.78.

    This weakens, without overturning, the one risk row §14 records as retired. Plate
    drift still needs no term in the model — even on 2C9 the plate component is about
    0.05 in log2fc units, well inside a single reading's own error — but the claim now
    rests on a measurement of the residual rather than of the raw spread.

14. §10 justified the dead ST-RAE decision layer by saying actives determine the metric
    while inactives contribute almost nothing. The first half holds. Decomposed by
    activity zone, inactives carry 22 / 16 / 12 / 19 % of the numerator and 34 / 30 / 20 /
    25 % of the denominator, their bands are visibly asymmetric (median a 0.997 against b
    0.702 on CYP1A2), and on CYP3A4 the intermediate zone carries more of the numerator
    (43 %) than the actives (38 %). The conclusion rests on the measurement, which is
    unaffected; the mechanism connecting it to a cause is now an admitted gap.
15. §12 claimed random five-fold is about twice as optimistic as the cluster split in R².
    It is not: 0.314 against 0.320 on average, and stricter on three enzymes of four. At
    threshold 0.35 there are 4703 clusters for 4905 compounds, so the cluster split is
    almost all singletons. Its effect is in the tail, not the median — the share of
    held-out compounds keeping a close relative in training falls from 0.033 to 0.005.
16. §12's third level of split strictness, the pseudo-test fold, cannot exist as a split.
    Leave-one-out — the most generous partition there is — gives median nearest-neighbour
    similarity 0.450, and the test sits at 0.587, above that ceiling. As a *stratum* it is
    constructible (328 compounds, median 0.643) but lands 267 of them on CYP3A4 against
    38 / 36 / 32 elsewhere, so it is a real check on 3A4 and not a check at all on the
    other three.
17. §14's summary table was captioned "all of it is five-fold Butina cross-validation".
    Six of its twelve rows are not: label counts and interquartile range are column
    statistics, reliability is a variance decomposition, and calibration E, h and the
    screening rank correlation are fits over the whole sample.
18. §14's tripwire for the series prior named pseudo-test folds, which §12 had just shown
    cannot be built; and its yardstick for "local estimate disagrees with the leaderboard"
    was the gap between random and cluster splits, now measured at essentially nil, so it
    would have fired on anything. Both replaced with measured thresholds.
19. `f11_pka.py` had never run either: two path substitutions were left as the literal
    string `'" + D + "'`, `src/` was missing from `sys.path`, and it read the atom index
    where it wanted the class name, which printed two of its four numbers as 0.0 %.
    Repaired, it reports that the most basic centre is an aromatic nitrogen in 55.8 % of
    the training set and that the caffeine configuration covers 2.2 %.
20. §13's "the rank always improves" is true of the macro only. Two-stage inversion drops
    Spearman on CYP2D6 from 0.398 to 0.362 and on CYP3A4 from 0.749 to 0.747; only the
    blend lifts all four.
21. §11's median |Δ pIC50| over pairs at Tanimoto 0.6–0.85 is 0.33–0.68, not 0.24–0.55,
    and three of the four medians rest on 15 to 30 pairs.

Smaller ones, all corrected in the same pass: the TDI fill rate is 21–58 % and not 24–58;
the "curves were run on whatever the screen flagged" rule holds for three enzymes, since
all 530 curve compounds without a screening reading carry a CYP3A4 label and nothing else;
CYP2D6 sits 0.007 below the trivial line, not two hundredths; CYP3A4's E spread is 0.067
and CYP2D6's 0.091; `f8_mnar` selects the top 28 % and not thirty; "understates its own
error by 17–21 percent" was the ratio 1.18–1.21 read as a percentage; the TDI rule was
checked on 2334 of CYP3A4's 3584 labels, the other 1250 having no direct arm and all being
negative; "hundreds of false positives" is zero with the guard and 56 on CYP2D6 without it;
the MCC figure caption and body gave 0.356 and 0.358 for the same quantity, a threshold-grid
artefact whose exact value is 0.3574; +0.028 MCC is the top of a +0.018–0.029 range over
calibration seeds and is not commensurable with the ST-RAE feature ablations; and the 0.015
decision-layer ceiling subtracts two rows computed with different Monte-Carlo settings.

## Scripts that had never run

Five, all broken by the migration off the original machine rather than by anything in the
analysis: `feats.py` (hard-coded dimorphite path), `f12_cvhard.py` (`fRES`), `f11_pka.py`
(three separate faults), `docs/tex/figs.py` (`/tmp` paths, then a name collision with a
pivot table), and `g1_calib.py`, which worked only because an import inside a loop body
happened to leak. All five run now, and the numbers they produce are in the document.


## The `h*` and `k*` groups

Written after the `f*`/`g*` sweep, against two questions the document had left open — how
far the test set really is from our cross-validation, and whether post-hoc rescaling of the
predictions is worth doing — plus one idea from a chemist about reactivity. The numbers
below were produced under scikit-learn 1.8.0 and pandas 3.0.2, the first inside the range
`pyproject.toml` documents as bit-for-bit reproducing and the second outside the `<3` pin,
so they want one confirming run under the lock: `make verify-extra`.

**22. The third level of split strictness cannot be built, and the reason is not that our
split is too easy.** `h1_geometry.py`: leave-one-out — every compound against the entire
rest of the training set, the most generous partition there is — gives median
nearest-neighbour similarity 0.450; random five-fold gives 0.437 to 0.438 over four seeds;
the cluster split gives 0.435; and the test set sits at 0.587. The test is closer to the
training set than the training set is to itself, so no re-slicing reaches it. The cluster
split earns its keep in the tail rather than the median: the share of held-out compounds
keeping a close relative (>0.7) in training is 0.005 against 0.033 for random and 0.101 for
the test. This is the measurement behind item 16.

**23. Structural reactivity alerts say nothing about the TDI label and a great deal about
the shift behind it.** `h2_tdi_alerts.py` scores fourteen classical mechanism-based
inactivation motifs against `is_TDI` and finds nothing — the largest absolute MCC on CYP3A4
is about 0.013, which is noise. `h3_alerts_delta.py` scores the same motifs against
Delta = pIC50(TDI) - pIC50(direct), conditioned on the compound being potent at all
(pi > 4), and finds a strong signal: on CYP3A4 cyclopropylamine gives a median Delta of
0.574 against 0.289 without it (p = 0.0007) and any-alert 0.369 against 0.277
(p < 0.0001); on CYP2D6 methylenedioxyphenyl gives 0.432 against 0.139 (p < 0.0001).

The two results are consistent, and their difference is the point: the label is
"potent AND shifted", the alerts are about the second half only, and the threshold on the
second half is a knife edge. The median Delta over CYP3A4 actives is +0.294 against a
cutoff of log10(2) = 0.301. Half the actives sit within a hundredth of the line, so a
predictor of Delta can be good while a predictor of the label looks worthless.

**24. Shrinkage toward the mean is a real effect measured by the wrong instrument.**
`k1_shrink.py`. Out-of-fold least squares independently recovers the attenuation —
`a = (1-b)*mu` to three decimals, with b = 0.789 / 0.898 / 0.743 / 0.999 — so the model
really is over-dispersed. But the gain is almost entirely ST-RAE's: MAE on CYP3A4 gets
*worse* (0.5348 to 0.5409) and R^2 does not move. On the top quartile by activity shrinkage
worsens every enzyme. It is off by default in `src/submit.py`.

**25. The "shrinkage centre" is not a shift of the predictions, and the gain does not come
from where it looked like it came from.** `k3_center.py`. Since
`c + l*(p - c) = [mu + l*(p - mu)] + (1-l)*off`, the predictions move by `(1-l)*off`, which
is +0.136 / +0.102 / +0.024 / +0.173 for the fitted optima — the centre offset of about
+0.40 quoted elsewhere is the same thing inflated about fivefold by the parameterisation.
The family also degenerates: `(c, l)` maps onto `a + b*p` with `a = c*(1-l)`, so `c` is
unidentifiable at l = 1 and unstable near it.

The obvious mechanism — the offset buys a free lunch inside the wide bands of censored
inactives — is wrong. The share of predictions falling inside the band *drops* when the
offset is applied (CYP3A4: 37.8 % to 34.5 %). Decomposed by activity zone, the numerator
change on CYP3A4 is -62.3 in the intermediate zone and -11.8 on the actives against +35.8
on the inactives; the same signs hold on all four enzymes. The offset corrects the
intercept in the middle of the range and pays for it with the inactives.

Macro out of fold: raw 0.7673, shrink to mu 0.7330, shrink to mu + off 0.7149. On an
activity-enriched evaluation the *oracle* affine map — fitted on that very subset — only
reaches 0.90 to 1.00, which is the constant-mean baseline. Within the active range our
predictions carry almost no ST-RAE-visible skill.

**26. The test set is shifted, and by much less than the stress test assumed.**
`k4_enrich.py` takes the label of each compound's nearest training neighbour as a proxy and
reports a large enrichment (+0.45 / +0.83 / -0.02 / +1.16). The proxy is confounded and the
script says so: test compounds have higher nearest-neighbour similarity than training
compounds have to each other, so their neighbours are not comparable draws. Only the
direction survives that objection — noise pulls a neighbour's label toward the mean, not
away from it.

`k5_shift.py` does it properly, applying one model to both sets and comparing the
distributions of its own output. The shift is +0.014 / +0.147 / -0.190 / +0.437, the
standard deviation rises on all four, the share of predictions above 5.5 rises on three
(CYP3A4 2.1 % to 14.5 %), and every KS p-value is at most 7e-3. CYP2D6 goes the other way,
as it has in every other diagnostic in this repository.

**27. Correcting for the shift I measured changes nothing — but I measured the wrong
marginal, and `src/reweight.py` disagrees for a good reason.** `k6_shift1d.py` reweights
the training set by the density ratio along the *prediction* axis (effective sample size
796 to 1444). The optimal `(off, l)` barely move — +0.95 to +1.10 and 0.82 to 0.84 on
CYP3A4 — and the ordering of strategies does not move at all. Read on its own that retires
the "shrinkage does not transfer" claim, since `k1`'s top-quartile stress test moves the
evaluation mean by +1.1 to +1.3 while the shift measured here is far smaller.

Read against `src/reweight.py` it does not, and the disagreement is the useful part. That
script tilts the *label* marginal by delta and finds the optimal centre tracking delta
closely: +0.4 at delta = 0, +0.8 at delta = 0.3, +0.9 at delta = 0.5, with shrinkage at the
training mean ceasing to help around delta = 0.5. Two reweightings, opposite conclusions,
and the reconciliation is attenuation. The model is over-dispersed by b = 0.789 / 0.898 /
0.743 / 0.999 (item 24), so a label shift of delta shows up in prediction space as about
b * delta. Inverting the measured prediction shifts gives an implied delta of +0.018 /
+0.164 / -0.256 / +0.437, macro +0.09 — the very bottom of the range `reweight.py` scans,
where its own table also says the optimum has not moved yet.

So the two agree on the mechanism and differ on where the test sits, and neither pins it
down. b measured in-distribution is itself an upper bound on how much an out-of-distribution
shift propagates into the predictions — a model that mostly interpolates would show almost
none of it — so +0.09 is a lower bound on delta, and `reweight.py`'s anchor-percentile
argument gives +1.05 as an upper bound and +0.3 to +0.6 as defensible. The honest statement
is that delta is bracketed between roughly +0.1 and +0.6, that the centre should be near
+0.4 at the bottom of that bracket and near +0.8 at the top, and that nothing available
before the intermediate leaderboard narrows it further.

Both corrections share one assumption neither can check: that p(y | yhat) is unchanged on
the test set. That is exactly what recalibration is supposed to test, so this is not a free
lunch either way.

**28. Two independent orderings of the enzymes coincide.** The strength with which the
screening channel tracks pIC50 (0.936 / 0.896 / 0.862 / 0.828) ranks the enzymes exactly as
the shift measured on the blinded test set does (+0.437 / +0.147 / +0.014 / -0.190), and
inversely as the effect of the joint likelihood at lambda = 3 does (-0.027 / -0.025 /
+0.047 / +0.120). Spearman +1.000 and -1.000. A reading that fits: the 750 compounds were
selected by the screen, so the enzymes whose screen is informative are the ones whose
actives got enriched — which would mean the coupling is partly learning the selection rule,
the middle factor of the factorisation in the model section.

This is four points. The permutation p-value for a perfect rank correlation at n = 4 is
1/24, and the screening correlation covaries with everything else that differs between
enzymes. It is one bit of evidence and the script says so. The within-enzyme version —
degrade the screening channel with noise in steps and watch whether the effect moves
monotonically — turns it into a dose-response curve on each enzyme separately, and has not
been run.

**29. The blind fraction was computed on the wrong derivative.** The joint-likelihood
section quoted 0.32 / 0.09 / 0.21 / 0.44 as the share of compounds where the fixed
calibration is uninformative about pi. Those are computed on |dI/dpi|, the derivative of the
inhibited fraction I = E/(1 + 10^(h(pC0 - pi))). What the loss actually compares against the
screen is g(pi) = log2(1 - I), so the informativeness that matters is |dg/dpi|, and the
fractions are 0.218 / 0.068 / 0.137 / 0.397. The ordering is identical under either
derivative, which is why the substitution never surfaced: the argument survives, the numbers
did not. A second slip travelled with it - |dg/dpi| is NOT maximised at pC0 but 0.19 to 0.35
units above it (4.529 / 4.494 / 4.657 / 4.600); only |dI/dpi| peaks exactly at pC0.

**30. The discriminator between saturation and bias absorption collapses.** The two stories
for why the single-latent arm is damaged were supposed to be separable by the per-enzyme
ordering, because the spread of E ranks 2D6 (0.091) above 3A4 (0.067) while the blind
fraction ranks them the other way. But 0.091 and 0.067 are q95 - q5 of the measured
per-compound `{CYP}_EmaxVsPosCtrl_direct_inhibition` column in
`data/cyp-challenge-TRAIN_Emax.csv` - a different quantity, in a different file, on a
different sign convention (measured Emax sits near -1; the calibration E is positive and
bounded in [0.2, 1.2]), over a different population (3A4: n = 2335 against n = 1805 in the
fit). `verify/g1_calib.py` computes no spread of E at all; it cannot be the source. The
fitted E's own spread is 30 to 60 times smaller: bootstrap sd 0.0026 / 0.0053 / 0.0031 /
0.0033 at B = 500, ordering 2C9 > 3A4 > 2D6 > 1A2, which is the opposite of the claim on the
decisive pair.

The quantity that does measure misfit of a single fixed (E, h) is the Hill residual sd that
g1_calib already prints: 0.231 / 0.185 / 0.497 / 0.612. It correlates with the measured-Emax
spread at Spearman -1.000 - the two candidate "misfit" proxies order the enzymes exactly
oppositely - and it puts 3A4 first, the same place saturation puts it. Under the right
proxy both stories predict the same ordering, so the per-enzyme breakdown does not separate
them and the noise-injection dose curve of item 28 is the replacement, not an extra.

Corroborating: inverting g per compound to get an implied E gives CYP3A4 an interquartile
range of 1.865 against 0.093 to 0.211 elsewhere, with 56.5 % of its compounds implying an E
outside the fit's own [0.2, 1.2] bounds against 9.6 to 12.9 % elsewhere. On that reading 3A4
is the worst-described enzyme, not the best.

**31. Out-of-fold isotonic does not preserve rank the way the argument needs - but the
distortion is common-mode.** Within one fold's map isotonic is exactly monotone: zero strict
inversions across all 20 enzyme x fold applications. "Ties aside" is not a small aside,
though - each fold's 242 to 496 compounds map onto 24 to 59 distinct values, and 95.9 to
97.0 % of compounds land in a collapsed tie. On the glued out-of-fold vector five different
maps are in play and rank genuinely moves: Kendall tau_b between raw and recalibrated is
0.917 to 0.959, 1.86 to 3.71 % of all pairs strictly reverse, and the worst compound shifts
143 to 223 positions even after breaking every isotonic tie in its favour. Spearman against
the LABELS - the thing the instrument is supposed to leave alone - falls on all four
enzymes, by 0.007 to 0.015. ST-RAE under isotonic also changes sign by enzyme: -0.0156 /
+0.0075 / -0.0504 / +0.0100, macro -0.0122.

What saves the two-arm comparison is that the distortion is nearly common-mode: +0.0132 /
+0.0135 / +0.0149 / +0.0070 on FP+DESC+MECH against +0.0139 / +0.0163 / +0.0141 / +0.0064
on FP+DESC, so the arm DIFFERENCE moves an order of magnitude less than either arm does.
Comparing two arms after out-of-fold isotonic stands; the sentence "it preserves rank" does
not.

**32. The affine pair's advantage is a delta = 0 advantage and reverses under tilt.**
Fitting (off, lambda) jointly out of fold beats the fixed +0.40 slice at delta = 0, 0.7150
against 0.7227. Under the label tilt the margin evaporates and turns over: 0.7910 against
0.7908 at delta = 0.3, and 0.9119 against 0.9061 at delta = 0.6. This is what a purely
SCALE advantage does when the mean moves, and the joint fit's is one. The family is not at
fault: refitting the pair UNDER the tilt gives 0.7579 at delta = 0.3 and 0.7892 at delta =
0.6, better than either delta-agnostic variant at every delta. The live question is
therefore not "affine pair or fixed offset" but "at which delta to fit it", and nothing
before the intermediate leaderboard answers that. `src/submit.py` fits at delta = 0, which
is not the best choice but is the only one that does not require guessing delta.

Related, on reading the reweight table: 0.7227 and 0.7831 do not come from the same column.
0.7227 is the fixed +0.40 variant at delta = 0, and it coincides with the oracle-offset
column there only because the oracle happens to pick +0.4 at delta = 0. 0.7831 is the oracle
column at delta = 0.3; the fixed +0.40 variant reads 0.7908 there. Quoting the pair as one
variant's movement understates it - the deployable movement, post-processing fitted at
delta = 0 because delta is unknown, is +0.0681 for the fixed slice and +0.0760 for the joint
fit against the oracle's +0.0604.

**33. `src/trunk.py` was silently overwriting one arm with the other.** Predictions were
saved under the key `{seed}|{lambda}` with no mode in it, and `--out` defaulted to the same
`results/preds/trunk.json` for both `--mode twohead` and `--mode calibrated`, so the second
run clobbered the first and nothing about the comparison survived on disk. Fixed: the key is
now `{mode}|{seed}|{lambda}`, the default path is `trunk_{mode}.json`, and a meta block
records the configuration and library versions the way `oof.meta.json` does. Every number
from the two arms is being recomputed; `src/trunkdose.py` consumes the result.

**34. Both mechanisms proposed for the single-latent arm's damage are refuted, and the
control is exact.** Both arms recomputed, four seeds, lambda in {0, 0.3, 1, 3}
(`src/trunkdose.py`). At lambda = 0 the two arms agree to a maximum absolute prediction
difference of exactly 0.0 on all four seeds - identity by construction, since the mode
branch sits inside `if lam > 0` and the weight seed is `hash((seed, fold))` with neither
lambda nor mode in it.

Per-enzyme loss at lambda = 3 against lambda = 0, mean over four seeds, calibrated arm:
1A2 +0.235, 2C9 +0.007, 2D6 +0.629, 3A4 +0.006. Saturation predicts the damage should
follow the blind fraction (0.22 / 0.07 / 0.14 / 0.40) and so be worst on 3A4; 3A4 is the
LEAST damaged and Spearman is -0.400, the wrong sign. The lambda slope agrees
independently: saturation wants the steepest response on 3A4 and the flattest on 2C9, and
the measurement is +0.178 on 2D6 against +0.006 on 3A4 and +0.005 on 2C9. Bias absorption
fares no better - Hill residual sd (0.231 / 0.185 / 0.497 / 0.612) gives Spearman -0.200,
also the wrong sign.

What does fit, at Spearman -1.000, is the screen's rank correlation with pIC50 (0.936 /
0.896 / 0.862 / 0.828) - the same factor that orders the two-head arm's effect. The fixed
calibration appears not to have a failure mode of its own: it multiplies the existing one
about fivefold on 1A2 and 2D6, and cancels the gain the two-head arm makes on 2C9 and 3A4
(-0.025 and -0.027 there, against +0.007 and +0.006 for the calibrated arm). n = 4 still,
so this is one bit in favour - but the two refutations are refutations, and a reversed sign
is evidence against rather than weak evidence for.

**35. The 0.006 gap is a small stable effect, not noise.** Out-of-fold isotonic on both
arms at lambda = 3: calibrated minus two-head is +0.0026 / +0.0042 / +0.0069 / +0.0097 by
seed, mean +0.0059, sign holding 4 of 4. Against the raw across-seed sd (0.021 at lambda =
0) it looks like a quarter of the noise; against the sd that survives recalibration
(two-head after isotonic, sd 0.0029, range 0.0068) it is 2.02 sd and 0.86 of the range.
Both arms are recalibrated by the time the gap is taken, so the second ruler is the right
one. The wording that follows is "the form of the coupling costs 0.006 of rank
information - small, but measurable", not "within noise".

**36. Trunk against boosting under the tilt: no reversal, but the gap narrows.** With the
same post-processing on both sides (the affine pair), boosting leads by 0.089 at delta = 0
and by 0.035 at delta = 0.6 - the direction expected if boosting's advantage is one of
scale and the trunk's is one of rank, but not enough to turn over inside the plausible
delta range. Note the trunk prefers a different repair: isotonic gives it 0.7506 at delta =
0 where the affine pair gives 0.8037, while for boosting the affine pair is much the better
of the two. Computed on split seed 0 only, and seed 0 is where the calibrated arm behaves
worst; this needs four seeds before it carries weight in a submission decision.

**37. The noise ladder: what it bought, and what it took away.** The screening channel was
degraded in five steps (eta = 0 / 0.5 / 1 / 2 / 4, rank attenuation 1.000 / 0.894 / 0.707 /
0.447 / 0.243), both arms, four seeds, lambda = 3. Noise is added to the standardised target
and the sum divided by sqrt(1 + eta^2), so only the correlation with pIC50 falls and the
target's spread does not - otherwise "less information" would be confounded with "more
weight". Control: at lambda = 0 the screening term is absent from the loss, and predictions at
eta = 2 reproduce those at eta = 0 to a maximum absolute difference of exactly 0.0 on all four
seeds.

What it bought: within-enzyme causality. Inside one enzyme nothing varies but how much the
channel knows, so a reproducible response there cannot be attributed to label count, band
width or Hill fit quality. Seven of eight curves hold their sign across all four seeds.

What it took away, and this is the larger result: **the single-latent arm's damage is not
informational.** Destroy 95-97 % of the channel's information and the damage does not vanish -
it grows, by +0.15 / +0.25 / +0.20 on 1A2, 2C9 and 3A4. The surviving few percent cannot do
more work than all of it did. The document's previous explanation - "the fixed calibration has
no disease of its own, it multiplies the existing one about fivefold" - is refuted by
measurement, and the sign is the refutation. Separately, and negatively: the two-head arm shows
no per-enzyme structure at all (mean deviations from a single curve at most 0.009 against a
spread of 0.033, indistinguishable from a permutation null). All of the per-enzyme structure
sits in the rigid coupling, not in the channel.

The reading that survives is single-mechanism: the channel pulls predicted potency toward a
pseudo-label derived from the screen; noise sets the label's QUALITY, which gives the
within-enzyme response, while the instrument sets the SCALE of the pull in pIC50 units, which
is fixed for the rigid arm and absent for the free head.

**38. The four-enzyme ordering is not identified, and one competitor has been eliminated.**
Three quantities order the four enzymes identically, pairwise Spearman exactly +-1: the screen's
rank correlation with pIC50; the mean confidence-band width, which is what the metric computes
error from; and the gain of the calibration curve, i.e. how many pIC50 units of spurious pull
one standard deviation of screening residual produces. No test in the joint-likelihood section
separates them, and noise cannot: it moves informativeness and moves neither the band nor the
instrument map.

Band width is now eliminated, cheaply. Recomputing ST-RAE from the same saved predictions with
the band equalised across enzymes, and again with no band at all (plain RAE), leaves the
per-enzyme structure essentially intact - 2D6 worst by a wide margin, 1A2 second, 2C9 and 3A4
near zero. Spearman moves from -1.000 to -0.800, and the pair that swaps is exactly 2C9/3A4,
which separate by 0.002 and were never distinguishable. So the effect lives in the predictions,
not in how the metric reads them.

The discriminating intervention is named and not run: swap (E, h) between CYP2D6 and CYP3A4 and
rerun the single-latent arm on the same four seeds. Data, labels, bands, split, initialisation
and the screen's informativeness all stay literally the same; only the instrument map moves.

**39. Three claims in section 14 were overstated and are now weakened.** The permutation
probability of a perfect ordering at n = 4 was quoted as "about one in twenty"; it is 1/12
two-sided. The two-head arm's perfect ordering at eta = 0 rests on a single run - per seed the
CYP2D6 effect is +0.609 / -0.040 / -0.121 / +0.034, and dropping seed 0 turns Spearman -1.000
into +0.400 - while the 2C9/3A4 pair differs by 0.002 against a standard error four times
larger; exactly one comparison in that arm is measured, 1A2 against the pair. And "two checks
out of two, the hypothesis fails" claimed too much from correlations of -0.400 and -0.200 at
n = 4, which sit in the middle of the permutation null: the honest statement is that the
predicted ordering is not observed. The single-latent arm's ordering, by contrast, is stable -
it holds under the median and without seed 0.

**40. Swapping the instrument calibration between two enzymes: the map does not order them.**
Three quantities order the four enzymes identically (item 38); the band was eliminated by
rescoring, leaving the screen's informativeness and the calibration gain. Noise cannot
separate those two, because it moves informativeness and leaves the map untouched. Swapping
(E, h) between CYP2D6 and CYP3A4 is the one intervention that does the reverse: the compounds,
labels, bands, split, initial weights, batch order and screening readout all stay literally
identical, and only the instrument map moves. `src/trunk.py --swap-cal CYP2D6,CYP3A4`, single
latent arm, lambda 0 and 3, four seeds, about fifteen minutes.

Two controls, both passed. At lambda = 0 the screening term is absent from the loss, so
g_of_pi is never called and the swap cannot reach the predictions: bit-identical on all four
seeds. And CYP1A2 and CYP2C9 keep their own calibrations throughout, so their movement
measures what a rerun costs: 0.019 and 0.009.

The prediction was that CYP2D6, handed CYP3A4's steeper map, should be hurt far less than its
usual +0.63, and CYP3A4, handed CYP2D6's shallower one, should go from near zero to tenths.
Neither happened. Excluding seed 0, CYP2D6 moved +0.003 and CYP3A4 moved +0.016 - both smaller
than the untouched CYP1A2 control. The damage stayed with the enzyme. The calibration gain is
therefore out, eliminated by intervention rather than by correlation, and what remains is the
screen's informativeness or something else in that enzyme's own data that has not been named.
The experiment establishes where the mechanism is not.

**41. Split seed 0 has now produced three separate false conclusions in the single-latent
arm.** *(Superseded by item 42: the cause is one compound, not the split. Kept because the
detection was right even though the diagnosis was wrong.)* In the swap above its CYP2D6 movement is +0.61 against +0.00 on the other three, and
the mean over four seeds reads +0.155 - fifteen times the control noise and apparently
decisive - while the median and the three-seed mean read zero. The same seed alone produced
the perfect four-point ordering of item 39, where dropping it turns Spearman -1.000 into
+0.400. And it is the seed on which the arm diverged to NaN at two noise rungs (item 37).
Three conclusions in one section have leaned on one split. The four-seed rule in README's "How
we work" exists for this, and it has now paid for itself three times; any per-enzyme claim in
the joint-likelihood work should be read per seed before it is believed. Why that split is
different has not been investigated.

**42. It was never the seed. It was one compound, and it overturns one of the two headline
orderings.** Item 41 blamed split seed 0 for three separate false conclusions. The detection
was right and the diagnosis was wrong. Seed 0's split is unremarkable: fold sizes, mean
CYP2D6 activity, spread, active fraction and the metric's own denominator all match the other
three seeds to two decimals.

What is remarkable is a single molecule. In the two-head arm at lambda = 3, `OCNT-2328942`
(true CYP2D6 pIC50 2.53) is predicted at -360 on seed 0 and -22 on seed 2, and behaves
normally on the other two. Its screening reading is ordinary (z = -0.58), so this is an
optimisation blow-up that lands on it, not an outlier in the data. ST-RAE is a sum of absolute
deviations, so one prediction at -360 contributes about 362 to a numerator whose denominator
is around 120: one molecule out of 1493 triples the enzyme's score.

Clipping predictions to the enzyme's label range plus or minus two units - a bound that
touches no compounds at all on three enzymes and half a compound on CYP2D6 - settles it:

  two-head 2D6, raw      +0.609 / -0.040 / -0.121 / +0.034   sign 2 of 4
  two-head 2D6, clipped  +0.035 / +0.027 / +0.027 / +0.034   sign 4 of 4
  single-latent, raw     +0.941 / +0.536 / +0.425 / +0.613   sign 4 of 4
  single-latent, clipped +0.575 / +0.574 / +0.587 / +0.596   sign 4 of 4

The two arms come apart. **The two-head arm's ordering does not survive**: its CYP2D6 damage
is +0.031 rather than +0.120, CYP2D6 stops being the worst enzyme (CYP1A2 takes it at +0.047),
and the Spearman against screening informativeness falls from -1.000 to -0.800. The claim that
the channel hurts most on the enzyme the whole document is built around rested, in that arm,
on one molecule appearing on two seeds of four.

**The single-latent arm's ordering survives and gets stronger.** CYP2D6 stays worst by a
factor of two and a half over CYP1A2, the ordering is unchanged, the Spearman stays -1.000 -
and clipping *stabilises* it, dropping the across-seed spread from 0.516 to 0.022, a factor of
twenty-three. Everything the section says about the rigid coupling now rests on numbers that
do not depend on which seed you take.

Not a submission risk today: the submitted model is the boosting, which does not extrapolate
past its label range, and all 750 predictions in `results/submission/` sit inside it.
`src/submit.py` has no clip, though, so nothing in the pipeline would stop a -360 if the trunk
were ever the model shipped.

**43. Trunk against boosting, on four seeds instead of one: the gap is a tenth of what the
post-processing is worth.** Block 6 of `src/trunkdose.py` is the only place the neural model
and the boosting are compared under a tilted label marginal, and it is the only number in the
joint-likelihood work that bears on what gets submitted. It ran on split seed 0 alone - the
seed the runaway compound of item 42 lands on. Recomputed on four seeds, with trunk
predictions clipped to the label range plus or minus two units, and the same post-processing
(the affine pair) on both sides. The boosting side exists for all four seeds: seed 0 in
`oof.json`, seeds 1-3 in `oof_seeds.json` from `verify/f3_seeds.py`.

The previous conclusion was "boosting leads, and the gap narrows by about half under the
tilt". Neither half survives. The gap is not 0.089 but 0.000 to 0.005, two orders of magnitude
smaller, and it does not narrow - it grows: the models are indistinguishable at delta = 0 and
the boosting is 0.005 ahead by delta = 0.6. The direction is the opposite of what was expected
from "boosting's edge is scale, the trunk's is rank".

Not all of it is measured. Below delta = 0.3 the sign of the difference does not hold across
four seeds, so the two models are simply indistinguishable there. From 0.3 to 0.6 the sign
holds and the boosting is ahead, by 0.0034 to 0.0051.

The number worth carrying: the whole difference between the two model families never exceeds
0.0051, while the post-processing is worth 0.0512 - **ten times more**. In this range the
choice of architecture decides almost nothing and the choice of what to do with the
predictions afterwards decides almost everything. That is where the remaining effort belongs.

`fit_affine_oof` was vectorised to make four seeds cheap: the ST-RAE denominator does not
depend on (off, lambda), so minimising the metric is minimising its numerator, and the
numerator over the whole grid is one broadcast. About four times faster, and verified to pick
the same (c, L) and produce identical predictions on every fold.

**44. Where to fit the affine pair, given that delta is only bracketed - and it matters twenty
times more than the model choice.** Item 43 left the post-processing worth 0.051 against 0.005
for the whole boosting-versus-trunk question. That post-processing has a free parameter we
cannot observe: the pair is fitted on our own label marginal, at delta = 0, while the test's
marginal sits somewhere in +0.1 to +0.6.

The object that answers it is a matrix, not a number: fit under an ASSUMED delta, score under
a TRUE delta, out of fold throughout, four seeds (`src/shrinkchoice.py`). Its diagonal is the
oracle; its top row is what `src/submit.py` does today. Over the bracket:

  fit at 0 (today)        worst 0.913   mean 0.816   at 0.1  0.734   at 0.6  0.913
  fit at the mid-bracket  worst 0.826   mean 0.776   at 0.1  0.747   at 0.6  0.826
  minimax over bracket    worst 0.800   mean 0.786   at 0.1  0.800   at 0.6  0.796
  oracle (unreachable)    worst 0.792   mean 0.764   at 0.1  0.731   at 0.6  0.792

Fitting at zero is the worst of the three by both criteria - 0.087 worse in the worst case,
0.040 worse on average - and its sign holds on all four seeds for every test delta at 0.2 and
above. It wins in exactly one place, the bottom edge of the bracket, by 0.013. So "fit at
zero" is not the choice that avoids an assumption; it is the choice that assumes the test
looks like the training set, which the same document spends a section showing it does not.

Per-enzyme the shifts differ, and on CYP2D6 the estimate is *negative*. Fitting each enzyme at
its own estimate would be worth up to 0.022 more than one global delta - but that number is
the per-enzyme oracle conditional on `verify/k5_shift.py` being exactly right, and it is an
upper bound, not a gain. Note also that the k5 estimate for 2D6 contradicts the +0.1 to +0.6
bracket outright: one of the two measurements is wrong there.

`src/submit.py` gains `--delta`, defaulting to 0 so nothing changes by itself. Verified the
default path is untouched: at delta = 0 the weights are all ones, the ST-RAE denominator does
not depend on the parameters, so the weighted-numerator objective picks the identical (off,
lambda) on all four enzymes - checked, not argued. Which delta to actually ship is a decision
for the team, not a side effect of an edit.

**45. The CYP2D6 contradiction was never a contradiction: the two numbers measure different
things, and the negative shift is real chemistry.** Item 44 left the anchor argument giving
+0.19 for CYP2D6 and `verify/k5_shift.py` giving -0.19, with one of them presumably wrong.
Neither is. `src/reweight.py` hard-codes anchor percentiles 93 / 98 / 62 / 90, and its own
prose says why 2D6 is 62: **that enzyme took no part in anchor selection**. Its number is the
internal control - an enzyme not used to pick anchors should show no enrichment, and it shows
none. It was never an estimate of 2D6's shift, and the +0.1 to +0.6 bracket comes from the
three enzymes that were used. Carrying that bracket over to 2D6 was our error, not a
disagreement between measurements.

The measurement then stands on its own, and it is tight: the prediction shift is -0.190 with a
bootstrap interval of [-0.234, -0.145] over the 750 test compounds. A model can move its own
output without the labels moving, though, so it needs a mechanism, and there is one that can
be checked without a single test label.

CYP2D6 binds through a salt bridge from an active-site aspartate and glutamate to a
protonated basic nitrogen; the other three bind by lipophilicity, planarity or an anion. The
test set is depleted in exactly that chemistry, by about half: the 2D6 pharmacophore feature
falls 0.148 to 0.083, the fraction basic at pH 7.4 falls 0.175 to 0.104, tertiary amines and
piperidines roughly halve, cation count falls 0.610 to 0.400. Every one of ten such features
moves down, all at |z| > 4.8.

And on our own labels, **CYP2D6 is the only enzyme where basic compounds are more active** -
by +0.55, against -0.24, -0.44 and -0.45 on the other three. So one composition shift lowers
2D6 and slightly raises the rest, which is what the measurement shows: the sign agrees on all
four enzymes. The magnitude does not - a single binary feature accounts for about a fifth of
2D6's shift and less elsewhere, where most of it comes from the activity enrichment the set
was built for. It is an estimate of direction, not a model of the shift, and `k7` says so.

The consequence for the submission is concrete: **the marginal shift is not one number for all
four enzymes.** Fitting the post-processing to a single global delta fits CYP2D6 in the wrong
direction outright. `verify/k7_2d6shift.py`.

**46. Per-enzyme fitting, and the global rule actively harms CYP2D6.** Item 45 established
that the marginal shift is not one number for four enzymes. The consequence is measurable
(`src/shrinkchoice.py`, block 5). Each enzyme gets its own plausible range for the true shift,
built from what was measured: the direct estimate is one edge, since it is a lower bound in
magnitude, and a margin in the direction of its sign is the other; for CYP2D6 the second edge
is zero, because the anchor control says there is no enrichment there. Rules are compared by
their WORST case inside each range rather than at a point - choosing a rule to match an
estimate and then scoring it at that same estimate would be one action, not two.

  worst case in own range   1A2     2C9     2D6     3A4    macro
  fit at zero (today)     0.911   0.878   0.913   0.762   0.866
  one global delta = 0.3  0.878   0.790   1.035   0.710   0.853
  per-enzyme              0.878   0.777   0.906   0.666   0.807
  chosen shift             +0.3    +0.4    -0.1    +0.8

Per-enzyme fitting is worth 0.046 over the best global rule and 0.059 over current behaviour.
The macro is not the interesting part. **The best global rule makes CYP2D6 worse than doing
nothing at all** - 1.035 against 0.913. A single positive delta does not merely underperform
on 2D6; it actively damages it, because the real shift there has the opposite sign.

The minima are interior, not artefacts of where the grid stops: rechecked on a grid extended
to +1.6, CYP3A4's worst-case curve descends to +0.8 and rises after it (0.663, 0.666, 0.674,
0.690), and CYP1A2's turns at +0.3 with both neighbours higher.

`src/submit.py --delta` now takes either one number or four comma-separated, defaulting to
zero so nothing changes by itself. Changing the default changes what gets submitted and is the
team's decision.

**47. The attempt to narrow the delta ranges mostly failed, and found two of our own errors on
the way.** The plan was to measure the propagation factor - what fraction of a label shift
reaches the model's predictions - and so turn the direct measurement from a bound into a point
estimate. Three subsampling designs gave 0.23-0.63 (random Butina clusters), 0.9-1.05 (sliding
windows along continuous descriptors) and 0.935 (fixed basic-amine fraction on CYP2D6). The
conclusion drawn was that propagation is essentially complete and the bracket collapses to the
bootstrap. That conclusion does not survive.

**The cluster design is `b` rediscovered, not a second lever.** For a random subsample the
forward slope is Cov/Var(y) and the reverse slope is Cov/Var(yhat), which is exactly the
attenuation b, and their product is R^2. Check: R^2/slope = 0.795 / 0.918 / 0.778 / 0.972
against b = 0.789 / 0.898 / 0.743 / 0.999. The two "independent" handles are algebraically one.

**Dividing by b was wrong, but not for the reason the document gave.** b is the share of the
OUTPUT deviation that is real - regression dilution, which section 11 already names. The share
of an INPUT shift reaching the output is the other slope, R^2/b = 0.29 / 0.41 / 0.22 / 0.59.
Dividing by b produces neither, so +0.018 / +0.164 / -0.256 / +0.437 are withdrawn, and the
claim they were a lower bound on |delta| was unfounded.

**Propagation is not a scalar.** Ten of twenty-eight window cells fall outside [0.9, 1.05], and
values of 1.22 and 1.55 are impossible for a "fraction that reaches the output" - the window is
a Wald ratio with the descriptor as instrument, and the model sees that descriptor directly, so
the exclusion restriction fails by construction. Disagreement across seven instruments is the
standard sign they are invalid. The tidy "0.9 to 1.05" summary was obtained by dropping the two
axes on which item 45 built the actual mechanism.

**The intervals should be wider, not narrower.** The test is anchors plus analogues, not 750
independent compounds: clustered at 0.50 it is 172 groups. A cluster bootstrap gives a design
effect on the variance of 4.4 to 6.7, so intervals are two to three times wider than the naive
ones: 1A2 [-0.111, +0.124], 2C9 [+0.025, +0.257], 2D6 [-0.281, -0.093], 3A4 [+0.295, +0.576].

Net effect on the ranges: they move rather than shrink. CYP2C9 is the only genuine narrowing.
CYP1A2 widens and now contains zero. CYP3A4 moves down about 0.1. **CYP2D6 is the one real
gain: zero leaves its range**, and it is the only enzyme whose interval excludes zero under
either bootstrap. Refitting the per-enzyme rule on the honest ranges: per-enzyme is worth 0.037
over the best global rule and 0.028 over today's behaviour, and a single global delta = 0.3 is
now worse than doing nothing even on macro (0.817 against 0.808), because the harm it does to
CYP2D6 outweighs the gain on the other three.

**48. A mask bug in `verify/k7_2d6shift.py` understated its own finding fourfold.** The basic
fraction was taken over all 4905 rows (0.175) while the activity contrast was measured inside
each enzyme's label mask. The masks differ enormously in composition, because dose-response
curves were run on the basis of the screen: the CYP2D6 mask is 0.355 basic against 0.114-0.122
for the other three. The test is fully labelled at 0.104. So the correct difference for 2D6 is
-0.251, not -0.071, the expected shift is -0.138 against -0.190 measured, and the salt-bridge
mechanism explains about **three quarters** of the observed shift rather than "about a fifth".
The reported figure for the other three enzymes stays at a few percent. Fixed.


**49. The propagation factor does not have to be estimated at all.** `k8_kernel.py`. Item 47
closed three designs for it and was right about all three, but every one of them was an attempt
to estimate a multiplier. The multiplier can be bypassed. The identity

    E_test[yhat] = INT E[yhat|y] * p_test(y) dy

holds exactly whenever the kernel E[yhat|y] is the same on both sets. It needs the *reverse*
regression — whose slope in the linear case is the R^2/b that item 47 identified as the correct
one — and here that kernel is estimated by isotonic regression, so linearity is not assumed
either. The training label marginal is then tilted until the implied prediction mean matches the
observed test mean. No instrument, so no exclusion restriction to violate; the reverse
regression, so no algebraic collapse onto `b`; and propagation never has to be a scalar, because
the direction is fixed by the tilt itself. The window cells at 1.22 and 1.55 are not evidence
against this — along a descriptor the model reads directly, propagation really is near one; that
is simply a different direction from the one the question is about.

| enzyme | d(yhat) | b | R^2 | R^2/b | *b/R^2 linear | isotonic kernel | 95 % clustered |
|---|---|---|---|---|---|---|---|
| CYP1A2 | +0.014 | 0.789 | 0.215 | 0.273 | +0.050 | **+0.045** | −0.059 … +0.156 |
| CYP2C9 | +0.147 | 0.898 | 0.364 | 0.405 | +0.364 | **+0.362** | +0.288 … +0.434 |
| CYP2D6 | −0.190 | 0.743 | 0.146 | 0.196 | −0.968 | **−0.917** | −1.151 … −0.715 |
| CYP3A4 | +0.437 | 0.999 | 0.592 | 0.593 | +0.737 | **+0.740** | +0.663 … +0.810 |

A null worth having: isotonic and linear agree to 0.05, so kernel nonlinearity contributes
nothing and the whole correction is the b -> R^2/b substitution.

Against the shifts item 46 chose by worst-case reasoning — +0.3 / +0.4 / −0.1 / +0.8 — this
agrees closely on the two enzymes that carry the decision (CYP3A4 +0.740 against +0.8, CYP2C9
+0.362 against +0.4, both by a completely different route) and disagrees sharply on CYP2D6
(−0.917 against −0.1) and mildly on CYP1A2 (+0.045 against +0.3, and its interval contains zero
exactly as item 47 says).

The disagreement is not a tie to be split, and CYP2D6 is where this estimate should be trusted
least rather than most. Its one assumption is kernel invariance, and item 45 is precisely a
finding that the test's composition changed in a chemically specific way — the basic-amine
depletion — on the enzyme where bases are the active class. If the test's low-activity compounds
are low for a different structural reason than the training set's, E[yhat|y] is not the same
function and the inversion is biased there. The fix follows from item 45 rather than contradicting
it: estimate the kernel separately inside and outside the basic stratum and mix the two at the
test's measured composition, which replaces the invariance assumption with a measured mixing
weight on the one axis known to have moved. Not yet run.

**50. Under an asymmetric cost the action is a quantile, not a midpoint.** With underestimation
costing 0.087 and overestimation 0.013 on CYP3A4, a piecewise-linear loss is minimised at the
quantile of level 0.087/(0.087+0.013) = **0.87** of the posterior for delta, not at its centre.
On item 47's range that is +0.539 rather than +0.435; on this file's bootstrap, +0.782 rather
than +0.741. This is the same Bayes-point argument the document already makes in section 10 for
the point prediction under ST-RAE, applied one level up, and it means the range does not have to
be narrowed before it can be acted on. The caveat is the shape of the loss: 0.87 is exact if
0.087 and 0.013 are slopes at comparable distances, and if instead they are costs at the ends of
the range with curvature in between, the expectation should be minimised over the measured curve
— the answer moves but stays well above the midpoint.

**51. One delta is not enough: the test is shifted *and* widened, and that partly pays the shift
back.** `k9_shape.py`. Exponential tilting matches the test's prediction mean by construction and
misses its spread: 0.625 / 0.504 / 0.517 / 0.718 against the test's 0.708 / 0.740 / 0.526 / 0.928,
with KS rejecting on three enzymes of four. Tilting can move a mean; it cannot add variance. Every
scheme in this repository that is parameterised by a single delta — `reweight.py` and
`shrinkchoice.py` included — is therefore incompletely specified.

This reverses one premise. "A sample narrower in activity has a smaller denominator, so ST-RAE
rises anyway" assumes narrowing; the test is **wider**, by 1.13 / 1.40 / 1.05 / 1.10 in the
standard deviation of the predictions, and at least that much in label space provided the kernel
noise is no larger on the test — which the similarity geometry supports, since the test sits
closer to the training set (0.587) than the training set does to itself (0.435).

Fitting both moments degenerates on CYP2C9 exactly as the weights in `k5_shift.py` did —
effective sample size 20 of 1285 — and those numbers are discarded rather than reported. The
conservative version widens the label distribution only as far as the predictions widened, a
lower bound, and keeps an effective size of 492 to 1355. Scored on the full ST-RAE, numerator
included:

| enzyme | delta = 0 | shift only | shift + widening |
|---|---|---|---|
| CYP1A2 | 0.8786 | 0.9028 | 0.8491 |
| CYP2C9 | 0.6908 | 0.8716 | 0.7692 |
| CYP2D6 | 0.9803 | 0.8104 | 1.0334 |
| CYP3A4 | 0.5194 | 0.8029 | 0.6623 |
| **macro** | **0.7673** | **0.8469** | **0.8285** |

Widening returns about a quarter of what the shift costs, not all of it. The expectation for the
intermediate leaderboard is therefore around **0.83 macro rather than 0.767**, and that is a
single number, checkable on 24 September, which is the cheapest test any of this has.

**52. Our delta ranges were on the wrong scale, and the fix nearly doubles the case for
per-enzyme fitting.** Item 46 built the ranges from a cluster bootstrap of the shift in
PREDICTIONS and used them as ranges for the shift in LABELS. Those differ by b/R^2 = 3.40 /
2.43 / 4.48 / 1.69, so all four were wrong and the chosen shifts derived from them were void.
Found by external review, though by a wrong route - the reviewer inferred we were still
dividing by b, when in fact we were not converting the scale at all.

Rebuilt on the kernel estimate of item 49, the picture is sharper than before:

  worst case in own range     1A2     2C9     2D6     3A4    macro
  range                    -0.1..0.2  0.3..0.4  -1.2..-0.7  0.7..0.8
  fit at zero (today)       0.874   0.829   1.007   0.762   0.868
  one global delta = 0.3    0.864   0.761   1.196   0.710   0.883
  per-enzyme                0.860   0.758   0.890   0.666   0.793
  chosen shift               +0.2    +0.4    -1.1    +0.8

Per-enzyme is now worth 0.090 over the best global rule and 0.075 over current behaviour,
against 0.037 and 0.028 on the wrong scale. And a single global delta = 0.3 is now clearly
worse than doing nothing (0.883 against 0.868), because CYP2D6 under a positive global shift
goes to 1.196.

**53. Both external scripts reproduce in this environment.** `verify/k8_kernel.py` and
`verify/k9_shape.py` arrived committed but had only been run in a mirror of the layout, not
through this repo's pinned interpreter. Run here they reproduce their reported numbers exactly:
kernel deltas +0.045 / +0.362 / -0.917 / +0.740, the linear limit d_yhat*b/R^2 agreeing
everywhere except CYP2D6 where R^2 is lowest, and the shape result - test predictions are
1.13 / 1.40 / 1.05 / 1.10 times WIDER than out-of-fold, KS rejecting on three enzymes of four.
A single delta is therefore an incomplete specification of the shift for every script that
uses one, `src/reweight.py` and `src/shrinkchoice.py` included. Expected intermediate
leaderboard macro is about 0.83 rather than 0.77, and the widening returns about a quarter of
what the shift costs rather than cancelling it.

**54. The stratified kernel on CYP2D6: nearly half of what we called a label shift was
composition.** Item 49's kernel needs one assumption, that E[yhat|y] is the same function on
the test as on the training set, and item 45 showed that assumption is weakest exactly where
the answer matters: the test carries 0.104 basic compounds against 0.355 in the CYP2D6 label
mask, and CYP2D6 is the enzyme that binds them.

Splitting the kernel by that stratum and mixing the two by the TEST's measured composition
rather than ours turns an assumption into a measurement. The strata really do have different
kernels on CYP2D6, by 0.35 to 0.58 across the whole scale - the model scores basic compounds
high there almost regardless of their true activity - while on the other three enzymes the two
kernels nearly coincide.

  enzyme    pooled    stratified   of which composition   of which labels   change
  CYP1A2    +0.045      +0.038            +0.003              +0.035       -0.007
  CYP2C9    +0.362      +0.369            +0.005              +0.365       +0.007
  CYP2D6    -0.917      -0.508            -0.138              -0.371       +0.409
  CYP3A4    +0.740      +0.742            +0.008              +0.734       +0.003

The control is the point: the three enzymes whose masks differ from the test by about one
percentage point move by 0.003 to 0.007, which is nothing. CYP2D6, whose fraction differs
threefold, moves by 0.409. Cluster bootstrap on the stratified estimate [-0.693, -0.325], zero
still excluded.

Consequence for the choice: CYP2D6's fitted shift goes from -1.1 to -0.3, per-enzyme fitting
is worth 0.053 rather than 0.090, and a single global delta = 0.3 is now level with doing
nothing (0.8425 against 0.8428) rather than clearly worse. The assumption has been narrowed,
not removed - invariance is still assumed within each stratum, and a shift along some other
axis relevant to CYP2D6 would not be caught.

**55. The bootstrap in both kernel scripts was half a bootstrap, and fixing it widens every
interval by two to five times.** `verify/k8_kernel.py` resampled training clusters while
holding the observed test mean fixed, so it captured the uncertainty of the kernel and threw
away the uncertainty of the target. Under inversion the target enters with a factor of one
over the kernel slope, and those slopes are 0.29 / 0.41 / 0.22 / 0.59, so the discarded part
was larger than the part kept. Our `k10_strat2d6.py` inherited the same construction. Both now
resample both sides, the test in blocks by analogue series.

  enzyme    delta    old 95%              honest 95%
  CYP1A2   +0.038   [-0.064, +0.153]     [-0.331, +0.437]
  CYP2C9   +0.369   [+0.288, +0.445]     [+0.055, +0.723]
  CYP2D6   -0.508   [-0.693, -0.325]     [-1.140, +0.093]
  CYP3A4   +0.742   [+0.666, +0.799]     [+0.484, +1.030]

**And the strongest claim in the section survives only barely.** Under the unstratified kernel
CYP2D6's honest interval is [-1.505, -0.505] and excludes zero, which is what the external
review predicted. Under the stratified kernel - the better estimate - it is [-1.140, +0.093]
and does not. Measured directly on 1500 draws, P(delta >= 0) = 0.039: the 90 % interval
excludes zero, the 95 % interval does not. "The shift on CYP2D6 is negative" remains the
firmest of the four statements; "clearly different from zero" is no longer the right wording
for it.

**56. Wider ranges make per-enzyme fitting more valuable, not less.** The prediction was that
honest intervals would shrink the case for per-enzyme fitting to "CYP2D6 and CYP3A4, zero
elsewhere". The opposite happened, and the reason is structural: a wider range makes the worst
case worse for any single fixed rule, while per-enzyme fitting gains room to hedge each enzyme
separately. Per-enzyme is now worth 0.091 over the best global rule and 0.107 over current
behaviour, against 0.053 and 0.053 on the narrow ranges.

The qualitative prediction inverted too. CYP2D6's range now straddles zero, so the worst-case
rule picks almost nothing there (-0.1), while the substantial shifts go to CYP1A2 (+0.4) and
CYP2C9 (+0.5). The enzyme that motivated per-enzyme fitting is now the one it barely touches.

Two smaller corrections accepted from the same review. The composition component of -0.138 in
item 54 is not independent confirmation of item 48 - it is the same product, the change in
basic fraction times the activity contrast, computed a second way. And the residual
uncertainty on CYP2D6 is set by conditioning rather than by which correction was applied: the
kernel slope there is 0.224 against 0.593 on CYP3A4, and 0.119 within the basic stratum, so
inversion amplifies error by 8.4 times where CYP3A4 amplifies it by 1.7.

**57. More than half the per-enzyme gain was the criterion, not the estimates.** Every rule
comparison so far picked by the WORST case inside each enzyme's range, which is insurance
against a bad leaderboard rather than a bid for the best expected score. The criterion was
inherited, never chosen. With 1500 posterior draws of delta per enzyme (`k10` now saves them
to `results/preds/delta_draws.json`) the same search can be run under the MEAN instead:

  per-enzyme gain over fitting at zero:   worst case +0.107   mean +0.045

The two criteria differ in the sign of their derivative with respect to range width - a wide
range pulls the averaging rule toward a single global number and pushes the worst-case rule
away from one - so this is not a detail. It is also a question about the goal, and it should
be answered out loud before 3 November rather than after.

**58. Our own rule carried the defect we rejected the global rule for.** A single delta = 0.3
was rejected because on CYP2D6 it is worse than doing nothing. Held to the same standard, the
worst-case per-enzyme rule picks +0.4 on CYP1A2 and is worse than doing nothing in **81 %** of
the posterior draws, average loss 0.061 where it loses. The defect had simply moved to another
cell, and it was harder to see there because it was in our favour.

  fraction of draws where the rule loses to doing nothing
                          1A2     2C9     2D6     3A4
  P(delta >= 0)          0.59    0.99    0.05    1.00
  worst-case pick        +0.4    +0.5    -0.1    +1.0
    loses to nothing     0.81    0.28    0.07    0.13
  mean pick              +0.1    +0.3    -0.3    +0.7
    loses to nothing     0.51    0.09    0.11    0.01

The defensible rule that survives both checks: **fit per-enzyme on CYP2C9, CYP2D6 and CYP3A4,
and leave CYP1A2 alone.** On CYP1A2 the sign of the shift is not determined at all - P(delta
>= 0) = 0.59 - so any non-zero choice there is a coin flip against doing nothing. The mean
criterion picks +0.1, gains 0.0001 by it, and loses half the time: that is added variance for
no expected return. Zero there is not caution, it is the better move.

Note this is neither of the two qualitative predictions on offer. It is not "per-enzyme on
CYP2D6 and CYP3A4, zero elsewhere", and it is not per-enzyme everywhere.

**59. The offset grid was clipping the fit at both ends, and one adopted number was wrong
because of it.** `OFFGRID` ran from -0.2 to +1.6. At assumed shifts below -0.3 the fitted
offset hit the lower edge and stopped moving, which made the objective look flat from -0.3 to
-0.6 when it was not; at +0.7 on CYP3A4 it sat on the upper edge. A fitted parameter resting
on a grid boundary is not a fitted parameter, and flatness beside one is an artefact. Widened
to +-3 in `src/submit.py` and `src/shrinkchoice.py`, and `submit.py` now prints a warning if
an optimum lands on an edge. The search is vectorised at the same time - the whole grid in one
broadcast - because the wider grid made the Python loop too slow to finish.

The correction changes one adopted number: CYP2D6's mean-criterion pick moves from -0.3 to
**-0.5**, improving that cell from 0.8335 to 0.8241. The other three are unchanged. Adopted
vector: **0, +0.3, -0.5, +0.7**, and `src/submit.py --delta` now defaults to it.

**60. Break-even is the better question, and it sharpens the warning rather than settling it.**
The posterior integrates sampling noise and contains nothing about systematic error - which is
exactly what a failure of kernel invariance would be. So the question to ask of a rule is not
"how likely is it to lose" but "how far can the estimate be wrong before this cell loses to
doing nothing":

  enzyme   pick   E[delta]   break-even   margin
  CYP1A2   +0.1     +0.037      +0.05      0.01
  CYP2C9   +0.3     +0.361      +0.16      0.20
  CYP2D6   -0.5     -0.518      -0.24      0.28
  CYP3A4   +0.7     +0.733      +0.40      0.33

CYP1A2's margin is 0.01, which is independent confirmation that zero is right there - any
non-zero pick sits on a knife edge.

The uncomfortable one is CYP2D6. **Stratifying the kernel moved its estimate by 0.41, and the
margin remaining is 0.28.** A correction we have already applied is larger than the room left
before the cell turns. That is not an argument for hedging - a partial shift of -0.2 instead
of -0.5 costs 0.007 expected and buys only 0.03 of margin, a bad trade - it is an argument for
finding out whether a second composition axis of comparable weight exists. The first was found
by chemistry in an evening; whether it is the only one has not been asked.

**61. The representation was never the bottleneck, but the data was.** Every row of the
ablation grid is Morgan counts plus RDKit descriptors plus the mechanistic block, so the
conclusion that model choice is worth 0.005 against 0.051 for post-processing was drawn inside
one representation and written as though general. Two experiments that had never been run,
both cheap, both with the control recomputed rather than read from `oof.json`.

**Encoder: negative.** The organisers published sixteen chemprop pretraining checkpoints;
`rdkit2d` was chosen because it is pretrained on the very descriptors our DESC block contains,
which is the encoder's best case. Same learner, same folds, same metric. The control
reproduces exactly - 0.7673 against 0.7673, difference zero - so the rows are comparable.
Embedding alone 0.7784, embedding plus mechanistic block 0.7765, against 0.7673 for
FP+DESC+MECH. The pretrained representation loses by 0.009. `src/embed.py` runs the encoder in
a separate environment and leaves an array, so chemprop's forty-one dependencies never touch
this repository's measured pins.

**External data: positive, and the correction made it worse.** The same organisation shipped
the baseline's training set - 8068 compounds with pIC50 on our four enzymes, curated from
ChEMBL, Apache-2.0. Overlap with the blinded test is **zero of 750**; overlap with our
training set is 64 and those are dropped.

  change against control, negative is better    1A2      2C9      2D6      3A4    macro
  external as they are                       -0.026   +0.022   -0.036   -0.011   -0.013
  shifted by the paired offset                -0.019   +0.028   -0.041   +0.003   -0.007
  shifted by the marginal offset               -0.009   +0.082   -0.043   +0.116   +0.037

The two label sets differ in scale - external values run +0.22 / +0.35 / +0.60 / +0.56 higher
on shared compounds, and up to +1.36 by marginal - and correcting for it makes things worse,
the crudest correction worst of all. The apparent gap is mostly real enrichment rather than
miscalibration: ChEMBL publishes what worked, and shifting those labels down corrupts values
that were right.

The gain lands where this document predicted it would. **CYP2D6 improves most**, and CYP2D6 is
the enzyme with the worst R^2 of the four, the one the document argued is short of information
rather than short of calibration. No post-processing repairs an R^2 of 0.17; half again as
many labels does. CYP2C9 is the only enzyme harmed, under every correction, and it also has
the fewest shared compounds - six - so its offset is the least known of the four.

**One split seed.** This repository's rule is four, and until then 0.013 is an indication
rather than a result. It is nonetheless twice everything the choice between boosting and the
neural trunk is worth, from a source untouched for two months.

**62. The zero overlap with the test says nothing about the test.** `k12_extneighbors.py`.
Exact overlap is the wrong check on its own: a compound one methyl from a test compound leaks
and matches nothing. Measured for every test compound its highest Tanimoto to the external
set — median 0.329, 75th percentile 0.382, **fraction above 0.7 exactly zero, above 0.9 exactly
zero**. Our own compounds sit at 0.450 from each other, so the external set is further from the
test than our training set is from itself. Leakage is closed by proximity, not only by identity.

The zero then invites an inference, and the obvious one is wrong. Overlap with our training
set is 64 of 4905, or 1.30%, which on 750 test compounds predicts about ten matches; the chance
of none is 5·10⁻⁵, so the zero is deliberate. Two readings survive and they point opposite ways.
*A-strong*: the test was assembled to avoid publicly known chemistry — then it is depleted of
public compounds and the model this document rests on, test = anchors plus their analogues from
the same pool, is weakened. *B*: the released external file was deduplicated against the blinded
test before publication — an ordinary act that says nothing about the test.

A-strong predicts the test must be **farther** from public chemistry than our training set is,
since that is the property it selects on. Measured: the test is **closer**, by +0.014 of median,
95% interval [+0.010, +0.021] over 2000 draws, none of them negative. Both sides are resampled
and the training side by Butina cluster, because holding one side fixed is the defect of item 55.
A-strong is refuted.

A third reading is not refuted and should not be. *A-weak* — the test was cleaned of exact
public matches only — predicts precisely what is observed, no exact overlap and an untouched
neighbour distribution, and is observationally identical to B. It is also harmless: under it at
most 1.3% of the test was removed, against the twofold depletion in basic amines that carries
the CYP2D6 shift. The orders of magnitude are not comparable and the δ analysis stands either way.

**63. Why CYP2C9 alone is harmed by the external data, and the answer is not their labels.**
`k11_exttransfer.py`. Item 61 offered an explanation that does not fit: that CYP2C9's paired
offset rests on six shared compounds and is unreliable. But CYP2C9 is harmed in the arm where
**no offset is applied at all**, and it has 2506 external labels, second most of the four.
Neither thin estimation nor scarce data can be the cause.

One run separates their labels from our merging — train on the external rows alone and predict
our compounds. Nothing is merged, so whatever appears belongs to their data.

    enzyme     ext   ST-RAE   our rho   theirs   ratio   effect
    CYP1A2    1343    0.973     0.496    0.351    0.71   -0.026
    CYP2C9    2506    1.236     0.597    0.397    0.66   +0.022
    CYP2D6    2540    1.237     0.403    0.340    0.84   -0.036
    CYP3A4    4773    1.247     0.765    0.535    0.70   -0.011

**CYP2C9 does not stand out.** By raw rank correlation it is second of the four, and by ST-RAE
it ties CYP2D6 and beats CYP3A4. The hypothesis that their CYP2C9 labels pool incompatible probe
substrates — plausible chemistry, since CYP2C9 IC50 depends strongly on whether the probe is
diclofenac or tolbutamide — is not supported by this measurement, and is not testable from the
published file at all, which carries SMILES and four label columns and no assay metadata.

Raw rank correlation is not comparable across enzymes, which are measured on different compounds
with different predictability. Dividing by our own out-of-fold rank on the same compounds removes
that, and the resulting **ratio orders the four enzymes exactly as the effect does**: CYP2D6
transfers best and gains most, CYP2C9 transfers worst and is the only one that loses. The
break-even ratio lies between 0.66 and 0.70. Four points, and a perfect ordering arises by chance
with probability 1/24 — so this is a hypothesis carrying a prediction, not a result: raise the
transfer on CYP2C9 and the sign must flip.

**64. One descriptor overflows, on the external compounds and nowhere else.** Both external runs
printed `overflow encountered in cast`. The cause is `Ipc`, which grows exponentially with
molecule size: our 4905 compounds reach 5.1·10¹⁴, the external 8004 reach **1.5·10³⁶**, and one
exceeds float32 outright. For the boosting this is harmless — bins are quantiles and the infinity
lands in the top one. For the neural trunk it is not: the standardising variance is computed in
float32, squaring 10³⁶ overflows, and the column is divided by infinity and **silently becomes
zero**. Both arms of the trunk experiment lose it equally, so the channel's contribution — the
quantity that experiment exists to measure — is unaffected; what it changes is the account of the
bookkeeping term, from "the normalisation moved" to "one feature was deleted".

The submission path was checked and is clean: the blinded test's largest `Ipc` is 4.1·10⁸, nine
orders below overflow. No change was made — guarding the descriptor in `feats.py` would move
published numbers for no gain.

**65. The gap between the two label sets is not uniform, and one indicator column beats every
offset.** `src/ablsrc.py`. Item 61 treated the gap as a bias and subtracted it three ways, the
crudest of them harmfully. But the gap is a **selection** effect — ChEMBL holds what people
chose to publish, which is what worked — and selection does not act uniformly across chemical
space. Subtracting one number per enzyme repeats the error item 51 diagnosed for the test set:
a shift is not the whole shape, and the test is widened as well as moved.

One extra column instead, one for external rows and zero for ours, trained on the union and
predicted with the column set to zero. The learner then decides region by region how much the
external labels say about ours.

    seed 0                    1A2      2C9      2D6      3A4    macro
    no external rows        0.8786   0.6908   0.9803   0.5194   0.7673
    external as they are    0.8526   0.7124   0.9444   0.5087   0.7545
    source indicator        0.8548   0.6905   0.9386   0.5179   0.7505

The comparison was pre-registered because this arm is a strict generalisation of both
extremes: split on the indicator at the root and it reproduces the first row, never split on
it and it reproduces the second. Landing between them would have meant the indicator bought
nothing a constant did not. It landed **below both**, so the borrowing genuinely has to differ
across regions.

The headline is the second column, not the macro. **CYP2C9 flips sign**, from +0.022 to
−0.0003, and is no longer the one enzyme the external data harms. That is the prediction item
63 attached to its own hypothesis — raise the transfer on CYP2C9 and the sign must flip — made
before this arm was run and satisfied on the first attempt. CYP2D6 improves further as well,
0.9444 to 0.9386. CYP1A2 and CYP3A4 give back a little against the raw arm and stay ahead of
the control.

One seed. The rule here is four, and this needs them before it enters the document as a number
rather than as a direction.

**66. The auxiliary head works, and it works on the enzyme the merge breaks.** `src/trunkext.py`.
The joint-likelihood trunk was built, measured and shelved because the screening channel it was
given carries little about pIC50. The machinery was never the problem; the channel was. The
external CYP labels are a channel with information in it, and feeding them to the auxiliary head
instead removes the question that the boosting route has to answer — how the two label scales
relate — rather than answering it. External rows enter with fold index −1, so they sit in every
training fold and no held-out one; their pIC50 entries are NaN and ours are NaN in the auxiliary
block, so neither head sees the other's rows. `trunk.py` is imported unchanged and its published
numbers are untouched.

The first attempt at a control was wrong and the run refuted it. It looked obvious that λ = 0 here
must reproduce `trunk_twohead.json` exactly, since external rows carry no pIC50 and so contribute
no gradient to the primary head. Maximum absolute difference: 6.49. The cause is not a leak in the
fold logic — feature standardisation is computed over the training rows and there are now 12909
instead of 4905, and the step count rises because epochs are fixed and batches are not. External
rows change the model without contributing a single gradient through either head. Item 64 names
one concrete part of it: `Ipc` overflows and the column is silently deleted.

So the effect decomposes into two terms, and conflating them would have made the whole measurement
noise:

    seed    baseline   lambda 0   lambda 3    channel    bookkeeping
       0      0.7672     0.7748     0.7684    -0.0064        +0.0076
       1      0.7837     0.7843     0.7748    -0.0095        +0.0006
       2      0.8099     0.7839     0.7627    -0.0212        -0.0260
       3      0.7623     0.7938     0.7606    -0.0332        +0.0315
    mean      0.7808     0.7842     0.7666    -0.0176        +0.0034

**The channel is worth −0.0176 and its sign holds on all four seeds.** The bookkeeping term — what
the rows do through normalisation and step count alone, with no information transfer — averages
near zero but ranges over 0.058, three times the effect it sits on top of. Any single-seed reading
against the published baseline would have measured mostly that.

Per enzyme the channel and the boosting route do not overlap. The channel gains most on CYP1A2
(−0.039) and CYP2D6 (−0.022); appending rows to the boosting gains most on CYP2D6 (−0.036) and
CYP1A2 (−0.026). The difference that matters is CYP2C9: the boosting route **harms** it by +0.022
and the auxiliary head helps it by −0.008. That is precisely the design argument — each source
keeps its own output and no assumption about the relationship between scales is needed — and the
enzyme that breaks under forced merging is the one it rescues.

An earlier reading of two seeds put the channel at −0.008 and concluded that a single indicator
column beat the whole shared latent. On four seeds it does not: −0.0176 against the indicator's
−0.0168. The comparison is still not settled, because every boosting arm here has one seed and
this has four.

**67. Item 63's hypothesis was never distinguishable from its rival, and the caveat named the
wrong risk.** Item 63 found that the external data's per-enzyme effect is ordered exactly by how
well the external labels transfer, and guarded the finding with "four points, a perfect ordering
arises by chance with probability 1/24". An outside reading proposed a different explanation for
the same ordering — not transfer but **identifiability of the source offset**, since the number of
compounds shared between the two label sets is 15 / 6 / 41 / 9 and CYP2C9 has the fewest.

Both explanations were checked against the same four numbers, and they are **rank-identical**.
Transfer ratio 0.71 / 0.66 / 0.84 / 0.70 has ranks [3, 1, 4, 2]; pair count 15 / 6 / 41 / 9 has
ranks [3, 1, 4, 2]. Both give ρ = −1.00 against the effect of appending rows. On this data no
measurement can prefer one, and the guard should have asked how many hypotheses produce the same
ordering rather than how often chance produces an ordering at all.

Worse for both: **neither explains the indicator column's own contribution**, which is what
actually flipped CYP2C9. Against the column's per-enzyme gain both give ρ = +0.20. And the
mechanism named in item 63 is not the mechanism that operated — an indicator column does not raise
transfer, it lets the model express a source correction — so the prediction "raise the transfer and
the sign flips" was satisfied by something else. Affirming the consequent, and it was written as a
confirmation.

**68. The pretrained encoder does add information; item 61's headline was wrong.** That item
concluded the representation was never the bottleneck, from a table in which the embedding only ever
**replaced** our features. The missing row is concatenation, and the reason it is the row that
decides is specific to this checkpoint: `rdkit2d` was pretrained to predict the very RDKit
descriptors our DESC block contains, which makes it the honest choice for the substitution question
and the **least** favourable one for the complementarity question. An encoder trained to reproduce
what we already have is the one least able to add to it.

    FP+DESC+MECH (control)   0.7673
    EMB (substitution)       0.7784
    EMB+MECH                 0.7765
    FP+DESC+MECH+EMB         0.7589

Concatenated it is worth **−0.0084** while as a replacement it loses 0.011. It helps CYP1A2 (−0.022)
and CYP2D6 (−0.039) and harms CYP2C9 (+0.009) and CYP3A4 (+0.019) — the same per-enzyme pattern the
external data shows. One seed.

**69. The scale gap is selection on two enzymes, quantity on a third, and neither on the fourth.**
`src/ablsrc.py`. If the gap between the two label sets is selection — ChEMBL holds what people chose
to publish — then the correction belongs on the WEIGHT, not the label: reweighting moves the label
distribution without altering a single label, where subtracting an offset corrupts every one of
them. Exponential tilting is the minimum-relative-entropy way to do it and was already written in
`src/reweight.py` for the delta work.

Two targets, because they answer different questions. Tilting to our own mean removes the whole
marginal gap, which is the quantity that as a subtraction was catastrophic. Tilting to the **paired**
offset removes only the source effect measured at fixed chemistry, leaving the enrichment in actives
intact as information.

    no external rows          0.7673
    external as they are      0.7545
    source indicator          0.7505
    tilt to our mean          0.7452
    tilt to the paired offset 0.7363

**Tilting to the paired offset is the best result this repository has produced**, −0.0310 against the
control, twice the indicator and six times the entire budget of the model-choice question.

The obvious alternative explanation is that tilting simply uses less external data — effective sample
size falls from 4773 to 106 on CYP3A4 under the marginal target. The control is to permute the same
weights across rows, which preserves effective sample size exactly and destroys the correlation with
the label. It separates the two cleanly:

    change from the raw arm      shape      quantity
    CYP1A2                      +0.0013      -0.0143
    CYP2C9                      -0.0363      +0.0028
    CYP2D6                      -0.0449      +0.0170
    CYP3A4                      +0.0268      +0.0103
    macro                       -0.0133      +0.0039

The two mechanisms are nearly orthogonal and distributed differently across enzymes. On CYP1A2 the
whole gain is **quantity** and the shape does nothing. On CYP2C9 and CYP2D6 the **shape** does
everything and reducing the data actively hurts. On CYP3A4 the marginal tilt is worse than random
down-weighting of the same strength, which is what the paired target then repairs — effective sample
goes from 106 to 3010 and the enzyme moves from 0.5458 to 0.5020. In the macro the two mechanisms
partly cancel, and without the permutation control the whole gain would have been credited to one.

**70. The noise floor at fixed seed is 0.007, and it was never measured before.** A falsifiable
prediction accompanied item 69: the paired and marginal targets coincide on CYP2D6 (−0.60 against
−0.59), so the two arms should differ there by nothing. They differ by 0.0071 — from a change of 0.01
in the tilt target, propagated through different split choices. On CYP2C9 the targets differ by 0.38
and the arms differ by 0.0053, less than that.

So 0.007 is this pipeline's chaotic sensitivity at a fixed seed and fixed folds, and differences below
it are not interpretable even before seed variation is considered. This is a smaller number than the
0.016 that item 7 measures for a change of split seed, and a larger one than several comparisons
recorded earlier in this log were resting on.
