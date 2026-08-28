# Verification

Twenty-two scripts in four groups. `f*` was a sweep over everything that had been computed
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
    current count is twenty-two and both files say so.)
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
