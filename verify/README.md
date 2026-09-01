# Verification

## Табло

Этот файл --- журнал дефектов, и читается он соответственно: верная идея получает один пункт,
неверная получает три (выдвижение, опровержение, поправка к опровержению). Соотношение при
чтении выходит три к одному в пользу провалов при положительном итоге. Табло существует, чтобы
состояние проекта не приходилось складывать в голове из ста шестидесяти пунктов.

**Где мы.** Ранг --- Спирмен с истиной, пара --- ST-RAE после аффинной пары. Больше ранг лучше,
меньше пара лучше.

    конфигурация                                 пара      ранг   прирост   сидов   пункт
    база FP+DESC+MECH, поферментно, HistGB     0.7150    0.5651         —       1     эталон
    ансамбль из четырёх членов                 0.6819    0.6009   +0.0358       4       120
    ансамбль из пяти (ЧТО ПОДАЁТСЯ СЕЙЧАС)     0.6824    0.6063   +0.0412       4   120, 121
    пять, мёртвая зона в четырёх членах        0.6653    0.6230   +0.0579       4       164
                                                                                 (нижняя оценка)

**Итого траектория: +0.058 ранга и -0.050 пары от базовой модели.** Шум одного счёта лидерборда
на 750 молекулах --- 0.08 пары (пункт 147), так что по паре проект пока внутри него; по рангу
--- вне, но ранг на лидерборде не показывают.

**Что стоит в конвейере и сколько стоит.**

    вклад                                  прирост ранга   сидов   пункт   статус
    мёртвая зона во всех членах                  +0.0167       4     164   в подаче нет
    ствол пятым членом                           +0.0054       4     120   в подаче есть
    механистический блок                         +0.0163       4      81   в подаче есть
                        он же в режиме теста     +0.0313       4     119
    пулирование контрастом                       +0.0141       4  84,132   в подаче есть
    GP и гребневая как члены                  своя ошибка      4  92,100   в подаче есть
    скрининг как мишень, поферментно             +0.0245       1     158   НЕ в подаче
    панель NCGC, поферментно                  считается        —     157   НЕ в подаче
    пятьдесят битов = 80 % фингерпринта                —       1 150,156   интерпретация

**Шумовые полы.** Макро 0.007 при фиксированном сиде (пункт 70); посидовый разброс макро 0.016
(f3). Поферментные полы **больше макро** и до сих пор нигде не были записаны (пункт 165):

    фермент   sd по сидам при неизменной руке
    CYP1A2                             0.0061
    CYP2C9                             0.0071
    CYP2D6                             0.0049
    CYP3A4                             0.0033
    МАКРО                              0.0036

Макро усредняет четыре фермента и потому тише каждого из них. **Поферментное заявление нельзя
мерить макро-полом.**

**Три главных открытых вопроса.**

  почему пулирование выигрывает на HistGB и разворачивается на обычных деревьях (158, 166);
  переносится ли что-либо из этого в тестовый режим --- проверить нечем (123, 129, 147);
  скрининг поферментно даёт +0.0245 на одном сиде и ждёт остальных.

**Что закрыто и переоткрывать не надо.** Постобработка сверх аффинной пары (77, 128 --- потолок
0.0076), предобученные представления (61, 117, 154 --- три чекпойнта), FCFP (101), kNN (107),
краевая регрессия (114), точное байесово действие (126), серийный слой (133), координата фермента
(154), logD и LipE (154), хи-квадрат-DRO (154), ChEMBL как внешний источник (77 --- все три
способа обращения), CYP2C19 как пятая изоформа (153).

**Не запущено, а не закрыто.** Квантовый блок: признаки ПОСЧИТАНЫ (`data/quantum.npz`, 4905x10,
homo/lumo/gap/dipole/q_basicN/q_aromN_min/fukui_minus/cone_free/n_arom_N/has_donor, ноль NaN), а
абляция никогда не запускалась. Отдельно `src/quantum.py` падает на `xtb-python`, которого нет в
реестре, но файл от более раннего прогона на месте. Это не результат ни в какую сторону.


Forty-three scripts in four groups. `f*` was a sweep over everything that had been computed
and written by that point; `g*` answers four questions raised against the document; `h*`
tests two claims the document made about geometry and about reactivity; `k*` began as a
group about post-hoc rescaling of the predictions and about how far the test set sits from
the training distribution, and `k20` onward is about a second question that grew out of the
first: whether the regime our cross-validation runs in is the regime the test will be scored
in, and which of this log's conclusions depend on the answer. All of them expect the data in `data/` and the organisers' metric code
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
| `k10_strat2d6.py` | the kernel estimate of delta, split along the one axis known to have moved | ~5 min |
| `k11_exttransfer.py` | train on the external labels alone: which enzyme fails to transfer, and how | ~15 min |
| `k12_extneighbors.py` | how close the external compounds get to ours, and what zero overlap with the test means | ~4 min |
| `k13_channels.py` | how many compounds carry all three observation channels, and whether the channels are distinct | ~1 min |
| `k14_design.py` | recover the organisers' rule for choosing which compounds get a dose-response curve | ~2 min |
| `k15_pooldelta.py` | re-estimate the test shift on the pooled base | ~6 min |
| `k16_modelspread.py` | how far apart three models put the same test shift | ~8 min |
| `k17_ensdelta.py` | the shift estimated on the model that is actually submitted | ~6 min |
| `k18_nbspace.py` | in which space is a neighbour informative — measured with no model at all | ~5 min |
| `k19_ens3delta.py` | the shift on the three-member ensemble, and whether a fourth family widens the spread | ~7 min |
| `k20_strat.py` | every model's rank, measured in the regime the test set actually sits in | ~4 min |
| `k21_borda.py` | average the members' ranks instead of their values | ~1 min |
| `k22_layerboot.py` | is the layer effect real, and does gating the ensemble on similarity pay | ~8 min |
| `k23_tilt.py` | reweight the ensemble toward pooling with one parameter, chosen out of sample | ~6 min |
| `k24_visible.py` | why pooling's advantage grows exactly where the test set sits | ~5 min |
| `k25_reweight.py` | re-score every saved ablation under weights matching the test set's regime | ~10 min |
| `k26_screen.py` | preconditions for the screening head: informativeness, novelty, and what it corrects | ~3 min |
| `k27_trunkens.py` | the trunk as a fifth ensemble member, with the lambda-zero control | ~4 min |

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

*Later.* The four-seed run was started and stopped after six of sixteen cells, because items 69
and the loss experiment produced arms worth 0.031 and 0.046 against this arm's 0.013, and it was
holding half the machine for eleven more hours to validate the weakest of the three. What it
reached is kept here rather than discarded: the control gives 0.7673 and 0.7690 on seeds 0 and 1,
the raw external arm 0.7545 and 0.7482, so the effect is **−0.0128 and −0.0208** and its sign
holds on both. Two seeds, and the second is the larger. The rule of four now belongs to the arms
that superseded this one.

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

The asymmetry is worth naming rather than assuming harmless. The trunk with external rows ran on
2294 live columns and every boosting arm on 2295, and the trunk's own published baseline keeps all
2295 because our training set tops out at 5.1·10¹⁴ and the sum of squares stays inside float32.
The lost column is a degenerate one and its loss sits entirely inside the bookkeeping term, but a
comparison between two families should not carry an unremarked difference in what they were given.

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

**68. The pretrained encoder does add information; item 61's headline was wrong.**
*(Withdrawn by item 117. The gain below is raw ST-RAE, which item 77 showed the affine pair
rewrites; after the pair it is +0.011 and by rank it is -0.0075. Item 61's headline stands and
this item does not. Kept in place because item 77's table counts it among the five.)* That item
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

**71. The three-channel proposal passes its size gate by two orders of magnitude, and fails a
different one.** `k13_channels.py`. A mechanistic proposal separates two latents the single pIC50
target confounds — affinity, and turnover into something reactive — using direct inhibition,
time-dependent inhibition and Emax as three views. The objection raised against it was that the
compounds carrying all three might number in the dozens, which would end it before any modelling.

That objection was a guess, and it was wrong by a factor of a hundred. The channels are essentially
co-extensive: every compound with a direct pIC50 also has a TDI curve and an Emax, giving 1412 /
1285 / 1493 / 2334 triples, **6524 compound-enzyme observations**, the whole labelled set. On
CYP2D6 and CYP3A4 the binary `is_TDI` covers them too, so there are four channels rather than
three; on CYP1A2 and CYP2C9 that label does not exist at all.

The constraint is elsewhere and it is real. Direct and TDI pIC50 are the same molecule under two
incubation protocols, and they correlate at 0.906 to 0.990. Everything the second channel says
about turnover lives in the residual between them, whose spread is 0.17 / 0.28 / 0.45 / 0.33 of the
signal's own. On CYP1A2 a turnover latent would have to be estimated from a sixth of the variance.
Emax is a genuinely separate axis, |r| about 0.4 against affinity, so a third dimension does exist.

One anomaly to explain before building rather than after: on CYP2D6 the sign of the Emax
relationship is **reversed**, +0.770 against −0.394 / −0.407 / −0.462 elsewhere, and the same flip
appears in the TDI arm. It is the same enzyme that is anomalous in the salt bridge, in the sign of
its shift, and in how well external data transfers to it.

**72. Emax has no dynamic range, which kills a proposal and one of my own claims.** `k13_channels.py`.
Item 71 recorded that Emax is a genuinely separate axis because it correlates with affinity at |r|
about 0.4. That was wrong, and wrong in the way this log keeps finding: a rank correlation says
nothing about whether there is anything to measure. The range does.

    enzyme      min      1%     50%     99%     max     IQR   above -0.5
    CYP1A2    -1.15   -1.10   -0.99   -0.82   -0.46   0.072       0.071%
    CYP2C9    -1.27   -1.19   -1.02   -0.73   -0.50   0.106       0.078%
    CYP2D6    -1.08   -1.06   -1.03   -0.91   -0.83   0.048       0.000%
    CYP3A4    -1.06   -1.04   -0.98   -0.84   -0.37   0.020       0.214%

**Every compound in the set inhibits essentially completely.** Ninety-nine percent of the mass sits
below −0.73 and not one of 6524 observations exceeds zero. The |r| ≈ 0.4 correlations are ordering
inside a band 0.05 wide. A latent for mechanism cannot rest on a channel that is nearly constant,
so the third view in the affinity-and-turnover proposal is not available, and that is a harder
obstacle than the collinearity of the first two.

It also refutes a prediction derived from one-site against multi-site binding: Emax should have been
bimodal on CYP3A4, where a large cavity with several sub-sites would leave partial inhibitors, and
unimodal at full inhibition on CYP2D6, where a single small site with a salt bridge would be
occluded completely. A Gaussian mixture prefers two components everywhere, but the components are
0.03 to 0.13 apart — it is splitting one narrow blob. And the direction is reversed: CYP3A4 has the
**narrowest** distribution of the four, IQR 0.020 against CYP2C9's 0.106.

**73. Rank agreement on four points does not single out a hypothesis, and this is the third time.**
Item 67 found that two explanations of the same ordering are rank-identical, and recorded it as "no
measurement here can prefer one", which is too strong: nothing about rank ordering can prefer one,
and rank ordering on four points is the weakest evidence available. The pattern has now failed three
times — the screening correlation against ST-RAE change, and transfer against pair count.

The rule to carry: **agreement of ranks at n = 4 does not select a hypothesis, it only fails to
refute it.** Hypotheses of this kind separate under intervention, not under more correlations. The
two candidates here make different predictions about what happens if more or better-matched external
labels are added on CYP2C9 — a transfer deficit predicts improvement, an identifiability deficit
predicts none, because the indicator column has already extracted what there was. Those are two
runs, not two correlations.

**74. Pre-registration: three orderings for where the L1 gain should land.** Switching the learner's
loss from squared to absolute error is worth **−0.0455** of macro on seed 0, 0.7218 against 0.7673.
That is larger than the external data with its best correction (0.031) and nine times the whole
budget of the model-choice question (0.005). The number is not in doubt; the mechanism is, and three
candidates make different predictions about which enzyme gains most. This is written before the
per-enzyme table exists.

*Outliers.* Squared error chases extreme residuals and absolute error does not. Measured: the worst
one percent of compounds carries 5.1 to 8.2 percent of ST-RAE, five to eight times uniform, and the
largest absolute error is 3.9. The −360 prediction cited earlier for this belongs to the neural
trunk, not the boosting, which cannot leave the range of its leaf values. Predicts the maximum on
**CYP2D6**, where the anomalies have lived.

*The inactive mass.* A large share of compounds sits low, 23.7 / 42.6 / 32.1 / 58.4 percent below
4.5, so for a typical compound the conditional distribution is one-sided and mean and median part
company. Note this is not censoring in the strict sense — the minimum label is 1.91 to 2.10, values
continue down to it, and only 0.1 percent sit exactly at the floor. Predicts **3A4 > 2C9 > 2D6 >
1A2**.

*The mean-minus-median gap.* Measured directly on the residuals: −0.049 / −0.027 / −0.009 / −0.010.
Predicts **1A2 > 2C9 > 3A4 > 2D6** — nearly the reverse of the previous one, with CYP1A2 first
rather than last.

None of the three currently accounts for the size. The mean-median gap is 0.05 of pIC50 at most
against band half-widths of 0.13 to 0.28, so it sits inside the region the metric does not penalise
at all; the outlier concentration is real but moderate. A fourth possibility is not about the
estimand but about training: with absolute error the gradient is a sign, so hard compounds cannot
buy a disproportionate share of the splits, and ST-RAE additionally makes error inside the band free,
so squared error spends capacity on differences that are not scored. That one predicts no particular
enzyme and is not separable with the losses this scikit-learn offers — there is no Huber in
`HistGradientBoostingRegressor`.

**75. The loss function is worth 0.046, and all three pre-registered mechanisms failed.**
`src/abloss.py`.

    seed 0             1A2      2C9      2D6      3A4    macro
    L2 on the label  0.8786   0.6908   0.9803   0.5194   0.7673
    L1 on the label  0.8104   0.6802   0.8872   0.5095   0.7218
    L1 on the bounds 0.8098   0.6813   0.9015   0.4906   0.7208

The control reproduces. Switching to absolute error is worth **−0.0455**, larger than the external
data with its best correction and nine times the model-choice budget, from one keyword.

The band-target arm confirms a prediction made in item 74 from the measured band asymmetry: pooling
`lo` and `hi` and fitting the median is worth **0.0010** against fitting L1 on the point label, six
times below the noise floor. Section 10's action is right and minimising it directly buys nothing,
because the bands are nearly symmetric so the pooled median sits on the point label. The whole gain
is L1 against L2. Without that middle arm the result would have been credited to the theorem.

The per-enzyme ordering is 2D6 0.0931 > 1A2 0.0682 > 2C9 0.0106 > 3A4 0.0099, and it refutes all
three candidates of item 74. The inactive-mass account predicted 3A4 first and 1A2 last and gets
ρ = −0.8, close to exactly inverted. The mean-minus-median account predicted 1A2 first, ρ = −0.2.
The outlier account predicted the maximum on CYP2D6 and got it — but the quantity it rests on, the
share of ST-RAE in the worst one percent, orders the enzymes 2C9 > 3A4 > 1A2 > 2D6 with **CYP2D6
last**, so the prediction was satisfied against its own evidence.

What does order the gain is how badly the model was doing: baseline ST-RAE 0.9803 / 0.8786 / 0.6908
/ 0.5194 against gains 0.0931 / 0.0682 / 0.0106 / 0.0099, a perfect match. By item 73's own rule
that selects nothing — it is rank agreement on four points, the weakest evidence there is — so it
is recorded as a description and not a mechanism. The mechanism is open.

**76. Monotonic constraints do nothing, and the way they do nothing is the finding.**
`src/ablmono.py`. Four features whose sign biochemistry knows, one to two per enzyme, constrained
per enzyme.

    seed 0             1A2      2C9      2D6      3A4    macro    tails
    free             0.8786   0.6908   0.9803   0.5194   0.7673   0.7007
    correct signs    0.8735   0.6942   0.9803   0.5155   0.7659   0.7026
    flipped signs    0.8765   0.6941   0.9803   0.5391   0.7725   0.7049

Correct signs gain 0.0014 and flipped signs lose 0.0052, both under the 0.007 noise floor. The
tails, where the whole extrapolation argument lives, get monotonically worse rather than better.
There is no effect to interpret at the macro level.

One enzyme is different and it is the wrong one. On **CYP2D6 the predictions are bit-identical
across all three arms** — maximum absolute difference exactly 0.00e+00 in both directions — while
the other three move by 0.29 to 1.09. The constrained feature there is `frac_prot_74`, the
protonated fraction, which is the feature this document's entire CYP2D6 story rests on: the salt
bridge, the stratified kernel, the composition half of the shift. The boosting never splits on it.

The dull explanation is ruled out: on CYP2D6 that feature is the most variable of the four, sd 0.434
against 0.29 to 0.30, with 551 distinct values, and no single fingerprint bit correlates above
0.405. The remaining explanation is joint redundancy — a ridge on the fingerprint block alone
reproduces it with R² = 0.653 out of fold, the highest of the four enzymes — so the tree can reach
the same chemistry through the bits and never needs the aggregate. That is partial: 0.653 against
0.50 to 0.62 elsewhere is not sharp enough to explain an exactly zero difference where the others
are nonzero. Recorded as explained in part.

Nothing here contradicts the document. Stratifying by a feature and splitting on it are different
uses, and the shift analysis needs the first. But the feature that carries the chemistry in our
analysis carries nothing in our model, and that is worth knowing before any more weight is put on it.

**77. Every gain measured today disappears under the post-processing the submission actually uses.**
The question was narrow — recompute the δ rule on the L1 predictions, since the rule in `submit.py`
was tuned against the residuals of a model that L1 replaces. The rule barely moved: under the mean
criterion CYP1A2 goes +0.1 to 0.0, CYP2C9 +0.3 to +0.4, CYP2D6 and CYP3A4 unchanged, and the value
of per-enzyme fitting is the same, +0.0486 against +0.0468. The shift work is not absorbed by the
loss change, which was the risk worth checking.

But the levels inverted, and following that inverted the whole session. Every number in this
repository's ablation tables is a raw out-of-fold ST-RAE, and the submission does not submit raw
predictions — it fits an affine pair per fold and applies it. Refitting that pair out of fold on
each arm gives:

    arm                       raw    +affine    gain raw   gain after
    L2 baseline            0.7673     0.7150     +0.0000      +0.0000
    L1 instead of L2       0.7218     0.7131     -0.0455      -0.0019
    source indicator       0.7504     0.7117     -0.0168      -0.0033
    tilt to paired offset  0.7363     0.7156     -0.0310      +0.0006
    external as they are   0.7545     0.7218     -0.0128      +0.0068
    encoder concatenated   0.7589     0.7263     -0.0084      +0.0113

**Not one survives.** All five fall under the 0.007 noise floor of item 70, and three of the five
change sign to harmful. The largest finding of the day, worth 0.046 raw, is worth 0.002 after — and
the mechanism is plain: the affine pair is a shrinkage fitted directly to the metric, per fold, so
anything that reduces prediction variance or removes a systematic offset is something it already
does. On the L2 baseline that pair is worth 0.0523; on the L1 predictions only 0.0087, because L1
had already taken most of what there was.

`fit_apply` was checked for contamination first: it fits on the training folds and applies to the
held-out one, so the comparison is honest.

What this costs is not five results but a method. **The ablation grid measures the wrong quantity.**
Comparing feature sets, learners and data sources on raw predictions overstates differences that the
submission pipeline erases, and the document's headline — model choice 0.005 against post-processing
0.051 — understates its own point: the post-processing does not add to the other improvements, it
substitutes for them.

Three qualifications, and none rescues the arms. This is seed 0, and at the noise floor no ordering
among the six is resolvable — but five independent arms collapsing into a band of 0.013 is a much
stronger pattern than any single comparison. The affine pair here is fitted at δ = 0; under the
posterior over δ the same reversal holds, L1 with its own δ scoring 0.7764 against the L2 baseline's
0.7647. And redundancy on our own label marginal is not proof of redundancy on the test's shifted
one — but that cannot be checked without test labels, and assuming it in our favour is the error
this log exists to catch.

*Four seeds, on the largest arm.* The single-seed objection was the obvious one, so the loss arm
was run on all four.

    seed    L2 raw   L2 affine    L1 raw   L1 affine    gain raw   gain after
       0    0.7673      0.7150    0.7218      0.7131     -0.0455      -0.0019
       1    0.7690      0.7183    0.7242      0.7156     -0.0449      -0.0027
       2    0.7661      0.7164    0.7199      0.7118     -0.0462      -0.0046
       3    0.7644      0.7124    0.7234      0.7139     -0.0410      +0.0014
    mean    0.7667      0.7155    0.7223      0.7136     -0.0444      -0.0019

The raw gain holds its sign on all four seeds and spans 0.041 to 0.046, so by this repository's own
standard it is a result. After the affine pair it is 0.0019 and the sign **breaks on seed 3**, so by
the same standard it is not one. The rule that has governed every comparison here for two months
answers this question by itself.

Incidentally the affine pair is worth 0.0512 on the L2 baseline averaged over four seeds, which is
the figure the document published by another route.

Practical consequence: **nothing in `submit.py` changes**, on any of today's findings.

**78. The two corrections do not add, and the full square shows why.** The question was whether the
loss change and the paired tilt combine. The missing cell was run — `src/ablsrc.py` gained a
`--loss` switch so the external arms can use absolute error — and the square closes.

    seed 0                      raw    +affine
    L2                       0.7673     0.7150
    L2 + paired tilt         0.7363     0.7156
    L1                       0.7218     0.7131
    L1 + paired tilt         0.7225     0.7238

**They do not add.** The tilt is worth −0.0310 on top of L2 and **+0.0007** on top of L1: L1 has
already taken everything the tilt was buying. Additivity would have predicted 0.6908; the observed
value is 0.7225, which is L1 alone.

With the affine pair applied the four collapse into a band of 0.0107, of which 0.007 is noise, and
the worst of the four is **both corrections together**. So there are three interventions here — a
metric-aligned loss, a reweighting of the external labels, and the affine shrinkage — and they are
substitutes rather than complements. Any one reaches about 0.713 and applying more than one gains
nothing or costs something.

That is the same statement as item 77 seen from the other side. The affine pair is a shrinkage
fitted to the metric per fold; L1 is a shrinkage built into the fit; the tilt is a reweighting that
happens to reduce the same variance. The raw column spans 0.046 across these four and the
post-processed column spans 0.011. Whatever we were measuring in the raw column was largely the
absence of a correction we always apply.

**79. Item 77 had an exception and missed it: the trunk's screening channel survives.** The
table in item 77 put five arms through the affine pair and found nothing left. It did not
include the sixth — the joint-likelihood trunk — because that arm lives in another file and was
measured raw, at $-0.0176$. Put through the same pipeline, with the same clipping the
document's own model comparison uses and the same affine pair on all three:

    seed    trunk lam=0   trunk lam=3   boosting    channel   trunk3 - boost
       0         0.7397        0.7182     0.7150    -0.0214          +0.0032
       1         0.7446        0.7142     0.7183    -0.0305          -0.0042
       2         0.7376        0.7145     0.7164    -0.0231          -0.0019
       3         0.7435        0.7128     0.7124    -0.0307          +0.0003
    mean         0.7414        0.7149     0.7155    -0.0264          -0.0006

The channel is worth **−0.0264 after post-processing**, the sign holds on all four seeds,
t = −10.87, p = 0.002. It is the only intervention in this repository that survives. Clipping is
not optional here and it is not a thumb on the scale: without it seed 0 gives 0.9192, which is
the single compound of item 42 predicted at −360, and the document's existing comparison already
clips for that reason.

The second half does not follow, though. The trunk with the channel does **not** beat the
boosting: 0.7149 against 0.7155, the sign alternates across seeds, t = −0.39, p = 0.72. Parity,
not superiority. So the decision to shelve the trunk stands; what was wrong was writing the
channel off as carrying nothing.

**80. Why post-processing absorbs some interventions and not others, and it is one line.** The
affine pair, and isotonic calibration too, are **monotone**. A monotone map can move scale and
location and cannot change the order of predictions. So an intervention survives post-processing
if and only if it adds **rank** information, and that is measurable directly rather than
inferred:

    change in Spearman with the truth      delta      t        p    signs
    channel in the trunk (lam 3 - lam 0)  +0.0350   13.83    0.001    ++++
    L1 instead of L2                      +0.0016    0.97    0.405    -+++
    trunk lam=3 against boosting          +0.0016    0.53    0.632    -+++

    rank itself: trunk lam=0  0.5297   trunk lam=3  0.5646   boosting  0.5630   L1  0.5646

That closes three questions at once. The loss change collapses because it adds no rank — and
note the collapse is not, as first supposed, because L1 merely rescaled: Spearman between L1 and
L2 predictions is 0.82 to 0.96, so it reordered a great deal, and none of the reordering was
information. The channel survives because it adds 0.035 of rank. And the parity between the
trunk and the boosting has the same explanation as both: the channel lifts the trunk to exactly
the rank the boosting already had, 0.5646 against 0.5630, and no further.

The practical rule for every future ablation, and it costs nothing: **report the change in rank
correlation beside the raw score.** The raw column answers "did the intervention move the
predictions", which is not the question; the rank column answers "did it move them somewhere the
post-processing cannot reach", which is. A post-isotonic score is the same criterion expressed in
the metric's own units, and is the conservative bound, since isotonic spans a wider class of
monotone maps than the affine pair does.

**81. The mechanistic block survives post-processing, and the raw measure had been hiding it.**
The question asked was whether the block does anything on CYP2D6 at all, prompted by item 76:
`frac_prot_74` is never split on there, which looked like evidence that the chemistry explains
how the test differs from the training set without entering the predictor. Both arms already
existed on four seeds, so the answer cost nothing.

    block against no block, 4 seeds     raw       p    after pair       p     d rho       p
    CYP1A2                          -0.0016   0.600      -0.0016   0.372   -0.0019   0.375
    CYP2C9                          -0.0170   0.022      -0.0147   0.029   +0.0148   0.017
    CYP2D6                          -0.0151   0.088      -0.0221   0.000   +0.0454   0.002
    CYP3A4                          +0.0017   0.498      +0.0007   0.777   +0.0008   0.619
    macro                           -0.0080   0.069      -0.0094   0.006   +0.0148   0.003

The chemistry **does** enter the predictor, decisively on CYP2D6: the block adds 0.045 of rank
correlation there, more than the trunk's screening channel adds anywhere, and after the affine
pair it is worth −0.0221 with p = 0.0004.

Note the direction, which is the opposite of everything else in items 77 to 80. For every other
intervention the raw score **overstated** the gain. Here it **understates** it: raw gives
p = 0.088 on CYP2D6, not significant, and after post-processing p = 0.0004. That is exactly what
item 80's criterion predicts. The block's contribution is rank, and the raw comparison mixes it
with scale-and-location noise that the affine pair removes; strip that away and the signal comes
out cleaner and larger.

Two things this corrects. The document records "the mechanistic block's ST-RAE gain is
indistinguishable from zero, both per enzyme and on the macro" as a negative result. That was
measured raw and is wrong. And the reading of item 76 that this check was meant to test —
chemistry explains the data but not the model — is refuted; what item 76 shows is narrower than
it looked, namely that the block works on CYP2D6 through something other than the protonated
fraction the text credits. Which part is the subject of `src/ablmech.py`.

**82. Which part of the mechanistic block works, and it differs by enzyme in the way the
chemistry says it should.** `src/ablmech.py`. Item 81 showed the block earns its place; item 76
showed `frac_prot_74`, the feature the text credits, is never split on. Splitting the block into
three groups — 16 counts of nitrogen and acid functionality, 8 protonation-state features
including `frac_prot_74`, and 6 topological distances from a basic nitrogen to an aromatic ring
plus the explicit CYP2D6 pharmacophore flags — resolves both. Four seeds, scored as item 80 says
to score.

    change in rank correlation vs no block    CYP2D6       p    CYP2C9       p
    whole block (30)                         +0.0454  0.0017   +0.0148  0.0167
    protonation state only (8)               +0.0454  0.0076   +0.0097  0.0013
    functional-group counts only (16)        +0.0437  0.0130   +0.0111  0.0042
    geometry and pharmacophore only (6)      +0.0580  0.0007   +0.0030  0.5454

On **CYP2D6 the six geometric features carry it, and carry it better than all thirty**: +0.0127
of rank over the whole block, p = 0.009. The other twenty-four dilute. On **CYP2C9 geometry does
nothing at all** (p = 0.55) and the contribution comes from counts and protonation state, with no
group distinguishable from the whole.

That is the split the mechanism predicts. CYP2D6 binds through a salt bridge to Asp301 and
Glu216, which is a **geometric** constraint — a protonated nitrogen at a particular distance from
an aromatic ring — and geometry is exactly what those six features encode. CYP2C9's Arg108 binds
anions, which is a question of **what the molecule contains**, not where. A prediction of this
shape was recorded in the script's docstring before the run and holds on CYP2D6.

It also settles item 76 without contradicting it. `frac_prot_74` is not split on because the
geometric features carry the same chemistry more sharply, not because protonation is irrelevant:
the state group alone still gives +0.0454, the same as the whole block. The groups are largely
redundant with one another, and geometry is the sharpest of the three.

One honest limit. Geometry's advantage over the whole block is solid in rank (p = 0.009) and only
a trend in the metric after the affine pair (−0.0038, p = 0.090). "Six features beat thirty" is
established about rank and not about ST-RAE, and the two are not the same claim.

**83. The screening file holds 11509 unlabelled readings, and the reason nobody could use them is
the reason they exist.** An outside reading pointed out that the single-concentration file covers
4376 molecules on all four enzymes — 17504 rows, one concentration of 4.95·10⁻⁵ M, which is
exactly pC₀ = 4.305 — while only 6525 molecule-enzyme pairs have a dose-response curve. Every
count reproduces: 2964 / 3091 / 2883 / 2571 screened pairs lack a curve, 11509 in total, 1.8× the
labelled set. The four-vector is complete for every one of the 4376 molecules, which also
corrects item 71's reading of the selectivity question: 73.3% of compounds carry one enzyme *in
the curve files*, and none of that applies to the screen.

The proposal is to invert those readings into pseudo-labels through the fitted Hill calibration
of section 4, keeping only where the Fisher information is high. Two gates were measured before
building anything.

*Invertibility.* The inversion π = pC₀ − log₁₀(E/I − 1)/h needs 0 < I < E.

    enzyme    screen-only   invertible        Fisher > 0.3 max
    CYP1A2           2964    2553  86.1%       2041  68.9%
    CYP2C9           3091    2831  91.6%       2185  70.7%
    CYP2D6           2883     362  12.6%        265   9.2%
    CYP3A4           2571    2269  88.3%       1678  65.3%

**CYP2D6 collapses to nine percent**, and it is the enzyme the project is built around. The cause
is not saturation at the top but blindness at the bottom: 87.4% of its screen-only readings have
I ≤ 0, no inhibition detected, with a median of −0.298 — noise around zero. Both of the
prediction orderings offered with the proposal put CYP2D6 above CYP3A4; neither anticipated it
getting almost nothing.

*Selection, and it is severe.* Comparing the screen readings of compounds that have a curve
against those that do not:

    enzyme    log2fc with curve   without    active fraction with / without
    CYP1A2               -1.551    -0.197           88.5%  /  40.0%
    CYP2C9               -0.815    -0.602           82.5%  /  67.0%
    CYP2D6               -1.738    +0.376          100.0%  /   5.7%
    CYP3A4               -1.437    -1.932           78.4%  /  92.0%

The CYP2D6 row explains the whole structure: **curves were run according to the screen.** All of
its curve compounds are screen-active and 5.7% of the rest are, so "screened but unlabelled" is
by construction the screen-negative subset. On CYP3A4 the selection runs the other way, so there
is no single statement of the form "pseudo-labels pull the marginal down while the test pulls it
up" — the direction differs by enzyme.

That yields an objection the proposal does not carry and which we judge decisive as stated: **E
and h were fitted on compounds selected for being active**, and inverting readings from compounds
selected for being inactive applies that calibration outside the population it was estimated on.
The error enters all ~5900 pseudo-labels systematically and in the same direction. The permutation
control suggested alongside catches "the gain was extra rows rather than information"; it does not
catch a calibration biased outside its fitted range, because permuting preserves that bias exactly.

Not run. The idea survives on CYP1A2 and CYP2C9, is questionable on CYP3A4 until the reversed
selection is understood, and is dead on CYP2D6 — and any version of it needs the calibration
re-fitted or validated on screen-negative compounds first, which is a separate measurement.

**84. Pooling the four enzymes into one model survives post-processing, and the gain lands where
it was pre-registered.** `src/ablpool.py`. Every model in this repository is fitted per enzyme, so
CYP2D6 sees 1493 rows and nothing else, while the four label sets share most of their chemistry.
Stacking them into one table with a four-way enzyme indicator gives 6525 rows and lets the shared
part be learned on all of them.

    seed 0          raw     after pair      rho
    independent  0.7673         0.7150   0.5651
    pooled       0.7415         0.7062   0.5792

**−0.0088 after the affine pair**, above the 0.007 noise floor, with a rank gain of +0.0141. The
third intervention in this repository to clear the ceiling, after the trunk's screening channel
and the mechanistic block, and like both of them it clears it by moving rank.

    rank by enzyme   independent   pooled    gain
    CYP1A2                0.4957   0.5084  +0.0127
    CYP2C9                0.5972   0.6029  +0.0057
    CYP2D6                0.4027   0.4448  +0.0421
    CYP3A4                0.7646   0.7607  -0.0039

The docstring recorded "expect the gain on CYP2D6 above all" before the run, and it is there:
+0.042, level with what the mechanistic block gives on the same enzyme. CYP3A4 is the only loss,
which is what the argument predicts — it has the most labels and the best fit already, so it has
the least to borrow.

One correction to how this should be described. Item 71 established that 73.3% of compounds carry
a label for exactly one enzyme, so pooling does **not** give CYP2D6 the same molecules seen
through four assays. It gives it five thousand *different* molecules labelled elsewhere. The
mechanism is transfer across compounds, not multi-task on shared ones, and the version that would
be multi-task on shared ones — adding each compound's pre-incubation pIC50, which item 71 showed
is 100% co-located — is a separate arm and is not yet measured.

*Four seeds.* The comparison is paired on the same folds, so the relevant scatter is that of the
difference, not the 0.007 floor of item 70 which was measured for a different kind of comparison.

    seed    independent    pooled    difference    rank gain
       0         0.7150    0.7062       -0.0088      +0.0141
       1         0.7183    0.7120       -0.0063      +0.0160
       2         0.7164    0.7110       -0.0054      +0.0086
       3         0.7125    0.7063       -0.0062      +0.0159
    mean         0.7156    0.7089       -0.0067      +0.0137

Sign holds on all four seeds in both columns: −0.0067 after the pair with sd 0.0015, t = −9.06,
p = 0.0028, and +0.0137 of rank, p = 0.0043. Three interventions now clear the ceiling, at
comparable significance — the trunk's screening channel at −0.0264, the mechanistic block at
−0.0094, pooling at −0.0067.

The three are not equally useful, and the difference decides. The trunk's channel only lifts the
network to parity with the boosting, so it adds nothing to what is submitted. The block and the
pooling both live inside the boosting route, and item 85 shows they add rather than substitute.
Pooling is therefore the first thing measured in this stretch that could actually move the
submission.

**85. Pooling and the mechanistic block do not substitute for each other, and the four enzymes
fall into four different regimes.** The two interventions that clear item 77's ceiling both land
hardest on CYP2D6 and give it almost the same rank gain, which invited the question whether they
carry the same information. The 2×2 says no.

    CYP2D6, rank, seed 0     no block    block    block effect
    independent                0.3509   0.4027         +0.0518
    pooled                     0.3847   0.4448         +0.0601
    pooling effect            +0.0338  +0.0421

Together they are worth +0.0939 against +0.0856 for the sum of the separate effects, so they
complement rather than substitute: each is worth more in the presence of the other. That reads
sensibly — the block points at the CYP2D6 recognition geometry explicitly, pooling supplies more
molecules in which that geometry can be seen, and a pointer plus examples beat either alone. The
excess over additivity, +0.0083, is barely above the 0.007 floor on one seed, so "they add" is
firm and "they add with a bonus" is not.

Per enzyme there is no single answer at all:

    enzyme     block alone   pool alone   together   interaction
    CYP1A2         -0.0012      +0.0265    +0.0115       -0.0138   substitutes
    CYP2C9         +0.0138      +0.0038    +0.0195       +0.0019   independent
    CYP2D6         +0.0518      +0.0338    +0.0939       +0.0083   complements
    CYP3A4         +0.0010      -0.0044    -0.0029       +0.0005   both empty

On **CYP1A2 the block is worthless alone and harmful under pooling** — pooling by itself gives
+0.0265 and adding the block drops that to +0.0115. On CYP2C9 the reverse: the block works and
pooling does not. On CYP3A4 neither does anything, which fits — it has the most labels and the
best fit already.

That is an argument for per-enzyme feature sets, and it is also exactly the extra degree of
freedom `src/ablmech.py` pre-registered a warning about: picking the best of four configurations
per enzyme on the same data. On one seed only the main effects on three enzymes clear the noise
floor; the interactions do not. Four seeds are running before any of this goes near the pipeline.

**86. The rule for choosing which compounds get a curve, recovered — and it is not the one
proposed.** `verify/k14_design.py`. Curves exist for about a third of the 4376 screened molecules
and the dependence on the screening reading looks different on every enzyme: a hard threshold at
the active end on CYP2D6, softer on CYP1A2, flat on CYP2C9, and rising toward the WEAK end on
CYP3A4. The decile table reproduces an outside reading's exactly.

The hypothesis offered was that all four are one rule — run a curve where the single point is most
informative about π, which is the Fisher weight E·h·ln10·x/(1+x)², maximal at π = pC₀ — and that
the apparent reversal on CYP3A4 is an artefact of looking along the wrong axis. It predicts that
P(curve) plotted against Fisher weight collapses the four curves into one.

**It does not collapse.** Between-enzyme spread is 0.212 at matched rank and 0.201 at matched
weight, a ratio of 0.95. And the mechanism fails directly on half the enzymes: the median weight
of curved against uncurved compounds is 0.326 / 0.564 on CYP1A2 and 0.466 / 0.529 on CYP3A4 —
**backwards** — against 0.688 / 0.556 and 0.622 / 0.000 on CYP2C9 and CYP2D6.

Significance was tried next and also fails to collapse (spread 0.226), but it localises the
anomaly to one enzyme. The share of curves among screen-significant and non-significant compounds
is 46.5% / 2.4% on CYP1A2, 34.4% / 8.1% on CYP2C9, 43.9% / 0.0% on CYP2D6 — three enzymes where
**a non-hit essentially never gets a curve** — and 37.3% / **66.4%** on CYP3A4.

What CYP3A4 follows instead is the compound's profile on the *other* three enzymes:

    significant elsewhere      n     P(curve on 3A4)
    on none                  115               0.930
    on one                   693               0.632
    on two                  1496               0.429
    on all three            2072               0.298

and it holds within both levels of the compound's own CYP3A4 significance — 0.854 / 0.563 / 0.344
when 3A4 is significant, 0.985 / 0.722 / 0.489 when it is not. So on CYP3A4 curves were run on the
compounds that look **clean elsewhere**, most reliably on those clean everywhere.

That reads as ordinary practice rather than an optimal design: a compound that hits the other
three is already disqualified as a promiscuous inhibitor and its exact CYP3A4 IC₅₀ changes nothing,
while an apparently clean candidate has to be confirmed because a false negative there is
expensive. It also explains the "reversed selection" without inverting any axis — the uncurved
CYP3A4 compounds are the promiscuous strong inhibitors, which is why their median log2fc is
−1.932 against −1.437 for the curved ones.

Consequence for the pseudo-label proposal, which is why this was measured. The unlabelled pool is
not a random remainder on any enzyme, and its composition differs by enzyme in a way now known:
confirmed inactives on CYP1A2, CYP2C9 and CYP2D6, and promiscuous multi-enzyme inhibitors on
CYP3A4. Any use of those readings inherits that composition.

**87. The shift estimate is not invariant to the model, and the amount it moves is as large as
everything we had called uncertainty.** `verify/k15_pooldelta.py`. Item 84 switched the submission
to a pooled model, and the δ rule was re-chosen on the pooled predictions and found unchanged.
That check used the wrong input: `src/shrinkchoice.py` takes the posterior over δ as a fixed file,
and those draws came from inverting the kernel of the **per-enzyme** model. Re-estimating δ with
`src/covshift.py` on pooled predictions on both sides — the blind predictions recomputed with the
pooled learner rather than read from the file built by the other one — gives a different answer.

    enzyme     delta (per-enzyme)              delta (pooled)          amplification
    CYP1A2   +0.045 [-0.361, +0.435]   +0.344 [-0.004, +0.724]        3.4 -> 3.5
    CYP2C9   +0.362 [+0.065, +0.665]   +0.801 [+0.480, +1.160]        2.4 -> 2.6
    CYP2D6   -0.917 [-1.486, -0.434]   -0.405 [-1.109, +0.059]        4.5 -> 4.8
    CYP3A4   +0.740 [+0.485, +0.986]   +0.847 [+0.589, +1.088]        1.7 -> 1.8

The left column reproduces the published estimates exactly, so the harness is sound. A prediction
recorded before the run — that the pooled model's better rank would steepen the kernel and shrink
the amplification factors, most on CYP2D6 — **fails**: they rise slightly on all four. The rank
gain does not reach the inversion.

What moves instead is δ itself, by **+0.30 / +0.44 / +0.51 / +0.11**, all upward. The mechanism is
visible: the gap between the mean test prediction and the mean out-of-fold training prediction
grows under pooling on every enzyme, from +0.014 / +0.147 / −0.190 / +0.437 to +0.101 / +0.300 /
−0.078 / +0.477. The pooled model raises what it says about the test more than what it says about
the training set, and the inversion reads that as a larger label shift.

δ is meant to be a property of the **test's labels**, which do not depend on our model. So this is
a measurement of how badly the kernel-invariance assumption holds — the step section 11 flags as
unverifiable and which had never been quantified. Set against the bootstrap half-widths of 0.40 /
0.30 / 0.53 / 0.25, the model-induced movement is comparable everywhere and **larger than the
whole sampling interval on CYP2C9**, from two models of the same family.

The rule does change, contrary to what item 84 recorded. Under the pooled posterior the mean
criterion picks +0.3 / +0.8 / −0.4 / +0.8 against +0.0 / +0.3 / −0.5 / +0.7, and all four enzymes
pass the "not worse than doing nothing" check rather than three.

**What was adopted.** Not either estimate — neither can be preferred without test labels, and the
larger correction resting on the possibly-biased one is the more expensive mistake. The draws from
both are pooled with equal weight, which is a decision and not a derivation, and the rule is chosen
under that widened posterior: **0 / +0.5 / −0.4 / +0.8**, with CYP1A2 skipped again because its
shift is worth 0.003 of macro against a 36% chance of harm. Per-enzyme fitting is then worth
+0.066, between the +0.047 and +0.101 the two single-model posteriors imply.

This is the first time model uncertainty has entered the δ posterior at all. It makes the method
stricter rather than softer: the systematic term was previously not counted anywhere.

**88. With three models the systematic spread in δ grows rather than settles, and on two enzymes
it exceeds everything we publish as uncertainty.** `verify/k16_modelspread.py`. Item 87 measured
the spread from two models and called it a lower bound. A third — the same learner with absolute
error, structurally the most distant of the three since it estimates a conditional median rather
than a mean — makes the bound larger.

    enzyme    per-enzyme L2   pooled   per-enzyme L1    range   bootstrap half-width   ratio
    CYP1A2           +0.045   +0.344          -0.011    0.354                  0.399    0.89
    CYP2C9           +0.362   +0.801          +0.283    0.518                  0.308    1.68
    CYP2D6           -0.917   -0.405          -1.167    0.762                  0.509    1.50
    CYP3A4           +0.740   +0.847          +0.616    0.230                  0.255    0.90

The range grew on every enzyme against the two-model figures of 0.299 / 0.439 / 0.512 / 0.107, so
the quantity has not converged and three models remain a small sample of one learner family on
one feature set.

Two consequences worth stating separately from the number.

**On CYP2D6 the sign holds across all three and the magnitude does not.** All are negative, which
is the claim the document actually makes, but they run from −0.405 to −1.167. The published
−0.917 is one point in that interval rather than an estimate bracketed by ±0.5, and the honest
form of the CYP2D6 statement is the sign plus a range about twice as wide as reported.

**On CYP1A2 not even the sign survives**: +0.045, +0.344, −0.011. The rule already skips CYP1A2,
adopted in item 58 because the chosen shift was worse than doing nothing in half the draws. That
decision now has a second and independent reason that was not visible then — the direction of the
shift there is not determined by the data at all.

Separately, the seed-to-seed noise in the same estimate is sd 0.005 / 0.005 / 0.048 / 0.012, so
the model term is nine to eighty-eight times the fold-fitting term. The disagreement is
systematic, not sampling, which is why averaging the models is not a way to resolve it: it
produces a third model whose own bias is unmeasured, with a narrow interval that is not earned.

**89. Averaging the models is the right thing to do to the predictions and the wrong thing to do
to the shift estimate, and the two answers have the same cause.** The question was whether to
average posteriors, as item 87 did, or to average the models and estimate δ once on the result.

*As a model, the ensemble wins clearly.* The average of the per-enzyme and the pooled predictions
gives 0.6921 macro after the affine pair against 0.7155 for the per-enzyme model and 0.7089 for
the pooled one — **−0.0234** and **−0.0168** respectively, p = 0.0001 and p = 0.00004 on four
seeds, with rank up by +0.0291 and +0.0155. It clears item 77's ceiling for the same reason
everything that clears it does: it moves rank. The mechanism is ordinary — the two bases see
different training tables, their errors are partly independent, and averaging removes variance
before the shrinkage gets to it — but by the rank criterion it was not guaranteed, and five other
interventions failed the same test. `src/submit.py` now has three modes and submits this one.

*As a way to estimate δ, it is not.* Seed-to-seed noise in the estimate — folds change, test
predictions do not — has sd 0.005 / 0.005 / 0.048 / 0.012, while the model-to-model spread is
0.11 to 0.51. The disagreement is systematic by a factor of nine to eighty-eight, so collapsing
it into a third model buys a narrow interval around an object whose own bias is unmeasured.

The measurement settles it more sharply than the argument does. A prediction was recorded in
`k17_ensdelta.py` before the run: the ensemble's δ need not lie between its components', because
the inversion is nonlinear and an average of predictions is not an average of solutions.

    enzyme    per-enzyme   pooled   ensemble   between?
    CYP1A2        +0.044   +0.335     +0.192   yes
    CYP2C9        +0.351   +0.782     +0.577   yes
    CYP2D6        -0.507   -0.408     -0.686   **no, outside both**

On CYP2D6 — the enzyme the document is built around — the submitted model's own estimate is more
negative than either component. Averaging the models did not produce a compromise there; it
produced the most extreme of the four estimates now available.

*The rule.* Draws from all three are mixed with equal weight, which remains a decision rather
than a derivation, and the rule is chosen on the ensemble's out-of-fold predictions under that
posterior: **0 / +0.5 / −0.5 / +0.8**, CYP1A2 skipped. Per-enzyme fitting is worth +0.066 under
it.

A second pre-registration failed and is worth recording. The widened posterior was expected to
make the rule more cautious. It did the opposite on CYP2D6, −0.4 back to −0.5, because the
ensemble's outlying estimate moved the mixture's centre more than it widened its spread. The
prediction was wrong about which of the two effects would dominate, not about the mechanism.

**90. Morgan is the worst space in which to look for a neighbour, and it is the one used
everywhere.** `verify/k18_nbspace.py`. Two queued methods — neighbour labels as features, and a
Gaussian process with a similarity kernel — both need a definition of "close", and this
repository has always used Tanimoto over Morgan counts without testing it. The test needs no
model: give each compound the label of its nearest neighbour among the other folds and
rank-correlate with the truth.

    space                  1A2      2C9      2D6      3A4     mean
    Morgan / Tanimoto   0.1887  -0.0512   0.2126   0.0638   0.1034
    pharmacophore pairs 0.1965   0.1407   0.1538   0.2575   0.1871
    chemprop embedding  0.2109   0.1357   0.2273   0.3148   0.2222
    RDKit descriptors   0.2110   0.2774   0.1268   0.4289   0.2610

**Morgan is last of the four and carries nothing at all on CYP2C9** (−0.05). Descriptors win by
2.5× over it, and the best space differs by enzyme: the embedding on CYP2D6, descriptors on
CYP2C9 and CYP3A4. An outside prediction that pharmacophore similarity would beat substructure
similarity holds — 0.187 against 0.103 — but its consequence, that neighbours should be built in
pharmacophore space, does not: two other spaces beat it.

A one-nearest-neighbour rule is a deliberately weak predictor and these numbers are far below the
boosting's. The comparison is between spaces.

**91. The split protects in both spaces, twice as strongly in the one it clusters on.** The
finding above raises a question about the Butina split itself: it groups by Morgan similarity,
and activity information turns out to travel by descriptors. Measured against a random split of
the same fold sizes, on cross-fold nearest-neighbour similarity:

    quantile of the random split's own distribution     random   Butina   ratio
    Morgan, top 1%                                        1.0%     0.2%    4.0
    Morgan, top 5%                                        4.9%     4.0%    1.2
    descriptors, top 1%                                   1.0%     0.5%    2.1
    descriptors, top 5%                                   5.0%     4.1%    1.2

Medians are identical to three decimals in both spaces, 0.431 and 0.731, which reproduces item
22's finding that the cluster split earns its keep in the tail rather than the median — and
extends it to descriptors. In the extreme tail the split cuts Morgan neighbours by 4× and
descriptor neighbours by 2.1×.

So the split is not blind to the space that matters, but it protects there half as strongly.
Held-out estimates are therefore somewhat optimistic in descriptor terms. The size of the gap is
a factor of two in tail protection, not an order of magnitude, so this qualifies the split rather
than invalidating it — and the qualification is that our numbers describe generalisation to new
substructures better than generalisation to new property combinations.

A first attempt at this comparison was confounded and is worth recording as a caution: comparing
within-fold to cross-fold nearest-neighbour similarity finds the cross-fold neighbour *closer*,
because the training folds hold four times as many compounds. The control has to be a random
split of the same sizes, not the other half of the same split.

**92. A Gaussian process is worse than the boosting alone and improves it in the ensemble.**
`src/gp.py`. The proposal was a GP with a Tanimoto kernel over Morgan fingerprints, the
cheminformatics standard. Item 90 having found Morgan to be the least informative of four
neighbour spaces, the kernel goes on standardised RDKit descriptors instead, clipped at five
deviations because `Ipc` alone spans fourteen orders of magnitude and would otherwise define the
metric by itself. Exact inference — 1285 to 2335 rows per enzyme is one to three seconds of
Cholesky — with the lengthscale and noise chosen by marginal likelihood on the training folds.

    four seeds                pair     rank
    per-enzyme boosting     0.7155   0.5630
    pooled boosting         0.7089   0.5767
    GP alone                0.7237   0.5535
    boosting ensemble       0.6921   0.5921
    all three, equal        0.6845   0.5994

**The GP is the worst of the three on its own and the ensemble is better with it**: −0.0076
against the two-model ensemble, p = 0.0009, with rank up 0.0073, p = 0.0003. From the per-enzyme
baseline the three-way ensemble is worth −0.0310. Weighting hardly matters — equal gives −0.0076,
three-to-one −0.0075, five-to-one −0.0062 — so equal weight is used, being simplest and no worse.

This is the fourth intervention to clear item 77's ceiling and the first from a different model
family. The reason it works is the reason it was proposed: a GP posterior mean is a
similarity-weighted average of neighbouring labels, which is the local structure a tree cannot
represent, so its errors are decorrelated from the boosting's in a way the pooled and per-enzyme
models' are not. Item 22's measurement that the test sits closer to the training set (0.587) than
the training set sits to itself (0.435) is the setting that favours it.

The first run had to be discarded and the reason is worth recording: the lengthscale sat on the
grid's lower edge and the noise on its upper edge for three enzymes of four, so the grid chose
rather than the marginal likelihood. That is item 59's error in a new place. The grid now spans
0.0625 to 4 and 0.01 to 3, an edge check prints, and nothing sits on a boundary.

`src/submit.py` gains the GP as a third ensemble member, with `ансамбль-без-GP` preserved to
reproduce the previous state.

**93. The heteroscedastic decision layer fails its own pre-registered threshold.** `src/abhetero.py`.
The argument for it was a proof rather than a hope: the affine pair is fitted once per fold and is
therefore global, so a correction indexed by the individual compound's uncertainty cannot be
expressed by it. It was the only proposal on the table that escapes item 77's ceiling by
construction. Two conditioners, both fitted out of fold with nested cross-validation so that bin
edges and per-bin pairs never see the rows they are applied to.

    seed 0        raw    global pair   by band width   by model spread
    CYP1A2     0.8786         0.8101          0.8105            0.8142
    CYP2C9     0.6908         0.6546          0.6604            0.6583
    CYP2D6     0.9803         0.9052          0.8904            0.8996
    CYP3A4     0.5194         0.4853          0.4834            0.4838
    macro      0.7673         0.7138          0.7112            0.7140

The threshold recorded before the run was a macro gain of 0.005 for either conditioner, or a
favourable sign on three enzymes of four. Band width gives −0.0026 and two of four; model spread
gives +0.0002 and two of four. Neither passes, so the four-seed run is not made.

What is not zero and should not be buried: on **CYP2D6 the band-width conditioner is worth
−0.0148**, and CYP2D6 is the enzyme with the widest variation in band width and the worst fit. The
effect exists, is concentrated on one enzyme, and is cancelled by a loss of the same size on
CYP2C9. That is an open thread rather than a clean null.

Item 80's criterion predicted this and the prediction was recorded only after the run was
launched, which is worth admitting: a per-bin affine map is monotone inside each bin and moves
rank only between bins, so by the rank rule it should not have been expected to survive. Being
inexpressible by the global pair and being useful are different properties, and only the second
one matters.

**94. The mechanistic block should not be trimmed to its geometric part.** Item 82 found that on
CYP2D6 six geometric features beat all thirty by +0.0127 of rank, which raised the question of
whether the block should be cut down. Running the same split on CYP1A2 and CYP3A4 completes the
macro, and the answer is no.

On those two enzymes every subset of the block sits within 0.002 of every other and of no block at
all, which agrees with item 81 — the block does nothing there in any composition. Combining with
the CYP2D6 and CYP2C9 halves:

    macro after the affine pair, four seeds     whole block (30)   geometry only (6)
                                                          0.7156              0.7177

Geometry wins on CYP2D6 and loses on CYP2C9, and the macro comes out **0.0021 in favour of keeping
the whole block**. That is scenario two of the three `src/ablmech.py` pre-registered: a gain on one
enzyme, a wash on the macro, and nothing left but per-enzyme feature sets — which is the extra
degree of freedom the same file warned against choosing on this data.

**95. Choosing the shift by worst case across models agrees with the mixture centre.** An outside
reading suggested that after item 89 — where the ensemble's estimate landed outside its
components' range — the rule should be picked by the worst case over models rather than the centre
of their mixture, so that one nonlinear model cannot drag the answer. Compared on the three
posteriors available:

    enzyme    model medians           mixture centre   minimax
    CYP1A2    +0.04 +0.33 +0.19               +0.195      +0.2
    CYP2C9    +0.35 +0.78 +0.57               +0.574      +0.6
    CYP2D6    -0.51 -0.41 -0.69               -0.538      -0.5
    CYP3A4    +0.74 +0.84 +0.78               +0.781      +0.8

They agree within one grid step everywhere. **That conclusion is wrong**, and the caveat recorded
beside it is why: minimising the largest distance to a model's median is not minimising the
largest expected loss under a model's posterior, and the loss is asymmetric in δ. The proper
version was then computed — `src/shrinkchoice.py` gained a `--dump-ac` option for the
assumed-by-true matrices — and the two criteria differ on three enzymes of four.

    enzyme    minimax   mixture centre   pooled     GP   ensemble3
    CYP1A2       +0.3             +0.2     +0.3   +0.2        +0.2
    CYP2C9       +0.8             +0.6     +0.8   +0.5        +0.6
    CYP2D6       -0.4             -0.7     -0.4   -1.0        -0.8
    CYP3A4       +0.9             +0.9     +0.8   +0.9        +0.8

The gap is largest on CYP2D6, 0.3, which is where the posteriors disagree most — the pooled model
says −0.4 and the GP says −1.0. Measured in worst-case expected ST-RAE the mixture centre costs
**+0.0108** of macro against the minimax, and 0.0266 of that is CYP2D6 alone.

*Which to adopt is a decision and not a computation.* The same choice was made once already inside
a single posterior — item 57, worst case against mean, worth 0.05 — and the mean was adopted,
because the goal is the best expected score rather than insurance. Taking the minimax across
models while keeping the mean within one is inconsistent. Against that: inside a posterior the
uncertainty is sampling and averaging is natural, while across models it is structural, and there
is no reason the truth should be the average of what models believe. The mixture centre stays as
the default and the alternative is recorded with its price, because by item 57's own finding this
kind of choice moves more than the estimates do.

**96. Fingerprints are not ballast, and the reason completes item 90.** Item 90 found Morgan the
least informative of four neighbour spaces, carrying nothing on CYP2C9. The natural inference is
that 2048 of 2295 columns sit in the space where signal is weakest and could be dropped. The
answer was already in `oof.json` and it is no.

    set              columns      raw   +pair     rank
    FP                  2048   0.8333  0.7566   0.5249
    DESC                 217   0.8071  0.7520   0.5166
    MECH                  30   1.0350  0.9013   0.3054
    FP+DESC             2265   0.7729  0.7215   0.5487
    DESC+MECH            247   0.7996  0.7420   0.5330
    FP+DESC+MECH        2295   0.7673  0.7150   0.5650

Dropping the fingerprints costs **0.027 after the pair and 0.032 of rank**, more than almost
anything this log records as a gain.

The reconciliation is the point and belongs in the document. **Substructure is a poor metric and a
good feature.** Tanimoto over the whole Morgan vector is a weak measure of closeness; individual
Morgan bits are strong variables to split on. A tree asks "is this fragment present", not "how
similar are these two vectors", so the two facts are about different objects and neither implies
the other. Both are now measured, which is unusual.

**97. The GP's contribution shrinks as neighbours get closer, which is the opposite of why it was
added.** An outside reading argued that the GP is evaluated in a harsher regime than it will work
in — held-out compounds sit at 0.435 median similarity to the training set, the test at 0.587 —
so its ensemble weight should be raised above equal. Stratifying the held-out compounds by their
cross-fold nearest-neighbour similarity and measuring the GP's contribution in each stratum:

    stratum       n     ensemble of 2   ensemble of 3   GP contribution
    < 0.35     1086            0.7701          0.7540           -0.0161
    0.35-0.45  2553            0.7051          0.6967           -0.0084
    0.45-0.55  2090            0.6746          0.6699           -0.0047
    > 0.55      795            0.7016          0.7012           -0.0004

The contribution falls monotonically with proximity and is zero in the closest stratum — which is
where the test sits. So the extrapolation runs the other way: −0.0076 measured out of fold
**overstates** what the GP will deliver on the test, and its weight should if anything be lowered.

That also revises why item 92 works. The GP was justified as exploiting close neighbours; it
earns its place where neighbours are **far**, and there because everything is wrong there and a
differently-wrong model helps most. Its value is generic variance reduction in hard regions, not
similarity exploitation.

A second thing falls out. Accuracy is **not monotone** in proximity: 0.7701 / 0.7051 / 0.6746 /
0.7016. The closest stratum is worse than the middle one, which is the signature of activity
cliffs — near-identical structures with sharply different activity. That is a measured population
on which similarity actively misleads, and it is worth its own look.

The rightmost stratum holds 795 observations across four seeds, so its −0.0004 is noisy on its
own; the monotone trend across four strata is the evidence, not any single one.

**98. The diversity well is nearly dry for this family of candidates.** Item 92's success invites
a general rule — accept an ensemble member by how decorrelated its errors are, not by how accurate
it is — and that is measurable before including anything. Out-of-fold residual correlations across
six candidates, averaged over the four enzymes:

    mean correlation with the others, lower is more valuable
      GP          0.908
      pooled      0.909
      DESC+MECH   0.930
      FP+DESC     0.932
      L1          0.932
      per-enzyme  0.936

Every pair sits between 0.89 and 0.97. The GP is the least correlated, which agrees with it being
the member that helped, but its margin over the pooled model is 0.001 — the criterion ranks them
correctly and separates them barely. The most similar pair is the per-enzyme model and FP+DESC at
0.969, which is as expected since they differ only by the mechanistic block.

The reading is sobering rather than encouraging: within gradient boosting on these features the
candidates are near-duplicates, and further ensembling of the same kind has little left to
extract. Whether a structurally different member breaks out of the band is a prediction to check
against the ridge and kNN arms now running.

**99. The shrinkage grid was parameterised so that it must fail as the model improves.** Building
the submission on the three-member ensemble printed a grid-edge warning on CYP2C9: offset +3.00,
which is the grid's maximum, at lambda 0.90. A parameter sitting on a boundary was chosen by the
grid rather than by the data — the third occurrence of item 59's error, and this one inside the
submitted pipeline.

The cause is structural rather than a matter of range. The shift the predictions actually receive
is `(1 - lambda) * offset`, so for a fixed shift the offset needed grows without bound as lambda
approaches one. Better models are trusted more, lambda rises — 0.58 to 0.90 on CYP2C9 between the
per-enzyme model and the ensemble — and the same correction demands an ever larger offset. **A grid
over the offset is guaranteed to break as the model gets better.**

The fix is a reparameterisation, not a wider grid. The family is unchanged:

    c + lambda*(p - c)   with c = mu + offset   ==   lambda*p + (1-lambda)*mu + s,   s = (1-lambda)*offset

so gridding over the shift `s` covers the same set of transformations. But `s` is in pIC50 units
and bounded by anything one would believe about a shift, while the offset is bounded by nothing.
Verified identical to 1e-15 on synthetic data at three (offset, lambda) pairs, including the
edge case that triggered this.

One error inside the fix, caught before it ran. Applying the new form used the mean of the TEST
predictions where the fit had used the mean of the out-of-fold training predictions. Those differ
by exactly the marginal shift between the two sets — which is the quantity `delta` estimates
separately — so the substitution would have silently applied that shift twice.

`src/submit.py` now grids over the shift; `src/shrinkchoice.py` is untouched, since its published
numbers were produced with lambdas below 0.8 where the offset grid did not bind.

**100. Three more members tried, one accepted, and the decorrelation criterion picks the wrong
one.** `src/ablweak.py`. Item 92's regime — worse alone, better in the ensemble — suggests adding
members for diversity rather than accuracy, and an outside reading proposed choosing them by
residual decorrelation measured before inclusion. Three candidates, four seeds:

    alone                pair     rank
    random forest      0.7256   0.5602
    ridge on desc      0.7239   0.5599
    kNN on desc        0.7952   0.4816
    (per-enzyme GBM)   0.7156   0.5630

    added to the three-member ensemble    pair     rank   change        p
    current three                       0.6845   0.5994        -        -
    + ridge                             0.6819   0.6009   -0.0026   0.0005
    + forest                            0.6866   0.5981   +0.0021   0.0033
    + kNN                               0.6887   0.5959   +0.0042   0.0013
    all six                             0.6873   0.5968   +0.0028   0.0006

**Only the ridge helps.** A linear function of 247 descriptors, six seconds to fit, improves an
ensemble of two boostings and a Gaussian process — because it errs smoothly where all three err
locally. The forest and kNN both make it worse, and adding everything is worse than adding
nothing.

The criterion fails, and instructively. Ranked by residual decorrelation the order is kNN 0.892,
pooled 0.896, per-enzyme 0.899, ridge 0.908, GP 0.918, forest 0.930 — so **the most decorrelated
candidate is the most harmful one**. Decorrelation is necessary and not sufficient: kNN's rank of
0.4816 against 0.56 for the others makes it too weak for its independence to be worth anything.
The proposal did carry the qualifier "subject to individual accuracy above a floor", and this
measures how much work that qualifier does — without it the criterion inverts.

`src/submit.py` takes the ridge as a fourth member. Membership is now settled by measurement in
every case rather than by the principle that more diversity is better, which this item shows is
false as stated.

**101. FCFP loses as a replacement and does nothing as an addition.**
*(This item was first written from seed 0 alone and claimed a win for the concatenation. The
completed four-seed run refutes it; the corrected numbers are below and item 119 records how the
error happened.)*
`src/ablfcfp.py`. An outside reading proposed the pharmacophoric fingerprint on the ground that
CYP recognition is about donors, acceptors and aromatic character rather than about the exact
atoms. Three arms, four seeds:

    набор              пара     ранг      1A2      2C9      2D6      3A4
    ECFP (как сейчас) 0.7156   0.5630   0.4957   0.5917   0.4035   0.7613
    FCFP вместо ECFP  0.7275   0.5497   0.4755   0.5786   0.3898   0.7549
    оба               0.7156   0.5627   0.4890   0.5984   0.3985   0.7649

As a replacement FCFP is clearly worse — it loses 0.0133 of rank, consistently on every seed.
**Concatenated it is worth nothing**: the pair moves by +0.000025 and the rank by -0.00035, three
orders of magnitude under the floor.

The per-enzyme pattern is real even though the total is not. The concatenation helps CYP2C9
(+0.0067) and CYP3A4 (+0.0036) on every seed and hurts CYP1A2 (-0.0067) and CYP2D6 (-0.0051) on
every seed — and CYP2C9 and CYP3A4 are the pair that correlate at 0.688 (item 111). Whatever the
pharmacophoric view adds, it adds it to one half of the panel and takes it from the other, and the
two cancel.

An earlier draft of this item read the seed-0 numbers as a win and generalised from them, citing
item 68's learned encoder as a precedent for "worse as a replacement, better as an addition". Both
halves have since gone. Item 117 withdrew the encoder precedent — its concatenated gain was raw
ST-RAE and does not survive the affine pair — and the three remaining seeds withdrew this one.
**The rule that a representation should always be tried concatenated is still worth keeping**,
because the substitution arm genuinely measures the wrong thing; what is not supported is any
expectation that the concatenated arm will win.

**102. The aggregator suspects are chemistry, not a solubility artefact — and the flag proves it.**
`src/ablagg.py`. The concern was that highly lipophilic compounds precipitate at assay
concentration, so their curves are wrong and their rows are noise. The falsification proposed with
it was that a real artefact would improve all four enzymes at once. Four arms, seed 0:

    рука                 пара     ранг    1A2     2C9     2D6     3A4
    база               0.7150   0.5651  0.496   0.597   0.403   0.765
    флаг               0.7126   0.5697  0.504   0.596   0.415   0.763
    вес по структуре   0.7128   0.5695  0.503   0.595   0.418   0.763
    вес со скринингом  0.7126   0.5692  0.500   0.591   0.417   0.768

All three interventions help by the same amount, about 0.0045 of rank, and the gain sits on CYP1A2
and CYP2D6 while CYP2C9 and CYP3A4 do not move. The pre-registered test therefore fails: two
enzymes out of four, not four.

The coincidence of the three arms says more than their size. **If these rows were corrupted
labels, down-weighting would have to beat flagging** — a wrong label cannot be explained away by
handing the model an extra column, it can only be discounted. The two are equal to within 0.0002,
so what the suspect column carries is signal the model can use, not noise to remove. The flag
stays in as a feature; the down-weighting does not go in, because it throws away rows to buy the
same thing.

**103. The reparameterised grid does not reproduce the old submission bit for bit, and should not.**
Item 99 replaced the offset grid with a grid over the prediction shift. Re-running the submission:

    фермент   lambda    сдвиг   эталон   среднее   эталон
    CYP1A2      0.58   +0.130   +0.126     5.091    5.087
    CYP2C9      0.76   +0.240   +0.247     4.929    4.933
    CYP2D6      0.52   -0.190   -0.192     4.478    4.476
    CYP3A4      0.74   +0.430   +0.429     4.852    4.851

Every shift differs, by 0.001 to 0.007. This is resolution, not behaviour: the old grid stepped
the offset by 0.05, which at lambda 0.58 is a step of 0.05 x (1 - 0.58) = 0.021 in the shift, while
the new grid steps the shift by 0.01 directly. Each optimum has moved by less than one old step,
and the new grid is the finer of the two. It is recorded here because "the regression passes"
would have been the wrong words and the four decimals in the document come from these numbers.

**104. The geometry prediction failed, and the two halves of the block are not separable.**
`src/ablshape.py`, `--sets`, four seeds. The block was justified by a per-enzyme prediction:
planarity should help CYP1A2, volume and flexibility should help CYP3A4, and CYP2D6 should not
move. Measured:

    набор           пара     ранг      1A2      2C9      2D6      3A4
    без формы     0.7156   0.5630   0.4957   0.5917   0.4035   0.7613
    с формой      0.7132   0.5664   0.5017   0.5914   0.4122   0.7605

The block does help — 0.0034 of rank, above the 0.007 floor only marginally but consistent across
four seeds — and it helps on **CYP2D6 (+0.0087) and CYP1A2 (+0.0060), with CYP3A4 flat**. That is
the enzyme the prediction said would not move, and not the enzyme it said would. The block is kept
on measurement and its stated mechanism is withdrawn.

Splitting it does not rescue the story either:

    только индексы формы   0.7152   0.5649   0.5022   0.5904   0.4066   0.7603
    только расстояния      0.7148   0.5650   0.4989   0.5922   0.4076   0.7615

The two halves are indistinguishable at the macro level and each recovers about half the CYP2D6
gain, while together they deliver more than the sum. Item 91's attribution of the effect to the
3D geometric distances alone was drawn from a partial run and does not hold: neither half is the
driver.

**105. There is no censoring spike, and the real finding is larger than the one that was looked for.**
The proposal was that labels pile up at the instrument's floor, pC0 = -log10(49.5 uM) = 4.305, and
that such rows are left-censored and should enter as inequalities. The histogram says no: the
density passes through 4.305 smoothly in all four enzymes, the most frequent single value is shared
by two compounds, and 20.5 / 31.9 / 15.3 / 51.4 per cent of the labels lie *below* the floor. The
organisers did not clamp; they let the fit extrapolate past the top tested concentration.

Asking where resolution actually ends, rather than assuming it, gives a much sharper picture. An
unresolved curve does not announce itself with a spike — the fit still returns a number — it
announces itself by the confidence band widening. Median band width against potency, CYP2D6 as the
example:

    pIC50           n   медиана ширины   доля шире 1.0
    [0.00,3.50)   102        2.609          100.0 %
    [3.50,4.00)    27        2.494          100.0 %
    [4.00,4.30)   100        0.513           18.0 %
    [4.30,4.60)   378        0.263            0.0 %
    [4.60,5.00)   340        0.297            0.3 %

The same cliff appears in all four. Emax is flat at about -1.0 throughout, so these compounds do
reach full suppression; it is the inflection that is unpinned, not the depth. Across the four
enzymes **19.3 % of the 6525 labels carry a band wider than 1.0 pIC50** and 12 % wider than 2.0.

What that costs, measured on the saved out-of-fold predictions:

    фермент   доля строк   доля кв.потери   доля ST-RAE
    CYP1A2       14.6 %         44.7 %         20.0 %
    CYP2C9       17.0 %         33.8 %         12.6 %
    CYP2D6        9.9 %         41.2 %         12.8 %
    CYP3A4       29.4 %         44.5 %         11.3 %

**A squared loss spends about forty per cent of its effort on rows that supply about thirteen per
cent of the score.** On CYP3A4, 73.9 % of the predictions on those rows already land inside the
band — the error is already free and the objective is still pushing on it. On CYP2C9 the median
such row has an absolute residual of 0.605 and an ST-RAE of exactly 0.000. The misallocation
compounds, because these are also the weakest compounds and therefore the ones a squared loss
weights most heavily.

The censoring hypothesis is closed. What replaces it is the measured case for item 106's objective.

**106. Averaging ranks instead of values changes nothing, and the reason it fails is worth keeping.**
`verify/k21_borda.py`. Since only rank survives the affine pair, an outside reading proposed that
the ensemble average ranks rather than values — Borda — so that no member's scale drags the joint
answer around. Four seeds, ranks mapped back through the empirical quantile function of the labels:

    членов   значения (пара/ранг)   ранги (пара/ранг)   разница пары   разница ранга
        2      0.6921 / 0.5921      0.6964 / 0.5922        +0.0043        +0.0000
        3      0.6845 / 0.5994      0.6923 / 0.5991        +0.0078        -0.0003
        4      0.6819 / 0.6009      0.6921 / 0.6000        +0.0102        -0.0009

**The rank is identical** — every difference is an order of magnitude below the 0.007 floor. The
members' scales are commensurate and there was nothing to fix. The pair, meanwhile, degrades, and
degrades further with each member added. The mechanism is clean: the quantile map forces the
predictions to have the spread of the labels, which is exactly the over-dispersion that fitting
lambda below one exists to remove, and a two-parameter affine pair cannot undo a quantile
transform. Half an hour, question closed, and it doubles as a check that the shrinkage is doing
real work rather than compensating for something else.

**107. kNN does not recover where the neighbours are near — it gets worse, and that is not what
the rejection assumed.** `verify/k20_strat.py`. Item 97 found the cross-validation runs at a
nearest-neighbour similarity of 0.435 while the test set sits at 0.587, so every choice in this log
was made in a harder regime than the one it will be scored in. The natural first suspect was item
100's rejection of kNN: a nearest-neighbour method is exactly the method whose accuracy should
depend on how near the neighbours are, and it was judged at 0.435. Rank inside each layer, four
seeds, four enzymes:

    модель                  < 0.35   0.35-0.45   0.45-0.55      > 0.55
    n на слой                 1086        2553        2090         795
    бустинг поферментно     0.5280      0.5350      0.5872      0.5657
    бустинг пул             0.5143      0.5570      0.6061      0.6023
    GP                      0.5328      0.5317      0.5660      0.5315
    лес                     0.5278      0.5363      0.5827      0.5479
    гребневая               0.5380      0.5404      0.5719      0.5458
    kNN                     0.4626      0.4599      0.5038      0.4567

Columns are not comparable to one another — the layers have different label spreads — so only
comparisons down a column are drawn. kNN's deficit against the per-enzyme boosting **widens** from
-0.0655 in the far layer to -0.1090 in the near one. The rejection was not delivered in the wrong
regime; if anything it was delivered in the regime most favourable to kNN.

The reason is the interesting part, and it does not close the difference-model proposal that this
check was run to gate. kNN assumes the label of a near neighbour transfers unchanged. Item 96
already measured that similarity actively misleads in the middle layers, and item 99's activity
cliffs are by construction near pairs whose labels differ. So kNN fails hardest in the near layer
precisely because that is where the assumption Delta-y = 0 is worst — which is the one thing a
model of Delta-y is for. The precondition test does not decide the proposal; it relocates it. The
question to ask of the difference model is not whether neighbours are near but whether Delta-y
between near neighbours is predictable at all, and that is a different measurement.

**108. The layer effect is real and monotone; fitting weights on it loses in every layer.**
`verify/k22_layerboot.py`. The table in item 107 is eyeballed from about 200 compounds per enzyme
per seed in the right-hand column, where a Spearman carries a standard error near 0.06. Paired
bootstrap over compounds, 1000 resamples, the same resample indices given to both models:

    сравнение                   < 0.35            0.35-0.45          0.45-0.55            > 0.55
    пул - поферментно   -0.0136 [-.031,+.004] +0.0221 [+.011,+.032] +0.0191 [+.007,+.031] +0.0367 [+.011,+.063]
    GP - поферментно    +0.0044 [-.012,+.022] -0.0034 [-.015,+.007] -0.0210 [-.033,-.009] -0.0336 [-.060,-.007]
    гребневая - пофермент +0.0098 [-.007,+.028] +0.0053 [-.005,+.017] -0.0153 [-.027,-.002] -0.0202 [-.045,+.004]

Both principal effects hold and both are monotone across the four layers. **Pooling's advantage in
the layer the test sits in is +0.0367 of rank, against the +0.0137 of rank it is worth overall** —
2.7 times larger, and in the far layer pooling is actually harmful.

*(An earlier version of this paragraph set the +0.0367 against 0.0067 and called it five times
larger. That 0.0067 is pooling's gain in ST-RAE after the pair, which sits in the column next to
the rank gain in the very table it was read from; the two are different quantities and the ratio
was wrong. The same slip reached section 13 of the document and has been corrected there. The
direction and the monotone trend are unaffected — only the multiple was overstated.)* The Gaussian process runs the other
way and is worst exactly where the test lives. Item 97 was not a one-off.

Acting on it directly fails. Ensemble weights fitted inside each layer, leave-one-fold-out so no
compound contributes to the weights that score it:

    слой         равные веса   подогнанные   разница   средние веса
    < 0.35            0.7367        0.7559    +0.0192   пофермент .35 пул .17 GP .25 гребн .23
    0.35-0.45         0.6925        0.6968    +0.0043   пофермент .22 пул .38 GP .17 гребн .23
    0.45-0.55         0.6640        0.6687    +0.0047   пофермент .20 пул .49 GP .11 гребн .20
    > 0.55            0.7185        0.7345    +0.0160   пофермент .20 пул .62 GP .04 гребн .13

Worse in all four. And the weights are not learning the wrong thing — in the near layer the fit
puts 0.62 on pooling and 0.04 on the Gaussian process, which is precisely the ordering the
bootstrap established. **Four free parameters estimated on about 160 rows per fold cost more in
variance than the misweighting costs in bias.** The effect is real and this way of spending it is
not; equal weights survive on a measurement rather than on inertia.

That leaves the effect unspent rather than refuted, which is what `verify/k23_tilt.py` is for: one
parameter instead of sixteen, chosen by holding out a whole seed.

**109. The difference model on near neighbours has neither a training set nor a starting point.**
Item 107 relocated the proposal: the question is not whether the neighbours are near but whether
the label difference across them is predictable. Two counts decide it before anything is built —
how many co-labelled pairs exist above the threshold, and how the spread of Delta-y across them
compares with the boosting's own residual, since Delta-y = 0 is the naive rule the model would
have to beat.

    фермент       n   порог     пар   пар/мол   std dy   станд. остаток
    CYP1A2     1412    0.55      58      0.04    1.188            0.913
    CYP2C9     1285    0.55      48      0.04    1.176            0.624
    CYP2D6     1493    0.55      62      0.04    0.932            0.846
    CYP3A4     2335    0.55     996      0.43    1.014            0.698

Both preconditions fail at once. There is no training set — 48 to 62 pairs on three of the four
enzymes, against the "an order of magnitude more pairs than molecules" the proposal assumed. And
**the neighbour's label is a worse starting point than the model's own prediction in every
enzyme**: the spread of Delta-y across near pairs exceeds the residual the boosting already
achieves, by a factor of nearly two on CYP2C9. Correcting a neighbour's label starts further from
the truth than simply predicting.

This is the second proposal to die on the same count. Item 88 refuted the matched-molecular-pair
method by finding 4230 distinct transformations across 4265 pairs and none with five examples.
The cause is structural and was visible from the start: Butina at 0.35 gives 4703 clusters for
about 4900 molecules, which is a set with almost no near neighbours by construction. Any method
whose unit of learning is a *pair* is refuted by that number, and the class should be treated as
closed rather than re-proposed one member at a time.

The count did leave something behind. Twelve per cent of rows have some training neighbour above
0.55 while co-labelled pairs that near are essentially absent — so the near neighbours exist and
carry labels for other enzymes. `verify/k24_visible.py` follows that.

**110. Pooling does not win by borrowing neighbours, and its largest effect stays unexplained.**
`verify/k24_visible.py`. Item 108 left the sign flip in item 107's table without a mechanism:
pooling is worth +0.0367 of rank where neighbours are near and -0.0136 where they are far. Item
109 supplied a candidate. The label matrix is sparse, so a molecule's nearest structural neighbour
usually carries a label for a different enzyme; the per-enzyme model cannot see that row at all
and the pooled model can. On that reading pooling's advantage is the neighbours it reaches that
the per-enzyme model cannot, and it should vanish where there are none to reach.

The invisible fraction is large enough for the mechanism to have worked:

    фермент   меток   медиана nn все   медиана nn с меткой   доля с зазором > 0.02
    CYP1A2     1412            0.415                 0.371                  40.2 %
    CYP2C9     1285            0.452                 0.423                  36.1 %
    CYP2D6     1493            0.392                 0.341                  44.5 %
    CYP3A4     2335            0.461                 0.442                  20.2 %

It did not:

    зазор          n   преимущество пула
    ~ 0         4398             +0.0165
    0.02-0.08   1153             +0.0035
    0.08-0.18    756             +0.0125
    > 0.18       216             +0.0450

The prediction was zero in the top row and a monotone rise. What happens is neither: the
advantage is **large precisely where the per-enzyme model can already see the whole
neighbourhood**, then falls, then rises again. The widest-gap bin is consistent with borrowing
and holds 216 rows; the top bin holds 4398 and refutes it.

So pooling is not reaching for neighbours. What is left is transfer of *function* rather than of
neighbours — the four enzymes share enough structure-activity relationship (item 86's promiscuity)
that fitting all 6525 rows at once estimates the shared part better, with the indicator carrying
the differences. That reading is consistent with the numbers but it was not tested here, and it
does not explain the sign flip by itself.

**The largest single effect available to the submission is therefore unexplained, and that is
worse than it sounds.** An effect we cannot attribute is an effect we cannot argue will survive
the move to the test set; +0.0367 measured in a proxy layer is not the same claim as +0.0367 on
the real thing. This is now the most important open question in the log.

**111. Pooling does not transfer shared function either — the relation is perfectly inverted.**
Item 110 left transfer of the structure-activity relationship as the surviving reading. It makes a
prediction that costs nothing to check: the amount of function available to transfer is the rank
correlation between two enzymes' labels on the molecules carrying both, and the enzyme most
correlated with the rest should gain most.

    фермент   меток   средняя корреляция с прочими   выигрыш пула (ранг)
    CYP2D6     1493                        -0.002                +0.0374
    CYP1A2     1412                         0.288                +0.0157
    CYP2C9     1285                         0.330                +0.0071
    CYP3A4     2335                         0.375                -0.0057

Monotone, and inverted in every step: rho = -1.000 across the four. **Pooling helps most exactly
where there is no shared function, and hurts the one enzyme that shares the most.** CYP2C9 and
CYP3A4 correlate at 0.688 and pooling does nothing for either; CYP2D6 correlates with nobody —
-0.120 against CYP2C9 — and takes the largest gain in the whole log.

Both natural mechanisms are now refuted. What the inversion suggests instead is **contrast**: with
the enzyme indicator in the design a tree can learn "this split matters for CYP2D6 and not for
CYP3A4", which is strictly more than "this split matters for CYP2D6", and a per-enzyme model
cannot represent it at all. On that reading a partner is useful in proportion to how much its
behaviour *differs*, and a near-duplicate partner adds rows without adding information.

Four points cannot establish this, and worse, they cannot separate it from plain sample size:
CYP3A4 has both the most labels and the highest correlation, so the two candidate explanations are
themselves confounded across these four numbers.

Two arms separate them and both are stated here before they report.

**`src/ablpool.py --arms "пул слепой"`** — the same pooled table with the enzyme indicator zeroed
rather than removed, so the design matrix keeps its width and the arms differ only in what the
model is allowed to know. Transfer of shared function works without the indicator: a model that
cannot tell the enzymes apart still estimates their common part from 6525 rows instead of 1285.
Contrast cannot work without it at all — "this split matters for CYP2D6 and not for CYP3A4" is
unrepresentable when the two are indistinguishable. **So: if the blind arm still beats the
per-enzyme model, the gain is rows and regularisation and the indicator is incidental. If it
collapses to the per-enzyme model or below, the indicator carries the effect.**

**`src/ablpool.py --arms "пул+TDI"`** — the pre-incubation arm of the same four enzymes, which
more than doubles the table to 13063 rows with maximally correlated ones. Under sample size this
should be the largest gain available anywhere. Under contrast it adds rows and no information.

`src/ablpair.py` would have separated them per partner as well, but at 2048 seconds for the
cheapest of twelve arms it was stopped as too expensive for what it buys. Its one completed cell
is worth recording as a hint and not more: CYP1A2 pooled with CYP2C9 alone reaches rank 0.5136
against 0.4957 on its own — **more than the full four-way pool gives it, 0.5084** — so partners
are not additive and one well-chosen partner can beat all three. One cell, seed 0.

**112. Tilting the ensemble toward pooling: consistent, and under the noise floor.**
`verify/k23_tilt.py`. Item 108 established the layer effect and refuted the sixteen-parameter way
of spending it. One parameter instead: slide from equal weights toward the pooling-heavy corner
that the near-layer fit found, w(t) = (1-t)·equal + t·(0.20, 0.62, 0.04, 0.13), with t chosen by
holding out a whole seed.

    t     близкий > 0.55        всё   далёкий < 0.45
    0.00          0.7058     0.6819           0.7030
    0.30          0.6996     0.6799           0.7032
    0.70          0.6963     0.6814           0.7073
    1.00          0.6973     0.6857           0.7132

Held out, with t chosen on the overall metric rather than on the layer — the more honest of the
two selections, since it never looks at the thing it is scored on:

    сид   t      близкий   разница        всё   разница
      0   0.30    0.6948   -0.0073     0.6773   -0.0020
      1   0.30    0.6842   -0.0040     0.6826   -0.0024
      2   0.45    0.7086   -0.0100     0.6829   -0.0006
      3   0.30    0.7082   -0.0060     0.6771   -0.0027

**All four held-out seeds improve, in both the near layer and overall**, and t lands at 0.30 three
times out of four. The consistency is the whole of the evidence: the near-layer gain averages
-0.0068 and the overall gain -0.0019, both at or under the 0.007 floor. Selecting t on the near
layer instead reaches -0.0079 there but reverses on seed 1.

The weights this implies are per-enzyme 0.23, pooled 0.38, GP 0.18, ridge 0.21. **It is not going
into `src/submit.py`.** Four consistent seeds are worth something, but the size is below the floor
this repository set for itself, and the case for it rests on the layer standing in for the test
set — which items 110 and 111 have just shown we cannot explain. Changing the submission on an
unexplained effect measured under its own noise floor is the trade this log exists to refuse. It is
recorded so that the decision is visible rather than silently taken either way.

**113. Every conclusion in the log survives being re-scored in the test set's regime — and item
97's gap was overstated.** `verify/k25_reweight.py`. This is the check items 107 and 108 made
necessary: if the models rank differently where the test sits, how many of the log's conclusions
were reached in a regime that does not apply to them.

The shift here is identifiable in a way the label shift of item 94 is not. The variable is
similarity to the training set, computed from structures, and the test structures are published:
p_test(s) is counted rather than inferred, with no posterior and no model-choice systematic.

One confound had to be removed first. An out-of-fold row's similarity is a maximum over four
fifths of the training set and a test row's over all of it, and a maximum over a larger reference
set is larger for free. Taking the test similarities against random four-fifths subsets instead:

    OOF медиана 0.435, тест медиана 0.566

Item 97's 0.587 is a maximum over the whole training set, and setting it against the out-of-fold
0.435, which is a maximum over four fifths, compares two different reference sizes. Matched at
four fifths the test sits at 0.566 and **the gap is 0.131 rather than the 0.152 that pairing
implies.** Item 22 already carried the matched control and was not caught out by this: its
leave-one-out figure of 0.450 is taken against about 4900 molecules, the same as the test's 0.587,
and that gap is 0.137. So the gap is about 0.13 by either size-matched comparison, and what was
too large is the mixed pairing, not any single number. Items 107, 108 and 112 use the layer
boundaries rather than the median and none of their conclusions turn on this.

Weights by similarity bin run from 0.00 below 0.35 to 8.00 above 0.70, and the price is steep:

    фермент       n     ESS   ESS/n
    CYP1A2     1412     350   24.8 %
    CYP2C9     1285     471   36.7 %
    CYP2D6     1493     309   20.7 %
    CYP3A4     2335     918   39.3 %

**Two enzymes fall below the one-third rule this file states for itself**, so for CYP1A2 and
CYP2D6 what follows is directional and not decisive. Said here rather than in a footnote, because
the rule was written before the numbers were seen.

*(Item 123 shows that rule cannot be met here by any estimator, and rewords what the low numbers
mean. The caution stands; the framing of it as a shortfall does not.)*

Six of the seven ablations keep their arm ordering exactly:

    файл                  вывод при равных весах        под тест
    признаки              FP+DESC+MECH > FP+DESC > FP   тот же
    пул                   пул > независимо              тот же
    слабые                лес > гребневая > kNN         тот же
    форма                 с формой > без формы          тот же
    FCFP                  оба > ECFP > FCFP             тот же
    агрегаты              флаг > ... > база             ПОРЯДОК ИЗМЕНИЛСЯ

The one reordering is the aggregator file, whose arms differ by about 0.0045 — under the floor —
so a reshuffle there is what noise looks like, not a finding.

The sizes are the interesting part. The two large effects roughly double in the test's regime and
the two small ones disappear:

    эффект                   равные веса   под тест
    механистический блок         +0.0163    +0.0282
    пул                          +0.0142    +0.0291
    FCFP оба                     +0.0035    +0.0024
    форма                        +0.0033    +0.0001
    флаг агрегатов               +0.0047    +0.0001

So the regime question, which items 107 and 108 opened as a threat to the whole log, closes as a
reprieve for the conclusions and a warning about the marginal ones: **nothing we concluded is
wrong in the test's regime, and two of the things we concluded are worth nothing there.** The
shape block and the aggregator flag are both retained on a gain that vanishes at the similarity
the test set actually sits at.

**114. Quantile regression onto the band edges is a monotone reparameterisation of the model we
already have — and the metric is simpler than we have been describing it.** The proposal was two
quantile regressions onto `lo` and `hi` as targets in their own right, taking the action as the
median of the pooled predicted edges. Its argument was that this is a *different* function of the
features rather than a monotone correction, so unlike item 93 it could reorder compounds within a
group.

Item 105 makes that checkable without fitting anything. Isotonic regression from the label alone
onto the band width:

    фермент       n   rho(y, ширина)   R2 изотоники   ост. std   std доли метки в полосе
    CYP1A2     1412           -0.885          0.963      0.125                     0.057
    CYP2C9     1285           -0.899          0.928      0.155                     0.063
    CYP2D6     1493           -0.558          0.955      0.143                     0.087
    CYP3A4     2335           -0.928          0.970      0.159                     0.088

**The width is a deterministic function of the label to within three per cent of its variance**,
and the label sits at a nearly fixed relative position inside the band — the fraction (y - lo)/w
has a standard deviation under 0.09. So `lo` and `hi` are the label plus and minus a function of
the label. A model of the edges is a monotone reparameterisation of a model of the centre, and
item 93's rule closes the proposal along with it. Five minutes, nothing fitted.

The side effect is worth more than the closure. If the half-width is a known decreasing function
of potency, then **ST-RAE is a potency-weighted absolute error**: the forgiveness threshold is
about 1.25 pIC50 below a label of 3.5 and about 0.13 above 4.6. That is a much simpler statement
than "soft-thresholded error against a confidence band" and it belongs in the document.

**115. Where the score is actually made, and a correction to the sentence above.** Reading item
114 as "get the potent compounds right, the weak ones are nearly free" is the obvious inference
and it is wrong. Share of the ST-RAE numerator by potency quartile, on the submitted ensemble:

    фермент   кв.1 слабые   кв.2    кв.3   кв.4 сильные   полупорог кв.1   полупорог кв.4
    CYP1A2         37.1 %   9.6 %  12.5 %         40.8 %            0.537            0.098
    CYP2C9         27.9 %   9.9 %  12.8 %         49.5 %            0.546            0.128
    CYP2D6         29.2 %   9.6 %  10.2 %         51.0 %            0.245            0.105
    CYP3A4         16.0 %  20.9 %  22.4 %         40.7 %            1.233            0.069
    среднее        27.6 %  12.5 %  14.4 %         45.5 %

**The penalty is U-shaped, not monotone.** The potent quartile does dominate at 45.5 %, for the
reason item 114 gives — there is almost no forgiveness there. But the weak quartile is second at
27.6 %, and the two middle quartiles together supply only 26.9 %. The weak end is not free in
practice, because item 105 measured our median absolute residual there at 1.27 to 1.67 pIC50, and
that overruns even a threshold of half a log unit.

The forgiveness at the weak end differs fivefold between enzymes — a half-threshold of 0.245 on
CYP2D6 against 1.233 on CYP3A4 — and the share tracks it: CYP3A4, the most forgiven, is the only
enzyme whose weak quartile is cheap at 16.0 %, while CYP2D6 with the tightest bands pays 29.2 %
there despite being the enzyme where every model does worst.

So the effort goes to both ends and the middle can be left alone, which is a different
prescription from the one the previous paragraph implied.

It also raises a concern about item 105's dead-zone objective, and the concern is recorded here
before that run reports. Item 114 found the band width is a function of the label, so "wide band"
and "weak compound" should be nearly the same set, and they are: the Jaccard overlap between the
rows with a band above 1.0 pIC50 and the equally many weakest rows is 0.823, 0.738, 0.897 and
0.930 across the four enzymes. The dead zone is therefore close to a potency-based down-weighting
of the weak end — **and the weak end supplies 27.6 % of the penalty, not the near-nothing the
band widths alone would suggest.** Withdrawing effort there may well cost more than it saves.

If the dead zone loses, that is the reason, and it will be a better-understood loss than a win
would have been. The residual disagreement between the two sets is largest on CYP2C9, at a
Jaccard of 0.738, so that is the enzyme where the two readings can be told apart.

**116. The doubling in item 113 holds for the mechanistic block and is unverifiable for pooling.**
Item 113 reported that the two large effects roughly double in the test's regime, as a macro
average, and in the same breath stated that CYP1A2 and CYP2D6 fall below the one-third effective
sample size rule. Those are the same two enzymes the effects live on, so the macro average may be
reporting a doubling that exists only in the cells the file said not to trust. Rank needs no affine
pair, so the breakdown costs nothing:

    фермент   ESS/n   мех. блок равн.   под тест   пул равн.   под тест
    CYP1A2   24.8 %           -0.0012    +0.0081     +0.0127    +0.0564
    CYP2C9   36.7 %           +0.0138    +0.0283     +0.0057    -0.0053
    CYP2D6   20.7 %           +0.0517    +0.0739     +0.0421    +0.0636
    CYP3A4   39.3 %           +0.0009    +0.0024     -0.0039    +0.0017

    макро по всем четырём      мех +0.0163 -> +0.0282   пул +0.0142 -> +0.0291
    только ESS выше трети      мех +0.0074 -> +0.0154   пул +0.0009 -> -0.0018

**The mechanistic block's doubling survives the restriction** — +0.0074 to +0.0154 on the two
enzymes whose weights can be read. The claim stands on the cells it is allowed to stand on.

**Pooling's does not.** Restricted to CYP2C9 and CYP3A4 it is +0.0009 uniform and -0.0018
weighted, and CYP2C9 flips sign outright. The whole of the macro doubling comes from CYP1A2 and
CYP2D6, at effective sample sizes of 24.8 % and 20.7 %.

The right conclusion is narrower than either "confirmed" or "refuted". Pooling helps CYP1A2 and
CYP2D6 and does nothing for CYP2C9 and CYP3A4 **under uniform weights already** — that much is
solid and has four seeds behind it. Under test weights the helped enzymes appear helped more, but
those are precisely the enzymes whose reweighted estimate is too thin to read, so the doubling is
**unverifiable rather than established**. Item 113's ordering claim is untouched; its magnitude
claim for pooling is withdrawn.

This is also the pattern item 111 predicts. Pooling helps CYP2D6, which correlates with nothing,
and CYP1A2, the next least correlated; it does nothing for the CYP2C9-CYP3A4 pair that correlates
at 0.688. The per-enzyme breakdown and the mechanism disagree with each other nowhere, which is
some comfort about both.

**117. The concatenated encoder loses by rank, and item 68 was never updated to match the
document.** The chemprop embedding concatenated to FP+DESC+MECH was recorded in item 68 at a raw
macro ST-RAE of 0.7589 against 0.7673, read as a gain of 0.0084 and used to overturn item 61's
conclusion that the representation was never the bottleneck. Item 77 then found that raw ST-RAE is
rewritten by the affine pair, and the document's table of the five collapsed findings already
carries this row at +0.011 after the pair. **The journal entry was left saying the opposite of the
document's own table**, and has been marked withdrawn above.

Re-reading the saved predictions under the criterion that decides things now adds the part neither
had:

    фермент   сырое база   сырое +EMB   пара база   пара +EMB   ранг база   ранг +EMB
    CYP1A2        0.8786       0.8562      0.8128      0.8184      0.4957      0.4733
    CYP2C9        0.6908       0.6998      0.6559      0.6833      0.5972      0.5869
    CYP2D6        0.9803       0.9412      0.9052      0.8966      0.4027      0.4217
    CYP3A4        0.5194       0.5383      0.4860      0.5071      0.7646      0.7485
    МАКРО         0.7673       0.7589      0.7150      0.7263      0.5650      0.5576

    сырое -0.0084   пара +0.0113   ранг -0.0075   ранг под тест -0.0182

**The raw number and the pair number carry opposite signs**, which is the sharpest instance of
item 77 in the log — an apparent gain of 0.0084 that is a loss of 0.0113 once the post-processing
runs. By rank it loses 0.0075, and in the test set's regime it loses 0.0182, more than twice as
much. Item 61 is reinstated: on this checkpoint the representation was not the bottleneck, as a
substitute or as an addition.

Two things follow beyond the row itself. The first is procedural: item 77 invalidated five
findings, the document was corrected, and **one journal entry was left standing in contradiction
to it for a month**. Nothing catches that automatically; the only defence is re-reading old items
against new criteria, which is what this one did.

The second bears on item 10 of the outside reading, which proposes pretrained representations as
features on the ground that item 96 forbids carrying a metric result over to features. That
remains true, but the transfer has now been measured and it is negative: chemprop as a metric
placed second at 0.222, and as a concatenated feature block it costs 0.0075 of rank and 0.0182 in
the test's regime. A larger pretrained model is a different question, but it starts from a worse
prior than the proposal assumed.

**118. Retracted: the screening head is built, measured, and already in the log at -0.0264.**
This item first read item 5 of the outside reading as an unbuilt proposal, checked its
preconditions, found them all met and recommended building it. Every operative part of that was
wrong, and the file it should have checked first was this one.

`src/trunk.py` implements exactly the proposal — shared trunk, a pIC50 head and a screening head,
`lambda_scr` as the switch, with the control being the same architecture at `lambda_scr = 0` so
the parameter count and the weight initialisation are identical. It has been run on four seeds.
**Item 79 measured the channel at -0.0264 after the affine pair, with the sign holding on all four
seeds, t = -10.87, p = 0.002** — the largest single surviving intervention in this repository. Item
79 also found the second half: the trunk carrying the channel reaches 0.7149 against boosting's
0.7155, sign alternating, p = 0.72. **Parity, not superiority**, which is why the trunk is shelved
and why the channel is not.

So there is nothing to build and nothing to greenlight. The question item 5 actually poses, given
all of that, is a different and narrower one: the channel carries rank information the boosting
does not have, and the vehicle carrying it only ties. Whether that information can be moved into
the boosting — as an ensemble member, or as the trunk's latent concatenated to the feature matrix
— is open, and item 101's rule says the concatenated arm is the one to run.

A counting error of my own, corrected. This item claimed 4335 molecules carry a screening reading
and no curve, and called that a near-doubling of the molecule set. The set was built as a union of
per-enzyme complements, which counts a molecule that has a CYP1A2 curve but no CYP2C9 curve as
"having no curve". **The correct number is one.** Every molecule in the screening file but one
carries at least one pIC50, exactly as `trunk.py`'s docstring has said since it was written: the
screen is extra *columns* on rows already present, not extra rows.

What survives from `verify/k26_screen.py` is one measurement that is genuinely new, and it is a
characterisation rather than a recommendation. Per enzyme, the compounds that received a curve are
selected on their screening effect, and the severity of that selection differs sharply:

    фермент   медиана |log2fc| с кривой   без кривой   rho(есть кривая, |log2fc|)
    CYP1A2                        1.554        0.209                        0.655
    CYP2C9                        0.815        0.602                        0.152
    CYP2D6                        1.738        0.395                        0.819
    CYP3A4                        1.437        1.932                       -0.151

CYP2D6 and CYP1A2 are severely selected — their curves went to the compounds that moved the
screen, and the 2883 and 2964 without curves are the ones that did not. CYP2C9 is nearly
unselected. **CYP3A4 is selected in reverse**: its un-curved compounds show the *larger* screening
effect. The per-enzyme profile of the trunk's screening channel should be checked against this
ordering, since a channel that repairs a label selection ought to pay most where the selection is
worst. That check has not been run and is not claimed here.

Structurally the un-curved compounds sit only slightly outside: median similarity to the curve set
0.296 against 0.342 within it. So whatever this is, it is a selection on labels and not on
chemistry.

Two things to carry forward from the retraction itself. The precondition discipline that closed
items 109, 114 and worked well elsewhere has a failure mode: **it checks whether an idea could
work and not whether it has already been tried**, and this log is now long enough that the second
question needs asking first. And the union bug is the kind that produces a large, satisfying number
in the direction one is hoping for; 4335 was never sanity-checked against the 4905 rows sitting in
the same table.

**119. Four seeds for the marginal effects, and the third seed-0 mistake of the day.**
Item 113 compared four-seed uniform gains against seed-0 reweighted ones and concluded that the
shape block and the aggregator flag "go to zero" in the test's regime. The two halves of that
comparison were not the same measurement. Weighted rank needs no affine pair, so the four-seed
version costs nothing, with the weights rebuilt per seed because the out-of-fold similarity
depends on the split:

    эффект        сидов   равн. ранг   под тест   разброс по сидам
    форма             4      +0.0034    +0.0014             0.0045
    флаг агр.         1      +0.0046    +0.0002             ---
    FCFP оба          4      -0.0004    +0.0032             0.0016
    мех. блок         4      +0.0163    +0.0313             0.0028

**The mechanistic block's doubling is confirmed on four seeds** — +0.0163 to +0.0313, with a
seed-to-seed spread of 0.0028, an order below the effect. Together with item 116's per-enzyme
check this is the one large claim of the day that survives every restriction placed on it.

**The shape block's collapse is not established.** Its weighted estimate is +0.0014 with a
seed-to-seed spread of 0.0045 — larger than the estimate. That is not "goes to zero", it is "the
measurement cannot tell +0.0034 from 0", and item 113's wording is corrected to that. The flag
still rests on one seed and is left there.

**FCFP is refuted outright**, and by its own completed run. Item 101 was written from seed 0,
where the concatenation won 0.0034 of rank. Seeds 1, 2 and 3 all reverse it, and the four-seed
average is -0.0004 — nothing. That item has been rewritten above.

This is the third time in one working session that a conclusion was drawn from partial output and
had to be withdrawn when the run finished, after "one column beats the whole shared latent" from
two seeds of four and "the state wins" from three rows of four. The pattern is specific enough to
name: **the first seed is read while the rest are still computing, and it is read as though it
were the answer.** Nothing in the tooling encourages waiting, and the ablation scripts print
per-seed lines precisely so that progress is visible. The remedy is a rule rather than more care —
**no item is written from a run that has not printed its aggregate table**, and the three
withdrawals above are the argument for it.

It is worth noting what the same discipline bought when it was applied. Items 109, 114 and 118
each settled a multi-day proposal in under an hour by measuring a precondition first, and none of
them has needed correcting. The failures and the successes of the day differ by exactly one thing:
whether the number was finished before it was believed.

**120. The trunk is stranded no longer: as a fifth ensemble member it is worth more than the
ridge.** `verify/k27_trunkens.py`. Item 118 left the screening channel in an awkward place — worth
-0.0264 as an intervention (item 79), sign holding on four seeds, but carried by a model that only
ties with the boosting, so the vehicle was shelved and the channel with it.

A tie is exactly the situation item 100 was about. The ridge is *worse* than the boosting alone
and improves the ensemble anyway, because it errs smoothly where the boosting errs locally. A
model that ties while coming from a different family is a stronger prior for that than the ridge
was, and every prediction needed was already on disk. Four seeds, the two-head trunk at
lambda = 3, clipped to the enzyme's label range plus or minus two units exactly as
`src/trunkdose.py` does:

    состав          пара      ранг
    четыре        0.6819    0.6009
    пять          0.6758    0.6063
    ствол один    0.7153    0.5646

    добавление ствола   пара -0.0061  (t = -17.68, p < 0.001, знаков 4/4)
                        ранг +0.0054  (t = +16.06, p = 0.001,  знаков 4/4)

**And it helps every enzyme, four seeds out of four in each:** CYP1A2 +0.0062, CYP2C9 +0.0097,
CYP2D6 +0.0029, CYP3A4 +0.0030. No other member in this log improves all four.

The size is the point. Item 100 admitted the ridge to the submission at -0.0026 with p = 0.0005;
**the trunk is 2.3 times that and equally consistent**, and its pair gain of -0.0061 sits just
under the 0.007 floor while its per-enzyme signs are unanimous. By the rule item 100 set — that
membership is settled by measurement — the trunk belongs in the ensemble.

**The attribution was flagged as unearned before the control ran, and the control earns it.** The
trunk is also the only member from a different model family, and item 100 established that family
diversity alone buys something, so a gain here could have been the multilayer perceptron rather
than the channel. The separating arm is the trunk at `lambda = 0` — same architecture, same
parameter count, same initial weights, screening head receiving no gradient:

    состав                пара      ранг
    четыре              0.6819    0.6009
    пять (lam = 3)      0.6758    0.6063
    пять (lam = 0)      0.6811    0.6001
    ствол один (lam=3)  0.7153    0.5646
    ствол один (lam=0)  0.7416    0.5297

    семейство само по себе   пара -0.0008 (t = -3.28, p = 0.046, знаков 4/4)
    вклад канала             пара -0.0053

**The family contributes -0.0008 of pair and -0.0008 of rank — that is, by rank it slightly
hurts.** The channel supplies -0.0053 of the -0.0061 and the whole of the rank gain. Standalone
the same split is starker: the trunk's rank goes from 0.5297 without the channel to 0.5646 with
it. So the credit belongs where item 5 of the outside reading put it, and it is the screening
supervision rather than the architecture that the ensemble is buying.

Two things make the measurement itself trustworthy rather than merely large. `trunk.py` calls `butina_folds` from
`cypsplit.py` with the same seed, so the trunk's out-of-fold predictions, the boosting's, and the
affine pair are all on one split; nothing is being averaged across incompatible partitions. And
the clip is not a thumb on the scale — item 79 records that without it seed 0 is one compound
predicted at -360, and the document's own model comparison clips for that reason.

What this costs to act on is the honest caveat. The four current members are all fitted inside
`src/submit.py`; the trunk is a torch model that would have to be trained on the full training set
and run on the 750 test structures, which is new plumbing rather than a new average. The
measurement says it is worth building; the building is not done.

It also answers the question item 118 was left holding, and answers item 5 of the outside reading
in the affirmative after all. The screening channel does transfer out of the shelved model; it
transfers by averaging rather than by concatenation; and the vehicle that only ties with the
boosting is worth carrying anyway, for the one thing inside it that the boosting has no route to.
The remaining unexplored path from item 118 — the trunk's latent as a feature block — is no longer
the only one and is no longer urgent.

**121. The fifth member is wired in, behind a flag, with three guards.** Item 120 measured the
trunk as worth -0.0061 of pair and +0.0054 of rank and item 120's control attributed almost all of
it to the screening channel. Acting on that needed the trunk applied to the 750 blinded structures,
which the repository could not do: `src/trunk.py` only ever cross-validated.

`trunk.fit_predict_test` is that path and it is deliberately not new code. The test rows are
appended with their own fold index and NaN in both target blocks, so `masked_mse` gives them no
gradient and every statistic `run_fold` computes — the feature standardisation included — is still
taken over training rows alone. `run_fold` is then called unchanged, which makes the claim
checkable rather than asserted: holding out fold 1 through the ordinary route and through the
appended-rows route must give the same numbers. **Maximum absolute difference: 0.000e+00.**
`src/trunk.py --check-test-path` runs it.

Three guards, because each of these fails silently rather than loudly:

- **The feature split.** The trunk takes FP, DESC and MECH separately because it puts the
  fingerprint through `log1p`, and `submit.py` carries them hstacked. The widths are read from
  `feats.npz` rather than written as numbers; all three blocks are verified to round-trip.
- **The device.** The saved out-of-fold predictions were computed on `mps`. The affine pair is
  fitted on those and applied to test predictions computed here, so the two halves must not land
  on different arithmetic. The device is taken from the file's own meta block, with a warning if
  it is unavailable.
- **The split.** This is the only member read from a file rather than recomputed. If the split
  ever moves, every other member follows it and this one silently stays on seed 0's folds, at
  which point its "out-of-fold" predictions are nothing of the kind. The fold digest is checked
  against the same golden `2d93c19815e14261` that `tests/test_split.py` pins, and the guard is
  verified to fire on a deliberately altered fold vector.

`src/submit.py --mode ансамбль5` runs end to end and writes both files: 750 rows, no missing
values, predictions between 2.0 and 6.1. **The default is not switched.** The measurement supports
the member; the choice of what to submit belongs to the team, and the mode makes the change one
flag rather than one edit.

**122. The document build stopped manufacturing conflicts.** Not a model finding, but this
repository is worked from four machines and its binary artefacts cannot be merged, so it is the
same class of problem as the split digest.

Rebuilding the document rewrote all ten figure PDFs and the main PDF every time. Checked rather
than assumed: `docs/tex/fig/nnsim.pdf` before and after a rebuild is 31185 bytes both times and
identical outside a single `CreationDate` string. Ten unmergeable binaries were entering every
diff carrying no information at all.

`docs/tex/figs.py` now writes figures with `metadata={"CreationDate": None}` and `docs/build.sh`
pins `SOURCE_DATE_EPOCH`. Verified the way the rest of this file verifies things — two consecutive
builds, byte-compared: the figures match, and so does the main PDF. **A rebuild that changes
nothing is now no diff at all**, which means a rebuild that does change something is entirely
signal.


**123. The low effective sample size in item 113 is the phenomenon, not a defect of the estimator.**
Item 113 reweighted the out-of-fold rows to the test set's similarity distribution, got ESS/n of
20.7 to 39.3 per cent, and disclaimed two enzymes for falling under the one-third rule. That reads
as though a better density-ratio estimate would rescue them. It would not, and the bound is worth
computing before anyone tries.

For self-normalised importance weights the effective fraction converges to 1/(1 + chi2) where chi2
is the chi-square divergence between the two distributions, whatever produces the weights. Measured
on the same bins:

    chi2(тест || OOF) = 2.838      потолок ESS/n = 1/(1 + chi2) = 26.1 %

    обрезка веса     ESS/n: 1A2    2C9    2D6    3A4
    4                      26.3   37.3   22.0   41.6
    8  (использована)      24.8   36.7   20.7   39.3
    без обрезки            20.5   35.3   17.0   32.8

**The weights already sit at the bound.** The unclipped column is the honest one and it is the
lowest; clipping buys apparent ESS by biasing the ratio, which is a trade and not an improvement.
Per-enzyme numbers straddle the 26.1 per cent because that figure is the divergence over all rows
and each enzyme's subset has its own overlap — CYP3A4's rows are already the most test-like and
CYP2D6's the least, which is the same ordering everything else in this group produces.

So the wording in item 113 was wrong in a way worth naming. The one-third rule is a rule about
when a reweighted estimate can be trusted; here it **cannot be satisfied by any method**, because
the test set sits where our cross-validation has almost no mass. That is not a shortfall of the
measurement — it is the measurement. The right statement is that any test-regime estimate on
CYP1A2 and CYP2D6 carries about four times the variance of the unweighted one and always will, and
that is a fact about the challenge rather than about our arithmetic.

It also closes a line of work before it starts. A smoother ratio — logistic discrimination between
the two sets instead of histogram bins — would reduce binning variance and cannot move the bound.
Not worth the afternoon.

**124. The verification layer had never been verified, and two of its scripts were lying.**
A regression pass after a day of edits, run only to check nothing had broken, found four defects in
the checking code itself. None of them is about the model; all of them are about whether this file
can be believed.

**`verify/f1_formulas.py` has never run on the numpy this repository actually uses.** It calls
`np.trapezoid`, which exists only from numpy 2.0, and the environment has 1.26.4. The call has been
there since the first commit. `pyproject.toml` states that numpy is deliberately unpinned because
"1.26.4 and 2.5.2 both reproduce exactly" — true of the pipeline, false of this script, and the
script is the one the table at the top of this file calls "44 numerical checks of every formula in
the document". It crashed before reaching any of them. Now uses whichever name the installed numpy
has.

**Once it ran, one of the 44 failed — and the check was wrong, not the formula.** `I -> E as
C -> inf` was tested at C = 100 M with a tolerance of 1e-9. The residual is exactly
10^(h(pC - pi)), which at C = 100 comes to 6.3e-8: sixty times the tolerance. The limit is correct
and the test point could never satisfy it; it needs C above about 10^3.8 M. Tested at 1e6 M it
passes. **44 of 44 now.**

**`verify/f6_data.py` reported three discrepancies that do not exist.** It "reconciles 21 numbers
from the document against the source tables", and three of the 21 were the script's own fault:

- it looked for a column named `concentration` where the file has `concentration_M`, and printed
  "no such column" as a mismatch against a value that is in fact exactly what the text says;
- its expected interquartile ranges had CYP1A2 and CYP2D6 **transposed** relative to
  `docs/tex/s14.tex`. The document says 1.00 and 0.83, the data give 1.00 and 0.83, and the script
  expected 0.83 and 1.00. Two reported errors, both phantom.

**21 of 21 now.** Nothing in the document moved; the checker was stale.

The meta-point is uncomfortable and worth stating plainly. This file's authority rests on those
scripts, and the two whose whole purpose is reconciling the document against the data were between
them producing one crash and three false alarms. A false alarm in a checker is worse than no
checker: it trains the reader to discount the output. And the failure mode is specific — **nothing
runs the verification scripts on a schedule.** `make verify` exists, CI runs only
`tests/test_split.py`, and f1 has been dead since the first commit without anyone noticing.

The cheap fix is not more care. `make verify` and `make verify-extra` should run somewhere that
reports, or at minimum the fast ones should join the test suite; f1 takes forty seconds and f6
twenty, and both would have caught all four of these on the day they appeared.

**125. Doubling the pooled table with maximally correlated rows buys nothing, so sample size is
not what pooling is doing.** `src/ablpool.py --arms "пул+TDI"`, two seeds. The pre-incubation arm
covers the same four enzymes under a different condition and adds 6538 rows, taking the pooled
table from 6525 to 13063 — more than double, and about as correlated with the direct-inhibition
labels as any rows could be, since they are the same molecules on the same enzymes.

Under the sample-size reading this had to be the largest gain anywhere in this log. Rank, averaged
over seeds 0 and 1:

    рука          MACRO      1A2      2C9      2D6      3A4
    независимо   0.5627   0.4934   0.5940   0.3999   0.7636
    пул          0.5778   0.5134   0.6027   0.4385   0.7566
    пул+TDI      0.5751   0.5136   0.6042   0.4302   0.7526

**It adds nothing.** Against plain pooling it is -0.0027, which is under the floor and therefore
not a loss either — the honest statement is that more than doubling the table changed the answer by
less than the seed does. Per enzyme it is +0.0002, +0.0015, -0.0083, -0.0040: no enzyme gains.

That is decisive against sample size and it was the arm designed to be decisive. Pooling four
enzymes at 6525 rows is worth +0.0151 of rank; adding 6538 more rows of the same enzymes is worth
nothing. Whatever pooling supplies, it is not quantity, and the argument no longer rests on the
four-point inversion of item 111 alone.

Two mechanisms down by measurement — borrowing neighbours in item 110, shared function in item
111 — and now quantity as well. Contrast is what is left standing, and `src/ablpool.py --arms
"пул слепой"` tests it directly by zeroing the enzyme indicator: shared function and quantity both
survive without it, contrast cannot exist without it.

**126. The exact Bayes action under ST-RAE loses to the two-parameter pair it was meant to
replace, and the argument for it has a hole worth naming.** `src/bayesact.py`, seed 0.

The proposal was the strongest thing an outside reading has put forward, because it claimed to
escape the ceiling that killed five earlier findings *provably* rather than empirically. The chain
was: item 114 gives the band as a near-deterministic function of the label, so the loss
max(0, lo - p, p - hi) is known as a function of (prediction, truth); under a posterior for the
truth the optimal action is a one-dimensional convex minimisation per compound; that action
depends on the posterior's **spread**, and item 93's rule only reaches corrections that are
monotone in the point.

The premise checks out. The best monotone function of the point prediction explains 0.110, 0.114,
0.100 and 0.080 of the spread's variance across the four enzymes — about nine tenths of the spread
is information the point does not carry. That was measured before the file was written.

The conclusion does not follow, and this is the hole:

    **"the correction depends on X, and X is not a function of the point" does not imply
    "the correction is not a function of the point".** The dependence also has to be strong
    enough to survive the loss's curvature. Here it is not.

Measured, with the falsification criterion written into the script before it ran:

    рука                  ST-RAE      ранг
    сырое                 0.6979    0.6029
    аффинная пара         0.6793    0.6014
    байесово действие     0.6893    0.6012
    байесово + пара       0.6863    0.5998

    монотонность действия к входу:  1A2 +0.985  2C9 +0.996  2D6 +0.972  3A4 +0.997
    средний сдвиг:                      +0.127      +0.125      +0.118      +0.178

**The action comes out 97 to 99.7 per cent monotone in the point**, and reduces in practice to a
shift of about +0.13 — which is precisely what the pair's shift parameter is. The pair fits that
shift by direct empirical minimisation against the *true* bands on the training folds; the Bayes
action derives it through a model of the band, a model of where the label sits inside it, and a
model of the posterior, accumulating three approximations. The pair wins by 0.0100.

Three checks were run before closing it, because none of the three failures should be blamed on
the implementation.

*The cliff mechanism does not exist.* The obvious rescue is that the action should bite where the
posterior straddles the width cliff — forgiveness is about 1.25 below a label of 3.5 and 0.13
above 4.6 (item 115). Deviation from the best monotone fit, by predicted potency: 0.043 / 0.056 /
0.061 on CYP1A2 and 0.010 / 0.038 / 0.050 on CYP2D6. It is flat, and if anything largest where the
band is *narrowest* — the opposite of the mechanism. CYP2D6 being both the least monotone enzyme
and the only one where the action wins (-0.0027) is a coincidence.

*The construction is not degenerate.* Within a spread bin every compound receives the same
residual pool, so the action there is "point plus a constant", and the whole non-affine effect
lives in the differences between bins. Those differences are real — a range of 0.07 to 0.20 across
five bins — so the design does have something to work with. But the per-bin offsets are not
ordered by spread (CYP1A2: +0.080, +0.020, **+0.190**, +0.070, +0.110), which is what estimation
noise looks like rather than signal.

*The residual subsample was mine, not the idea's.* `NRES = 150` was a speed choice. Using every
residual and sweeping the bin count:

    корзин   байесово   пара
      3       0.6864   0.6793
      5       0.6866   0.6793
     10       0.6886   0.6793
     20       0.6944   0.6793

Dropping the subsample does help — 0.6893 to 0.6866 — which confirms the shortcut was adding
noise. And the trend settles it: **the fewer bins, the better the action performs, and one bin is
the affine shift.** Every degree of per-compound adaptivity the construction adds costs more in
estimation variance than it buys in decision quality. Extrapolated to its own optimum, the
principled rule converges to the thing it was built to beat.

The fourth arm is the summary. A pair applied on top of the Bayes action still improves it, from
0.6893 to 0.6863 — so the action has not replaced the post-processing, it has approximated it
badly. Had the construction been right, that arm would have found nothing left to fix.

This closes the decision layer. It does not touch the rest of the hierarchical proposal, and the
diagnosis is specific enough to be useful there: a layer whose benefit is per-compound has to
clear an estimation-variance bar that two well-fitted global parameters set surprisingly high.

**127. The selection layer is buildable on three enzymes and structurally impossible on the
fourth — which is the one that needs it most.** Precondition check for the MNAR correction, run
before building anything.

The labels are missing not at random: a compound received a dose-response curve if it looked
active in the single-concentration screen. The standard correction is a weight of 1/P(curve | x),
and it is estimable here because the screening column is complete and it is known who received a
curve. The target population is the full 4905 training compounds — a compound unlabelled for an
enzyme *is* an unselected one, and its screening readout exists.

The first diagnostic I ran was wrong and is recorded because the error is instructive. I computed
the chi-square divergence of the weights over the **selected** compounds, as item 123 does for the
similarity reweighting, and read CYP2D6's 0.017 as "cheap". It is not cheap; it is degenerate.
Chi-square over the selected set measures the spread of the weights among them, not whether they
can reach the target at all. When P(selected | x) goes to zero over a whole region, no weight
recovers it — that is a **positivity** violation, and an effective sample size cannot see it.

The correct check is overlap:

    фермент   отобрано   недостижимая доля цели   перекрытие по скрину
    CYP1A2       1412                     1.3 %                 98.7 %
    CYP2C9       1285                    15.4 %                 77.5 %
    CYP2D6       1493                    99.7 %                  0.3 %
    CYP3A4       1805                     1.4 %                 98.1 %

("Unreachable" is the share of the target population whose propensity falls below the 1st
percentile of the selected. Counts are over compounds that also carry a screening readout, 89.2 %
of rows.)

**On CYP2D6 the two populations are disjoint.** Median screening log2fc is -1.738 among the
selected and +0.376 among the unselected; the selection rule is effectively a hard threshold and
99.7 % of the target sits below anything the labelled set contains. Propensity weighting cannot be
done there at all. A parametric selection model can be *written* for it, but every number it
produces outside the support is the model's assumption rather than the data's, with nothing to
check it against — the Heckman situation, and it is known to be at its most fragile exactly when
overlap is this poor.

On CYP1A2 and CYP3A4 positivity is comfortable and the layer is worth building. CYP2C9 sits in
between at 15.4 %.

Two things this changes. First, the layer's coverage is inverted with respect to need: it is
available on the enzymes we already predict best, and unavailable on CYP2D6, which is both our
worst (ST-RAE 0.86 against 0.46 on CYP3A4) and the one whose selection is most severe. Second, the
mechanism is not common across enzymes and should not be modelled as if it were. The rank
correlation between the screening readout and getting a curve is **-0.642, -0.151, -0.821 and
+0.151** — on CYP3A4 the sign flips, and its unselected compounds inhibit *more* strongly than its
selected ones. Whatever chose CYP3A4's curves, it was not activity in the screen.

**128. The oracle gate: one measurement that bounds every post-hoc per-compound proposal, made
before any of them is built.** `verify/k30_oracle.py`. This is a **procedure**, and it is recorded
as one.

Item 126 closed a proposal that claimed to escape item 77's ceiling *by construction*. The general
form of that failure is worth more than the instance: **item 77's ceiling is not a claim about a
class of functions, it is a claim about magnitudes.** Showing that a correction lies outside the
affine family says nothing on its own; what decides is whether its non-affine part survives the
curvature of the loss. So a proposal of the form "a per-compound correction applied to the
predictions" needs an a-priori bound on the size of its non-affine part, not a proof that one
exists. This file computes that bound once, for all such proposals at once.

Take the correction a perfect per-compound layer would apply — move each prediction to the nearest
point of its own true band, `c = clip(p, lo, hi) - p`. It drives the numerator to exactly zero, so
it is unimprovable, and it is unachievable by construction since it uses the true label. Anything
a real layer does is a subset of it.

    фермент    пара   монот. потолок   бюджет   corr(оракул, остаток)
    CYP1A2   0.7777          0.7668   0.0109                   0.948
    CYP2C9   0.6097          0.6003   0.0094                   0.885
    CYP2D6   0.8686          0.8588   0.0098                   0.935
    CYP3A4   0.4716          0.4714   0.0002                   0.871
    МАКРО    0.6819          0.6743   0.0076

Two things fall out, and the second is the general one.

**The whole budget for better monotone post-processing is 0.0076 macro** — at the noise floor,
and on CYP3A4 it is 0.0002, meaning the affine pair is already at the monotone optimum there. The
ceiling is computed by majorise-minimise, isotonic regression onto `clip(m, lo, hi)` iterated,
which is the same reduction `src/abldead.py` uses applied to the monotone class instead of to a
learner. Fits are in-sample on purpose: that is generous to the monotone class, which is the
conservative direction for a gate.

**The oracle correction is the model's own residual, at a correlation of 0.87 to 0.95.** That is
the statement that closes the class, and it is stronger than any argument about monotonicity: the
target of a post-hoc per-compound layer *is the error*, and an error computable from information
available at prediction time is information that belonged in the model. The R-squared column
(0.001 to 0.004) says the correction is not a function of the prediction at all — which is not a
licence for a cleverer layer but the opposite, since it rules out every function of `p` alone, not
just the monotone ones.

One error inside this file, caught by the file itself and kept because it is the same mistake in
miniature. The first version computed the ceiling as `p + isotonic(c | p)` and got 0.6880 —
**worse than the affine pair**, which is impossible for a ceiling. With the oracle correction
nearly uncorrelated with `p`, that isotonic fit is nearly flat and the "ceiling" degenerates to a
shift, while the pair also shrinks. The metric column contradicting the R-squared column is what
exposed it. Reporting a gate in variance units alone would have hidden it.

The procedure, for future use: a proposal that adds a per-compound correction on top of the
predictions is checked here **before** it is built, against magnitude rather than against class
membership. If the monotone budget is at the floor and the oracle correction tracks the residual,
the proposal is closed regardless of how it is motivated.

**129. CYP3A4 is not one dataset. It is two, and that explains most of the ways CYP3A4 has
behaved oddly in this log.** Found while checking a proposal to estimate the series parameter tau
from the screen instead of from the curves.

The check was arithmetic. The screen covers 4375 of 4905 rows, so if membership were independent
of structure, a similar pair would have both members in it 0.892^2 = 79.6 per cent of the time.
Measured, by similarity band:

    сходство     пар   обе в скрине   ожидалось
    0.60-0.85    601         14.0 %      79.6 %
    0.50-0.60   1734         34.5 %      79.6 %
    0.45-0.50   1321         58.7 %      79.6 %
    0.35-0.45   6850         79.9 %      79.6 %

At low similarity the assumption is exactly right; as similarity rises, screen membership
collapses to 14 per cent. The analog pairs in our training set are systematically **not** from the
diversity screen.

Following that gives the composition:

    группа                молекул   близкий сосед >= 0.587
    в скрининговой библиотеке  4375                    4.9 %
    вне её                      530                   68.1 %

    размечены вне скрина:  CYP1A2 0,  CYP2C9 0,  CYP2D6 0,  CYP3A4 530 из 2335

**All 530 compounds outside the screen carry exactly one label, and it is CYP3A4.** Their names
form a block (OCNT-049xxxx and up), their pIC50 distribution is different — median 4.45 against
4.20, interquartile range 1.19 against 1.64 — and their bands are narrower, 0.315 against 0.412.
Compound-ID medians order as diversity screen 2313565, campaign 2395534, test 2535312, with only
11 of 750 test compounds inside the campaign's range; if the identifiers are sequential these are
three successive batches and the test set is the newest.

So the training set is a diversity screen of 4375 singletons assayed on all four enzymes, glued to
a **CYP3A4-only analog campaign of 530 compounds** with no screening data, two thirds of which
have a close neighbour.

This retro-explains a long list of CYP3A4 anomalies that had been recorded separately:

- 2335 labels against 1285 to 1493 for the others — 530 of them are a different campaign;
- item 22's stratum landing 267 of 328 compounds on CYP3A4, and this session's anchor-split counts
  coming out 79 / 72 / 69 / **448**;
- item 127's selection propensity having the **opposite sign** on CYP3A4, +0.151 against -0.642
  and -0.821 — the campaign compounds were never selected by screening activity at all, they were
  designed;
- CYP3A4 being the enzyme where pooling hurts (item 111), where the dead zone hurts (item 122),
  where the aggregator flag does nothing (item 102), and where item 128's monotone budget is
  0.0002 against about 0.010 elsewhere.

Split by campaign, the model itself behaves differently: rank 0.798 on the 530 against 0.752 on
the 1805, and the mechanistic block is worth +0.0043 on the campaign against +0.0004 on the
screen. **Every CYP3A4 number in this log is an average over two populations** that differ in
geometry, in label distribution and in band width.

Two consequences worth acting on. The analog-series validation this session declared
unconstructible is constructible after all — on those 530, for CYP3A4 only. That is precisely the
stratum item 22 found empirically without knowing what it was, and knowing what it is turns a
curiosity into a designed control. And the tau proposal that started this: the screen yields 84
usable pairs per enzyme, not the ~475 an independence assumption predicts, which is still 2.8 to
5.6 times what the curves give on CYP1A2, CYP2C9 and CYP2D6 — and *fewer* than the 510 CYP3A4
already has, because those 510 are the campaign.

**130. The leaderboard identifies delta, but only as well as the test's label spread is pinned —
and that ties the delta programme to the series parameter.** Checked before the 24 September
submission, because it decides what that submission should be chosen for.

The claim under test: the intermediate leaderboard returns four numbers, ST-RAE is sensitive to a
global shift, so knowing our own prediction distribution the shift is recoverable per enzyme. That
would be worth a great deal — item 94 measured the model-choice systematic in delta at 0.11 to
0.51, wider than the entire bootstrap interval, and nothing local has narrowed it.

ST-RAE is a ratio, so there are two ways to be wrong about it: the numerator depends on prediction
quality as well as on the shift, and the denominator depends on the spread of the **test** labels,
which nobody has.

**Quality is not the problem.** Degrading the predictions until the rank falls by 0.01 — one and a
half noise floors — changes the score by as much as a shift of 0.041 / 0.035 / 0.010 / 0.061 does.
Against a current delta uncertainty near 0.5 that is a contamination of 2 to 12 per cent. On this
axis the claim is right and the margin is eightfold at worst.

**The denominator is the problem.** It enters multiplicatively, so an error in it converts to an
error in the recovered shift through the same slope:

    фермент     счёт    наклон   ошибка знаменателя 10 %   20 %
    CYP1A2    0.7987   -0.2473                     0.323  0.646
    CYP2C9    0.6252   -0.3989                     0.157  0.313
    CYP2D6    0.8796   -0.1685                     0.522  1.044
    CYP3A4    0.4880   -0.3026                     0.161  0.322

Between our own folds the denominator varies by 2.5 to 5.4 per cent, but the test is a different
population — analog series, a different batch, and on CYP3A4 a different campaign again (item 129)
— so ten per cent is the optimistic figure rather than the pessimistic one.

So delta is recoverable on CYP2C9 and CYP3A4 at about 0.16, marginally on CYP1A2 at 0.32, and
**not on CYP2D6**, where the slope is flattest and a ten per cent error costs 0.52 — the whole of
the current uncertainty. CYP2D6 is, as usual, the enzyme where it matters most.

The useful part is what pins the denominator. It is the mean absolute deviation of the test
labels, and the test is series of analogs around anchors that sit in our training set. Its spread
therefore decomposes into the spread of the anchor labels, which we have, and the within-series
spread, which is **tau** — the parameter the series layer needs and item 129 found 84 screening
pairs per enzyme to estimate.

**Measuring tau is what makes the leaderboard readable.** The two programmes are one piece of work,
not two, and that promotes tau above the selection layer in the order: without it, 24 September
returns four numbers we cannot convert into the quantity we most need.

**131. Removing the enzyme indicator does not merely cancel pooling's gain — it drops the model
0.076 of rank below training each enzyme separately.** `src/ablpool.py --arms "пул слепой"`, two
seeds. This is the arm items 110, 111 and 125 were converging on, and it was pre-registered in
item 111 before it ran.

The design differs from the pooled arm in exactly one thing: the enzyme indicator is zeroed rather
than removed, so the matrix keeps its width and the model keeps its 6525 rows. Only what it is
allowed to know changes.

    рука (2 сида)   MACRO ранг      1A2      2C9      2D6      3A4
    независимо          0.5627   0.4934   0.5940   0.3999   0.7636
    пул                 0.5778   0.5134   0.6027   0.4385   0.7566
    пул слепой          0.4868   0.4309   0.5507   0.2690   0.6967

Pooling with the indicator is worth +0.0151; pooling without it is worth **-0.0759**. The
indicator is therefore not a feature that helps, it is the condition under which pooling is not
harmful at all: forced to fit one function to four enzymes whose structure-activity relationships
diverge — CYP2D6 correlates with the others at -0.002 and with CYP2C9 at -0.120 (item 111) — the
model finds a compromise worse than any of the four separate fits.

That completes the elimination. Pooling does not work by borrowing neighbours (item 110), by
transferring shared function (item 111), or by sample size (item 125), and it does not work at all
without the indicator. What the indicator supplies is what pooling is.

**What this does not yet establish**, and the distinction is the last one standing. The indicator
can be doing either of two things, and this arm removes both at once:

  *level* — giving each enzyme its own offset. Medians run from 5.13 on CYP1A2 to 4.27 on CYP3A4,
  so a model that cannot tell them apart must average those together.

  *contrast* — letting the trees learn different dependencies, "this split matters for CYP2D6 and
  not for CYP3A4", which a per-enzyme model cannot represent at all.

The damage is at least suggestive: it is threefold uneven, CYP2D6 losing 0.1309 against CYP2C9's
0.0433, and a pure level effect ought to cost every enzyme about its own offset rather than
tracking how much it disagrees with the rest. But across four enzymes the rank correlation between
the damage and the mean between-enzyme correlation is only -0.4, not item 111's -1.000, so this is
a pointer and not a result.

The separating arm is `--arms "пул центрированный"`: labels centred per enzyme, indicator still
zeroed, the enzyme mean added back at prediction. That hands the model the level for free and
still denies it the contrast. If the gain returns, the indicator was carrying levels; if it does
not, contrast is what is left. It is queued.

**132. Pooling works by contrast, not by level. The mechanism of the log's largest effect is
settled.** `src/ablpool.py --arms "пул центрированный"`, seed 0. This closes a question open since
item 110.

Item 131 showed the enzyme indicator is what makes pooling work, but it removes two things at
once. The indicator can be supplying *level* — each enzyme's own offset, and the medians do run
from 5.13 on CYP1A2 to 4.27 on CYP3A4 — or *contrast*, the ability to learn different dependencies
per enzyme. The separating arm hands the level over for free and still denies the contrast: labels
centred per enzyme on the training folds, indicator still zeroed, the enzyme mean added back at
prediction.

    рука (сид 0)             ранг      1A2      2C9      2D6      3A4
    пул                    0.5792   0.5084   0.6029   0.4448   0.7607
    независимо             0.5651   0.4957   0.5972   0.4027   0.7646
    пул слепой             0.4885   0.4309   0.5507   0.2690   0.6967
    пул центрированный     0.4647   0.3860   0.5314   0.2339   0.7073

**Giving the model the level for free recovers none of the gain**, and lands slightly below the
blind arm. Level is not the mechanism; contrast is.

The ordering of the four arms is itself coherent and worth reading, because it says the centred
arm is the cleanest of the three negatives rather than an anomaly. In the blind arm the labels are
still uncentred, so the enzymes differ in level and features correlated with level let the model
partially infer which enzyme a row belongs to — a crutch it uses despite the zeroed indicator.
Centring removes that too, leaving the purest possible no-contrast model, and it is the worst.
0.4647 is therefore not "level did not help" but "this is what pooling is worth when contrast is
fully unavailable".

The elimination is now complete, and every step of it was a measurement rather than an argument:

    одалживание соседей       опровергнуто, пункт 110 (выигрыш крупнейший там, где нечего брать)
    перенос общей функции     опровергнуто, пункт 111 (связь перевёрнута, rho = -1.000)
    объём выборки             опровергнуто, пункт 125 (удвоение таблицы не даёт ничего)
    уровень фермента          опровергнуто, здесь
    контраст                  единственное, что осталось

So the largest effect the submission relies on has a mechanism at last: the indicator lets a tree
learn "this split matters for CYP2D6 and not for CYP3A4", which a per-enzyme model cannot represent
at all, and which is worth +0.0151 of rank against the per-enzyme baseline and +0.1145 against the
same pooled table with the contrast taken away.

Two consequences follow directly. Item 129's finding that only 41 compounds carry all four curves
means this contrast is learned almost entirely *across* molecules rather than within them — the
model never sees the same molecule on two enzymes — which is a strong constraint on any joint
latent built later. And it is the argument for the screening channel in `src/ablcontrast.py`: the
screen is the one place where the same molecule is measured on all four enzymes, 4376 times over,
so it is the only source of within-molecule contrast the dataset contains.

**133. The series layer would make things worse, and this is decidable without a single label.**
`verify/k34_series.py`. Half an hour, nothing built, and it closes a bet that item 130 had already
promoted in the queue.

The series layer shrinks a test compound's prediction toward the level of its analog series. Item
129 established that its benefit cannot be measured locally except on CYP3A4, so it looked like a
bet to be taken on faith. It is not: the *direction* is checkable on the test predictions we
already have, with no test labels at all.

Two quantities. **tau** is how much chemistry allows members of a series to differ — the standard
deviation of the label difference across training pairs in the 0.60-0.85 similarity band, divided
by root two. **The model's within-series spread** is the standard deviation of our own predictions
across the members of each of the 132 test series. If the model spreads a series more than
chemistry does, shrinking is justified and its coefficient follows; if it spreads it less, the
model is already over-smoothed and shrinking will compound the error.

    фермент     пар   tau (химия)   sd модели   отношение
    CYP1A2       30         0.659       0.231        0.35
    CYP2C9       15         1.006       0.288        0.29
    CYP2D6       23         0.398       0.162        0.41
    CYP3A4      510         0.669       0.390        0.58

**The model is already two to three times over-smoothed inside series on every enzyme.** The layer
is closed.

The obvious objection is that these predictions are post-affine-pair, and the pair shrinks by
lambda = 0.52 to 0.76, so some of the smoothing is the pair's rather than the model's. Dividing it
out gives pre-pair spreads of 0.398, 0.379, 0.312 and 0.527 against the same taus — ratios of 0.60,
0.38, 0.78 and 0.79. **The pair makes it worse but does not create it**; the boosting is already
under-spread within series before any post-processing touches it.

That is worth stating on its own, because it is a property of the model nobody had measured: on
near-neighbour molecules a gradient boosting interpolates, and the label does not. It is the same
fact item 107 met from the other side, where kNN's deficit *widened* as neighbours got closer, and
item 99's activity cliffs are the extreme case of it. Three observations, one phenomenon: **near
neighbours in this data disagree more than any smooth model of structure will predict.**

One correction to the proposal that prompted this, repeated because it changes the cost of the
remaining route. The suggestion was to estimate tau from the screen instead of the curves, gaining
"about 475 pairs per enzyme instead of 15 to 30". That figure assumes screening membership is
independent of pair membership, and item 129 measured that it is not: of the 601 pairs in the
0.60-0.85 band, **84** have both members in the screening library, not 475. The similar pairs are
overwhelmingly the CYP3A4 analog campaign, which has no screening data at all. Eighty-four is still
2.8 to 5.6 times what the curves give on CYP1A2, CYP2C9 and CYP2D6, so the route is worth having
— but it does not deliver a measured tau, only a better-estimated one, and the layer it was meant
to de-risk is closed on other grounds anyway.

**134. What the gap to the screening column is made of, and why most of it is not recoverable
from structure.** `verify/k36_ceiling.py`. This is the item that was built to be able to say
"stop", and on two enzymes it does.

The framing first, because it is the uncomfortable part. Across 130 items no intervention moved
macro ST-RAE by more than about 0.03 against a floor of 0.007. The distance from our model to a
single calibrated screening column, recomputed here on the same folds with the same affine pair:

    фермент   модель+пара   скрининг один   разрыв
    CYP1A2         0.8128          0.3106    0.502
    CYP2C9         0.6559          0.3154    0.340
    CYP2D6         0.9052          0.5992    0.306
    CYP3A4         0.5057          0.2168    0.289   [исправлено, см. пункт 139]
    МАКРО          0.7150          0.3605    0.354

**The largest untapped gap is CYP1A2, not CYP2D6**, and it is 1.6 times CYP2D6's. CYP2D6's
screening-alone score of 0.599 is the worst of the four *even with a direct measurement in hand*,
which says CYP2D6 is intrinsically hard rather than badly modelled. This repository's entire
mechanistic narrative — `frac_prot_74`, the basic amines, Glu216, the stratification of item 85 —
is built around the enzyme with the second smallest remaining headroom.

So instead of another arm: take the compounds carrying the gap — the top 15 per cent by how much
the screening column beats the model — and ask what they are, on axes split into chemical and
behavioural **before** looking.

    фермент  выигрыш   сходство    logP   колец   осн. N   полоса   метка
    CYP1A2     1.053     +0.031  -0.397  -0.200   +0.015   +0.982  -1.176
    CYP2C9     0.597     -0.004  -0.108  +0.017   +0.000   +0.253  -0.059
    CYP2D6     0.818     +0.000  +0.193  +0.177   +0.142   +0.178  +0.683
    CYP3A4     0.744     -0.009  +0.248  +0.018   +0.030   -0.083  +0.125

The answer is per-enzyme, which was pre-registered as readable, and it splits two ways.

**[ЧАСТИЧНО ОТОЗВАНО, см. пункт 135. Подписи остатка ниже верны как описание верхних 15 %;
причинное прочтение «разрыв на CYP1A2 в основном невосстановим» проверено срезом по порогу
прибора и не подтвердилось.]**

**On CYP1A2 the gap is behavioural, and that is a ceiling.** The residue's chemistry is flat —
similarity, size and basicity do not move — while its band is three times wider (1.418 against
0.436) and its label sits at 3.955 against 5.131, below the resolution floor item 105 measured.
These are compounds whose own pIC50 is an extrapolation from a curve the instrument never
resolved. The screen beats the model there because it **measures** "this compound is inactive" at
49.5 uM, while the model has to predict a number for a curve that does not exist. That is not
missing chemistry and no representation recovers it: **most of the 0.502 is a measurement beating
a prediction on labels that are themselves extrapolated.**

This also refutes, where it is testable, the proposal that prompted the analysis. CYP1A2's gap was
suggested to be shape selectivity — its site is a narrow planar slot and it takes flat aromatics.
The residue has **fewer** aromatic rings (-0.200) and **lower** logP (-0.397): the opposite of a
planar-aromatic chemotype. Whatever the shape story is worth elsewhere, it is not where this gap
lives.

**On CYP2D6 the gap is chemical, and it points at ground the project has already broken.** The
residue is more basic (+0.142 aliphatic nitrogens), more aromatic (+0.177 rings), more lipophilic
(+0.193) and **more potent** (+0.683) than the rest — potent basic amines, which is exactly the
pharmacophore the mechanistic block was built for. So the block is aimed correctly and is
under-reading: there is a chemotype here the model does not resolve, and it is the one the
repository already believes in. That is the most actionable single sentence in this file.

CYP2C9 and CYP3A4 show nothing on either group: their residues are indistinguishable from the
rest, so their gaps are diffuse rather than concentrated, and neither a chemotype nor a ceiling
argument applies to them from this evidence.

The practical reading. Macro-level effort against the macro gap is misdirected: half of the
largest component is a measurement advantage nobody can model away. Effort has an address, and it
is CYP2D6's potent amines.

**135. The causal half of item 134 is retracted: cutting the unresolved labels does not remove the
gap.** The test was proposed against item 134 within a day of its being written, and it is one
line.

Item 134 found that the compounds carrying most of the CYP1A2 gap have bands three times wider
than average and labels a full unit lower, and read that as: the gap is a direct measurement
beating a prediction on labels that are themselves extrapolations below the instrument's
resolution. That reading has a hard consequence -- restrict both columns to labels above
pC0 = 4.305 and the gap must collapse on CYP1A2 while barely moving on CYP2D6, whose share of
sub-floor labels is smallest.

    фермент   доля ниже pC0   разрыв весь   разрыв y > pC0   схлопнулся
    CYP1A2           20.5 %         0.502            0.431         14 %
    CYP2C9           31.9 %         0.340            0.398        -17 %
    CYP2D6           15.3 %         0.306            0.250         18 %
    CYP3A4           51.4 %         0.289            0.567        -96 %

**The prediction fails on both halves.** CYP1A2 gives up 14 per cent of its gap and CYP2D6 gives up
18 -- the enzyme that was supposed to move least moves most. And on CYP3A4 the gap nearly *doubles*
above the floor, so there the screen's advantage lives among the compounds whose curves resolved
perfectly.

The error, named plainly because it is one this file catches in other people's proposals and did
not catch in its own. **I took the signature of the top 15 per cent of the residue and attributed
the whole quantity to it.** The signature is real -- those compounds do have wide bands and weak
labels -- but removing them removes a seventh of the gap, so the bulk of CYP1A2's 0.502 lives in
the resolvable range and has no explanation yet.

The test is not tautological, and its design is what makes the failure informative: CYP3A4 has
*twice* the share of sub-floor labels (51.4 per cent) and the *smallest* gap, so "fraction of
unresolved labels" could never have explained the ordering by itself. That was checkable in
advance and would have weakened the claim before it was written.

What of item 134 survives:

  the gap decomposition, 0.502 / 0.340 / 0.306 / **0.289**, recomputed on the same folds
  (the CYP3A4 figure corrected in item 139 -- the two columns were scored on different masks);
  CYP1A2 having the largest untapped gap and CYP2D6 the second smallest, so the repository's
  mechanistic narrative is aimed at the enzyme with less headroom than the one it ignores;
  CYP2D6's screening-alone score of 0.599 being the worst of the four **with a measurement in
  hand**, which still says CYP2D6 is intrinsically hard rather than badly modelled;
  the residue signatures as descriptions of the extreme tail.

What does not survive: any claim that the CYP1A2 gap is largely unrecoverable, and with it the
"stop" conclusion the file was built to be able to reach. It reached it and the reaching was
wrong.

The refutation of the shape hypothesis weakens with it. The CYP1A2 residue does have fewer
aromatic rings and lower logP, which is the opposite of a planar-aromatic chemotype -- but since
that residue carries only a seventh of the gap, it is evidence about the tail and not about the
gap. The shape question on CYP1A2 is open again.

**136. The screening contrast is predictable and still adds nothing, so the idea closes rather
than its implementation.** `src/ablcontrast.py` measured the four predicted contrasts at -0.0015 of
rank against a base of 0.5651, and the *level* control at +0.0040 -- the control beating the
proposal, which was the pre-registered falsification. That left one ambiguity: a channel fed with
an unpredictable quantity is empty by construction, and a null would then close the implementation
rather than the idea.

Measured out of fold, structure onto screen:

    мишень            sd мишени   sd остатка      R2    ранг
    уровень               0.673        0.441   0.570   0.757
    контраст CYP1A2       0.604        0.498   0.320   0.580
    контраст CYP2C9       0.441        0.355   0.350   0.577
    контраст CYP2D6       0.834        0.686   0.324   0.506
    контраст CYP3A4       0.846        0.632   0.442   0.669

The contrast is recoverable from structure at an R-squared of 0.32 to 0.44 and a rank of 0.51 to
0.67. **It is not an empty channel.** So the model was handed a real, if noisy, estimate of which
enzyme a molecule prefers, and did nothing with it -- the idea is closed, not the wiring. The level
being better predicted (0.570) is consistent with its being the arm that moved.

The structural argument behind the proposal stands and is worth keeping separate from its failure:
the screen is still the only place in the dataset where the same molecule is measured on all four
enzymes, 4376 times against the curves' 41, and the contrast still carries 64.9 per cent of the
screen's variance. That the boosting cannot use it is a fact about the model, not about the
chemistry.

**137. The pairwise arm was restricted to distinguishable pairs, so what failed is the proposal
and not a weaker version of it.** Recorded because the distinction was raised. `pairwise_grad` in
`src/ablpairloss.py` builds its pair mask as `lo[i] > hi[j]` -- a pair enters the gradient only
when the bands are disjoint and i is strictly above j, and overlapping pairs receive exactly zero
weight. On CYP2C9, where 56.5 per cent of pairs are distinguishable, that discards nearly half of
them by construction. The smoke test's result -- pairwise 0.5991 against the same learner's squared
0.6103 -- is therefore about the restricted objective that was proposed.
**138. Feature subsampling does not free the mechanistic block, and there was nothing to free.**
Proposed as the mechanism behind a fifty-line booster beating the pinned one, and as the link to
the oldest open thread in this file: `frac_prot_74` takes part in no split on CYP2D6 while the
fingerprint reconstructs it at R-squared 0.653, so the basicity signal reaches the model
indirectly. If subsampling columns removes those bits from consideration at some nodes, the
explicit column should start being chosen.

Split counts over 200 trees, fold 0, at `max_features` 1.0 against 0.3:

    фермент     mf       FP     DESC     MECH   frac_prot_74
    CYP1A2     1.0   35.5 %   62.5 %    2.1 %        0.350 %
    CYP1A2     0.3   35.6 %   62.1 %    2.3 %        0.419 %
    CYP2C9     1.0   33.3 %   62.3 %    4.4 %        0.772 %
    CYP2C9     0.3   34.1 %   61.4 %    4.4 %        0.748 %
    CYP2D6     1.0   36.8 %   59.0 %    4.1 %        1.054 %
    CYP2D6     0.3   37.6 %   58.0 %    4.4 %        0.731 %
    CYP3A4     1.0   29.9 %   66.5 %    3.6 %        0.665 %
    CYP3A4     0.3   33.3 %   63.0 %    3.8 %        0.647 %

**On CYP2D6 the column is chosen thirty per cent less often with subsampling, not more**, and the
block shares are flat to within a percentage point everywhere. The hypothesis is refuted on the
enzyme it was made for.

The reason it could not have worked is in the same table and is worth more than the refutation.
**The dense blocks already dominate: DESC takes 58 to 66 per cent of splits while being 9.5 per
cent of the columns, and FP takes 30 to 38 per cent while being 89 per cent of them.** A
sixfold over-representation. There is no drowning of informative columns among two thousand sparse
bits, because the trees barely look at the bits to begin with -- which also explains why
`frac_prot_74` not being chosen is not a story about competition for attention.

Two consequences. The oldest thread stays open: `frac_prot_74` is unused on CYP2D6 for some reason
other than being crowded out, and item 134's residue says the compounds it should describe are
exactly where the model loses most. And if the own-booster's advantage survives four seeds, it is
**not** feature subsampling -- what remains is depth 5 against 31 leaves, 200 trees against 300, or
exact trees against histogram ones, and the `max_features` sweep now running separates the first
from the rest.


**139. The CYP3A4 gap is 0.289, not 0.269, and the cause is a mask that differed between the two
columns being compared.** `verify/k37_gap.py` scored the model on every labelled compound and the
screening column only on compounds that have a screening reading. Three of the four enzymes are
unaffected because their labelled sets are subsets of the screened set. CYP3A4 is not: item 129's
530-compound analog campaign carries curves and was never screened.

    фермент   модель на ВСЕХ   модель на общих   скрининг   разрыв кривой   разрыв верный
    CYP1A2            0.8128            0.8128     0.3106           0.502           0.502
    CYP2C9            0.6559            0.6559     0.3154           0.340           0.340
    CYP2D6            0.9052            0.9052     0.5992           0.306           0.306
    CYP3A4            0.4860            0.5057     0.2168           0.269           0.289

The correction makes CYP3A4's gap larger, not smaller, so nothing that depended on it being the
smallest changes -- but the number is wrong wherever it was quoted and is fixed above.

**140. The own booster, four seeds: the reference learner is neither the fastest nor the best, and
the mechanism is column subsampling.**

    рука                        MACRO пара   MACRO rho     1A2     2C9     2D6     3A4    время
    HistGB (эталон, сид 0)          0.7150      0.5651   0.496   0.610   0.408   0.760   ~2400 с
    свой бустинг, квадрат mf1.0     0.7278      0.5522   0.480   0.583   0.389   0.755   456-716
    свой бустинг, квадрат mf0.3     0.7204      0.5615   0.486   0.601   0.405   0.754   119-179
    свой бустинг, квадрат mf0.1     0.7171      0.5644   0.493   0.602   0.405   0.757    47-67
    свой бустинг, попарно mf0.3     0.7172      0.5692   0.519   0.597   0.424   0.738   200-262

Three readings, in order of how well they hold.

**Column subsampling is worth +0.0123 of rank and is monotone in it** (0.5522, 0.5615, 0.5644 as
max_features goes 1.0, 0.3, 0.1). The sign is the same on all four seeds and the size is nearly
twice the 0.007 floor. The pinned scikit-learn's HistGradientBoostingRegressor **has no
max_features at all** -- checked -- so this is not a setting the journal could have swept; changing
the learner was the only route to it, which is the one good reason to have written one.

**The pairwise objective is worth +0.0077 over the same learner's squared loss**, sign consistent
across all four seeds (+0.0046, +0.0090, +0.0049, +0.0124), size at the floor. It is not a uniform
improvement and the trade is the interesting part: **+0.033 on CYP1A2 and +0.018 on CYP2D6,
-0.016 on CYP3A4 and -0.004 on CYP2C9.**

**HistGB's slowness is an algorithm-to-data mismatch, not a misconfiguration**, and this was worth
checking before four seeds rather than after. Measured: 300 iterations actually run, `n_iter_`=300,
`do_early_stopping_`=False, **2378 ms per tree** on 1028 rows by 2295 columns. Against 15 ms for an
exact tree at mf0.1 and about 40 ms at mf0.3 -- 161-fold and 60-fold. Nothing pathological is
happening. A histogram booster builds a histogram over every feature at every node at a cost of
O(features x bins) that barely depends on the number of rows; it is built to win when rows are
many, and 1028 rows against 2295 columns is the shape where the overhead is the whole cost and
buys nothing. **The consequence is throughput, not accuracy**: an arm that took an hour takes two
minutes, and that is a different working regime for the nine weeks that remain.

**141. Rank as this file measures it counts pairs the metric will never pay for, and CYP3A4 is
where that bites.** Item 140's pairwise arm loses CYP3A4 rank while winning overall, which is the
wrong shape for a loss that only ever discards pairs the bands cannot order -- so the suspicion
falls on the criterion.

After the affine pair, ST-RAE pays for an order only when the two confidence bands are disjoint:
if they overlap, no prediction can be penalised for putting the pair either way round. Spearman
against the point label pays for **every** pair. CYP3A4 has 51.4 per cent of its labels below the
instrument's resolution floor -- the highest of the four -- so a large share of its rho is earned
on pairs the competition cannot see, and an objective that deliberately drops those pairs must
look worse under rho while being no worse under the metric. Consistent with that, the pairwise
arm's macro *метрика* (0.7172) is level with the best squared arm (0.7171) while its rho is 0.0048
higher and its CYP3A4 rho 0.0165 lower.

Item 80 made rank the criterion because raw ST-RAE was untrustworthy, and that was right. The
qualifier it needs is that rho over-counts on wide-banded enzymes. **[ИЗМЕРЕНО И ОПРОВЕРГНУТО,
см. пункт 146 --- квалификатор не нужен, критерий в порядке.]** `verify/k42_visiblerank.py`
recomputes the file's headline comparisons under a Kendall tau restricted to distinguishable
pairs -- the same criterion with the invisible pairs removed -- and is queued. Pre-registered: if
CYP3A4's loss shrinks under tau, the loss was the criterion and every per-enzyme comparison here
involving CYP3A4 carries the same bias; if it survives, the trade is real and the arm is a
per-enzyme choice.

**142. The pooled arms sit on both sides of a threshold inside the learner, so item 125 compared
two different learners.** The pinned scikit-learn defaults `early_stopping='auto'`, and 'auto'
means on above 10000 samples and off below.

    рука                всего строк   обучающих (4/5)   ранняя остановка
    независимо (2D6)           1493              1194   нет
    пул                        6525              5220   нет
    пул+TDI                   13063             10450   ДА
    пул+скрининг              18030             14424   ДА

`пул` trained on all its rows for all 300 iterations. `пул+TDI` held out ten per cent as a
validation set and stopped when that set stopped improving. The margin is 450 rows -- 4.5 per cent
over the line -- which is why it went unseen. This does not show item 125's conclusion is wrong;
it shows the measurement did not test what it was written to test, and that matters because
`пул+TDI` is the arm that closed "more supervision helps", a question items 143 and 144 both
reopen with much larger tables that would cross the same line. `verify/k41_earlystop.py` reports
`n_iter_` and both regimes on one fold and is queued; every new arm sets `early_stopping=False`
explicitly.

**143. Item 136 closed the screening idea on a test that had no branch on which it could have
succeeded.** `src/ablcontrast.py` handed the model the level and contrast **predicted from
structure** as extra columns and measured -0.0015 of rank. Item 136 pre-registered that an
unpredictable channel would be empty by construction, measured predictability at R-squared 0.32 to
0.44, found it non-zero, and concluded the idea rather than the wiring was closed. But both
branches lead to a null:

    непредсказуемо из структуры -> канал есть шум       -> ноль
    предсказуемо из структуры   -> канал есть функция X -> избыточен -> ноль

so the null carried no information about the screen. What was never given to the model is the
screen's **measurements**, and the preconditions for doing so are the strongest in this file:

    фермент   кривых   скрининг без кривой для этого фермента
    CYP1A2      1412                                     2963
    CYP2C9      1285                                     3090
    CYP2D6      1493                                     2882
    CYP3A4      2335                                     2570
    итого       6525                                    11505

1.76 times the labelled table, on compounds that table never mentions for that enzyme. And the
second number is the one that makes this the natural next arm rather than one more idea: item 132
settled that pooling works by **contrast**, and the curve table carries all four enzymes for **41**
molecules while the screen carries all four for **4375** -- a hundred and seven fold. The place
where contrast is actually measured has never been in the training table.

`src/ablaux.py` pools curve rows and screen rows with an enzyme indicator and a source column,
calibrating log2fc to pIC50 by an isotonic fitted **on training folds only**, and its decisive arm
is a permutation control: the same rows, the same target marginal, the molecule-to-measurement
link destroyed. A gain that does not survive the permutation is regularisation from extra rows and
not the screen. Pre-registered per enzyme: the gain should fall as 2C9 (2.4x), 1A2 (2.1), 2D6
(1.9), 3A4 (1.1); a gain concentrated on 3A4 refutes the reading. Queued.

**144. The NCGC panel is merged, does not touch the test set, and cannot be calibrated -- so it
will not be.** `src/ncgcmerge.py` folded the six fetched assays into 54177 rows over 13126
structures, `src/ncgcfeats.py` built their features through `feats.build` with the descriptor and
mechanistic columns pinned by name.

**Zero of the 750 blinded test molecules appear in the panel**, so the question of whether public
prior measurements on test rows may be used does not arise. It was measured first and printed
first precisely because it is not a question a script may settle.

The calibration precondition failed, and this is the finding rather than a setback. Molecules
carrying both a fitted NCGC AC50 and one of our curves:

    фермент   общих   rho    сдвиг NCGC-наш   sd разности
    CYP1A2       32  0.72            +0.762         0.577
    CYP2C9       11  0.40            +0.448         0.710
    CYP2D6       48  0.74            +0.442         0.484
    CYP3A4       21  0.49            +0.868         0.843

The offset is **larger** than the ChEMBL offset items 60 to 63 rejected, and eleven molecules at
sd 0.71 give it a standard error of 0.21 against a quantity of 0.45. So the arm that subtracts a
fitted shift is not built. What is built is item 63's actual lesson -- a **source indicator**,
which lets the model learn the offset from all 54177 rows rather than from eleven, and as an
interaction with chemistry rather than as a constant.

Two things the merge had to get right and did. Panel rows are mapped to our folds by canonical
SMILES and dropped when their molecule is held out, because 183 to 319 panel structures are also
ours and their measurement would otherwise be a measurement of a test row. And censored rows --
7258 on CYP3A4 alone -- are kept as a marked arm rather than discarded, because dropping a "not
active up to the top concentration" is rebuilding by hand the selection that killed the ChEMBL
merge. `src/ablncgc.py` is queued, with a permutation control and a CYP2C19 arm that asks whether
a fifth isoform we will never be graded on still teaches the shared trunk.

**145. What the intermediate leaderboard on 24-25 September is for, decided in advance.** This
file has twice proved that validation in the test regime cannot be built here: reweighting is
forbidden by the chi-square of 2.838 between the test and out-of-fold similarity distributions,
which caps the effective sample at 26.1 per cent (item 123), and construction is forbidden by
composition, 93.6 per cent Butina singletons (item 129). The leaderboard is therefore **the only
sample from the test distribution that will exist before the close**, and spending it purely on
"submit the best guess" wastes the one measurement money cannot buy.

The largest open question it could answer is whether pooling's contrast gain survives the shift.
Pooling is the biggest effect the submission relies on, its mechanism is now identified (item 132)
but rests on the enzyme-conditioning structure being stable between distributions, and nothing
internal can test that. So:

  if more than one submission is allowed, one per-enzyme and one pooled. The **difference** between
  their scores is worth more than either rank, because it is the only external measurement of the
  assumption the submission is built on. **[Item 147 strengthens this: the difference is not merely
  worth more, it is the only quantity the leaderboard measures cleanly.]**

  if only one, submit the best and pre-register now what leaderboard score would falsify the
  out-of-fold estimate, so that the answer cannot be re-read afterwards.

Written before the date rather than after it, which is the whole point.

**146. The rank criterion is safe: restricting it to the pairs the metric can see changes nothing.**
`verify/k42_visiblerank.py`. Item 141 suspected that Spearman over-counts on wide-banded enzymes
and that CYP3A4's loss under the pairwise arm was an artefact of the criterion. Both halves are
now measured and both are wrong.

The statistic is a Kendall tau over ordered pairs with disjoint bands -- the same criterion as rho
with the pairs ST-RAE cannot pay for removed. Visible-pair shares are 73.6, 56.6, 73.5 and 72.5
per cent.

    рука                          MACRO rho   MACRO tau   разность
    эталон FP+DESC+MECH              0.5651      0.5209     0.0442
    свой бустинг, квадрат mf1.0      0.5522      0.5074     0.0448
    свой бустинг, квадрат mf0.3      0.5615      0.5164     0.0451
    свой бустинг, квадрат mf0.1      0.5644      0.5196     0.0448
    свой бустинг, попарно mf0.3      0.5692      0.5226     0.0466
    L1 по метке                      0.5619      0.5167     0.0452
    мёртвая зона x1                  0.5923      0.5453     0.0470
    пул                              0.5792      0.5327     0.0465

**The offset is 0.044 to 0.047 for every arm.** At this resolution the two criteria are an affine
reparametrisation of each other: they order the arms identically, and every difference the file
has recorded under rho survives under tau at nine tenths of its size. Item 80's rule needs no
qualifier.

The specific prediction fails too. CYP3A4 under the pairwise arm against the same learner's
squared loss: **-0.0164 of rho and -0.0175 of tau.** The loss does not shrink when the invisible
pairs are removed, it grows slightly. So the pairwise objective really does give up CYP3A4 order
that the metric would have paid for, the trade of item 140 is real, and that arm is a per-enzyme
choice rather than an improvement.

One thing the table shows that was not being asked. Per enzyme the *level* gap between rho and tau
is large and uneven -- CYP1A2 -0.062, CYP2C9 **+0.009**, CYP2D6 -0.047, CYP3A4 -0.077 -- and CYP2C9
is the only enzyme where the model does better on the pairs that count than on all pairs, while
having the fewest visible pairs of the four (56.6 per cent). Levels are not comparable across
enzymes under either statistic; differences between arms are. That was already the rule and now
has a measurement behind it.

**The incidental result is the one worth acting on.** Read down the macro column: the dead zone at
+0.0272 of rho and +0.0244 of tau is **nearly twice pooling's** +0.0141 and +0.0118, and it wins on
three enzymes out of four (+0.039 CYP1A2, +0.034 CYP2C9, +0.044 CYP2D6, -0.007 CYP3A4). Pooling is
the effect this repository's submission is built around and the one whose mechanism took four
eliminations to find. The largest measured effect in the file is now something else.

**147. A rigorous bound on the leaderboard score is vacuous, and the reason says what to submit.**
`verify/k43_lbpredict.py`, written before 24 September so the threshold cannot be chosen once the
score is known.

Three sources of movement, separated so they can be argued about apart:

    фермент   ST-RAE вне фолда   выборка на 750   числитель под chi2   знаменатель под chi2
    CYP1A2              0.8128      0.734-0.908          0.000-2.387    -0.395 .. 1.379  ВЫРОЖД
    CYP2C9              0.6559      0.555-0.774          0.000-2.289    -0.368 .. 0.996  ВЫРОЖД
    CYP2D6              0.9052      0.847-0.976          0.000-2.654    -0.400 .. 1.238  ВЫРОЖД
    CYP3A4              0.4860      0.420-0.560          0.000-1.542    -0.313 .. 1.360  ВЫРОЖД
    МАКРО               0.7150      0.639-0.804          0.000-2.218

**The chi-square bound is vacuous on all four enzymes, and it is the metric's shape that does it.**
ST-RAE's denominator is a constant predictor at the mean put through the same soft threshold, so
for a large share of compounds it is exactly zero -- the mean already lands inside the band. Its
standard deviation therefore exceeds its mean, and the Cauchy-Schwarz bound at radius 2.838 admits
a reweighting that drives the denominator to zero, leaving the ratio unbounded above. A bound that
excludes nothing cannot falsify anything, however correctly it was derived. This is the third
independent route to the conclusion items 123 and 129 reached: the test regime is not checkable
from here.

**The sampling interval is the number nobody had.** With no shift at all, a 750-molecule draw moves
macro ST-RAE over 0.639 to 0.804 -- **plus or minus 0.08, five times the seed spread of 0.016 and
eleven times the 0.007 floor.** Every effect this file has ever measured is smaller than the noise
in a single leaderboard score. That is not an argument against the effects; it is an argument about
what a leaderboard number can be read as.

**What follows for the submission, and it is sharper than item 145 had it.** The denominator is a
property of the test set, not of the submission, so two submissions scored on the *same* test set
share it exactly and it cancels in their difference. The difference between two submissions is
therefore the only quantity the leaderboard measures cleanly -- free of the unknown denominator,
free of the sampling draw, and paired. So if more than one submission is allowed, the second one is
not a spare guess: it is the measurement. Per-enzyme against pooled tests the assumption the whole
submission rests on, and nothing internal can test it.

If only one is allowed, the pre-registered statement is the weaker one that remains honest: a score
outside 0.639 to 0.804 (widened by the seed spread) falsifies **the named assumption** -- that the
test's label spread and band widths resemble the training set's -- and not the model. A score inside
confirms nothing at all.

One methodological note, because the first version of this script got it wrong and the error was
the quiet kind. `cluster_ids` returns a tuple `(cid, n_clusters)`, and wrapping it in a
`try`/`except` made `np.asarray` raise and fall back silently to bootstrapping over molecules
rather than clusters. The interval that produced was 0.641 to 0.799 -- narrower, plausible, and
wrong for the right-looking reason. The fallback is now gone and the cluster count is asserted
against the split's own: 4703. That the correction moved the interval so little is itself the
composition result of item 129 showing through -- 93.6 per cent of clusters are singletons, so
clusters and molecules are nearly the same unit here.

**148. The dead zone reproduces on all four seeds at +0.033 of rank, and its pre-registered
mechanism is refuted by its own numbers.** `src/abldead.py`, seeds 1, 2 and 3 against seed 0's
+0.0272:

    рука               MACRO пара   MACRO rho     1A2     2C9     2D6     3A4
    L2 по метке            0.7157      0.5623   0.4956  0.5899  0.4038  0.7602
    L1 по метке            0.7137      0.5655   0.5073  0.5851  0.4093  0.7603
    мёртвая зона x1        0.6872      0.5985   0.5364  0.6330  0.4590  0.7657

**+0.0362 over the squared loss on seeds 1-3, +0.0272 on seed 0**, four seeds agreeing in sign and
roughly in size, against a floor of 0.007. It is the largest reproducible single-model effect in
this file, larger than pooling's +0.0141.

And the mechanism story is wrong. `abldead` pre-registered: "the gain must be uneven across
enzymes: the share of wide bands runs from 9.9 per cent on CYP2D6 to 29.4 per cent on CYP3A4, so
CYP3A4 must gain most." Measured:

    фермент   доля полос шире 1.0   прирост ранга
    CYP2D6                  9.9 %         +0.0552   <- меньше всех полос, больше всех прирост
    CYP1A2                 14.6 %         +0.0407
    CYP2C9                 17.0 %         +0.0431
    CYP3A4                 29.4 %         +0.0056   <- больше всех полос, почти ноль

**The ordering is inverted, not merely different.** The explanation is in `abldead`'s own docstring
two paragraphs above the prediction, unused: on CYP3A4 **73.9 per cent of the predictions on
wide-banded rows already land inside the band**. Where that is true, `clip(p, lo, hi)` equals `p`,
the majorised target is the current prediction, the gradient is zero and the pass does nothing. A
wide band means the error is already free -- so there is nothing there to win. The gain comes from
**narrow** bands, where the model sits outside and reprojection actually moves the target.

The intervention is unaffected and the reading of it is inverted: the dead zone is not a way to
stop wasting effort on wide bands, it is a way to stop paying for the last tenth of a pIC50 on
narrow ones.

**149. The dead zone applied to every ensemble member gives 0.6203 macro rank -- the best number in
this file -- and applying it to the boosters alone makes things worse.** `src/abldzens.py`:

    рука                          MACRO пара   MACRO rho     1A2     2C9     2D6     3A4
    базовый ансамбль                  0.6793      0.6029  0.5341  0.6407  0.4561  0.7807
    мёртвая зона в бустингах          0.6799      0.6001  0.5310  0.6439  0.4474  0.7781
    мёртвая зона везде                0.6599      0.6203  0.5463  0.6661  0.4710  0.7978

The arm was split because the majorise-minimise reduction is exact only under an absolute loss: on
the boosters it is L1 against the reprojected target, while the Gaussian process and the ridge get
a dead zone with a *square* outside it, a different estimator. The pre-registered question was
whether that approximation survives. It does, and by a margin:

    пул          0.5792 -> 0.5916   (+0.0124)
    GP           0.5534 -> 0.6044   (+0.0511)
    гребневая    0.5595 -> 0.5778   (+0.0183)

**The member with the inexact reduction gains four times what the exact one gains.** So the split
was worth making and the answer is the opposite of the caution that motivated it.

Members gain +0.027 on average while the ensemble gains +0.0174, which is the diversity cost the
script warned about, measured: reprojecting every member onto the same bands does make them err
more alike, and it costs about a third of the per-member gain. It is still the largest ensemble
number here. Seeds 1-3 are queued.

**150. The fingerprint reduces to about fifty substructures, and that is the first chemical
sentence this file can say.** `verify/k40_topk.py`. Item 96 measured the whole fingerprint block at
+0.032 of rank; item 138 measured that the trees spend only 30 to 38 per cent of their splits on 89
per cent of the columns. Ranking bits by how often they are chosen **on the training folds only**
and keeping the top k:

    k        MACRO пара   MACRO rho   доля эффекта фингерпринта
    без FP       0.7411      0.5330                          --
    top-20       0.7290      0.5554                      69.8 %
    top-50       0.7228      0.5589                      80.7 %
    top-100      0.7182      0.5624                      91.6 %
    top-200      0.7176      0.5633                      94.4 %
    все 2048     0.7150      0.5651                     100.0 %

**Fifty bits out of 2048 carry four fifths of the effect, and twenty carry seven tenths.** The
fingerprint is not working in bulk.

**[ДВЕ ПОПРАВКИ ИЗ ПУНКТА 156, обе сужают заявление.** Во-первых, `k40_topk` отбирает
пятьдесят битов НА ФЕРМЕНТ и внутри каждого фолда, а не пятьдесят на всю задачу; списки
ферментов пересекаются на 7-11 битов, общих у всех четырёх --- три, так что объединение
около ста семидесяти. Во-вторых, «подструктуры» --- преувеличение: тринадцать из полусотни
имеют радиус 0, то есть это типы атомов на трёх с половиной тысячах молекул, а не алерты.
Расшифровка в пункте 156.**]

**151. Pooling reverses on a depth-5 tree with 30 per cent column subsampling, which is a warning
about the submission and not yet a result.** `src/ablaux.py` ran a pooled arm on the own booster of
item 140 and got **0.5420 against the per-enzyme 0.5687** -- pooling *losing* 0.0267, where on
HistGradientBoostingRegressor it wins 0.0141.

The gate was pre-registered ("if pooling does not reproduce on this learner there is nothing
further to read") and it failed, so the rest of that run is read below only in its internally
controlled comparisons.

There is a mechanical hypothesis and it is sharp. The enzyme indicator is **one column out of
2300**. At `max_features=0.3` it is absent from seventy per cent of split decisions. A pooled model
whose entire mechanism is conditioning on that indicator (item 132) is therefore starved of it,
while a per-enzyme model has no indicator to lose -- and item 140 measured that the same
subsampling is worth **+0.0123** to per-enzyme models. The two findings interact: column
subsampling helps every model that does not need a specific column and cripples the one that does.
Queued at `max_features=1.0`, where the indicator is always available. If pooling returns there,
the reversal is the learner; if it does not, pooling's gain is HistGB-specific and the submission
rests on something that survives one learner out of two.

**152. The screening measurements carry real information -- the permutation control says so by
+0.036 -- but on this learner it only buys back what pooling lost.** The three pooled arms of
`src/ablaux.py` share a learner and are comparable to each other whatever item 151 says:

    рука                          MACRO rho     1A2     2C9     2D6     3A4
    независимо                       0.5687  0.4965  0.6103  0.4078  0.7604
    пул                              0.5420  0.4754  0.5655  0.3901  0.7371
    пул+скрининг                     0.5727  0.5049  0.6207  0.4047  0.7605
    пул+скрининг, перемешанный       0.5365  0.4895  0.5500  0.3859  0.7204
    пул+скрининг w0.3                0.5688  0.4973  0.6178  0.4045  0.7557

**The decisive comparison is +0.0362**: the same 11505 rows, the same target marginal, the
molecule-to-measurement link broken, and the score falls by that much. The regularisation
explanation is dead -- extra rows with broken links are *worse* than no extra rows at all
(0.5365 against 0.5420). What the screen adds is measurement.

The per-enzyme ordering was pre-registered as 2C9 (2.4x more screening than curves) > 1A2 (2.1) >
2D6 (1.9) > 3A4 (1.1), with a gain concentrated on 3A4 as the falsifier. Measured gain over the
pooled arm: **2C9 +0.055, 1A2 +0.030, 3A4 +0.023, 2D6 +0.015.** The first two land as predicted,
the last two swap while both stay small, and the falsifier does not fire.

What cannot be claimed. Against the *per-enzyme* reference the screen is worth **+0.0040**, below
the floor -- it recovers pooling's loss and a hair. Whether it adds anything to a model that is
not paying pooling's penalty is exactly what the queued `независимо+скрининг` arm asks, and it is
the arm that should have been in the first design.

**153. The NCGC panel hurts as a source of rows, its measurements are real, and its censored rows
are the half worth keeping.** `src/ablncgc.py`, 18711 fitted-AC50 rows against our 6525:

    рука                          MACRO rho     1A2     2C9     2D6     3A4
    база (пулированные кривые)       0.5462  0.4811  0.5669  0.3933  0.7435
    +NCGC активные                   0.5286  0.4905  0.5550  0.3447  0.7243
    +NCGC активные, перемешанный     0.5141  0.4799  0.5364  0.3209  0.7193
    +NCGC активные+2C19              0.5204  0.4856  0.5412  0.3324  0.7225
    +NCGC всё, с цензурой            0.5365  0.4854  0.5781  0.3511  0.7315

Every arm is below the base, so as built the panel is a cost. Three things are nonetheless
measured rather than guessed.

**The measurements are real**: +0.0145 over the permutation control, the same test that vindicated
the screen. The panel knows something; it is the merge that is wrong, not the data.

**The censored rows help, by +0.0079 over the fitted-only arm.** Keeping "not active up to the top
concentration" as a row was the right call, and dropping them -- the obvious tidy-up -- would have
rebuilt ChEMBL's selection bias by hand and cost more than it saved.

**CYP2C19 teaches nothing transferable**: adding the fifth isoform is 0.0082 *worse* than not. That
is the only testable form of a proteochemometric enzyme coordinate available here, and it is
negative.

The design error is visible in hindsight and is mine: at weight 1.0 the panel is 74 per cent of the
training rows, so a foreign protocol with an offset of +0.44 to +0.87 dominates the objective and
the source indicator has to undo it from a minority position. The screen arm above had a w0.3
control and this did not.

**154. Four ideas close, and one of them saves three days.**

`src/ablcoord.py` -- the enzyme as a two-column MDS coordinate of the measured screening
correlation matrix, against the four-column one-hot: **0.5799 against 0.5792**, and both together
0.5830. Seven ten-thousandths, an eighth of the floor. Geometry is not what the indicator carries.
The pre-registered consequence was explicit: the sequence-based variant, three to four days of
assembling external isoforms, **should not be started**. It will not be.

`src/ablpre.py` with MoLFormer-XL -- база 0.5651, embedding alone 0.4370, base + embedding 0.5386,
base + embedding through a 64-component PCA fitted inside each fold **0.5637**. The width story
from ChemBERTa reproduces: the full embedding hurts and PCA-64 restores the baseline, but restores
is all it does. Pretrained SMILES representations add nothing to this matrix, on two checkpoints
now, and that is the first clean answer to a question the file had left open three times.

`src/abllipo.py` -- logD at pH 7.4 as a feature is +0.0051, below the floor; LipE as a target is
-0.0174. The basic centre exists in 81.8 per cent of molecules and the median logD shift among
bases is **0.00**, which is why: at pH 7.4 our bases are mostly not protonated enough for logD to
differ from logP, so the new column is nearly the old one.

`src/ablrobust.py` -- chi-square DRO over similarity strata costs rank as soon as it is switched
on, and keeps costing it: **0.5651, 0.5574, 0.5456 at radii 0, 0.25 and 0.5**, and it costs rank
under the test-weighted criterion too (0.5963, 0.5866, 0.5599). The run was stopped after three of
its five radii, with 1.0 and the measured 2.155 unrun: the decline is monotone in both criteria and
a larger ball can only widen the worst case it optimises against, so the remaining radii could
confirm the conclusion but not reverse it, and the machine was needed for the NCGC re-weighting.
Recorded as a truncation rather than left to be inferred from a short log. The stratum shares are the striking part: out of five quantile strata the
test puts **77.1 per cent of its mass in the highest-similarity one** against our uniform 20 per
cent, and the radius measured for this partition is 2.155. Robustness to the whole ball is the
wrong ask when the shift is that concentrated in one direction.

**155. Under the anchor protocol the mechanistic block's gain is a third of what the cluster
protocol reports.** `verify/k32_anchor.py` on item 129's 530-compound CYP3A4 campaign, the one
place a test-like split is constructible. The construction check passes: median similarity of the
held-out set to training rises to 0.61-0.62 across three rotations.

    выигрыш мех. блока    кластерный    якорный
    ST-RAE                   -0.0262    -0.0172
    ранг                     +0.0056    +0.0018

Both are small and the rank figures are inside the floor, so this is not a refutation. But it is
the first time the file can *see* rather than infer that a number depends on the protocol, and the
direction is the uncomfortable one: the more test-like the split, the smaller the mechanistic
block's contribution.

**156. The fifty bits decoded: about five nameable fragments, and they are the textbook CYP
pharmacophore.** `verify/k44_bits.py`. Item 150 promised a chemical sentence; this is the sentence,
and getting to it required withdrawing two thirds of the promise first.

**What the top fifty actually consist of.** Morgan environments are hashed into 2048 buckets, so
before any chemistry two things had to be counted: how many distinct environments land on each bit,
and what radius they are.

    радиус   битов   медиана молекул   чистых (главное окружение >= 0.8)
    0           13              3538                                 10
    1           26               618                                 17
    2           11               147                                  4

**Thirteen of the fifty are radius zero** -- a single atom with its invariants, carried by three and
a half thousand molecules out of 4905. On a *count* fingerprint that is an atom census, not an
alert, and it duplicates what MolWt, NumAromaticRings and the heteroatom counts already give the
descriptor block. Twenty-six more are radius one, a single atom and its bonded neighbours: `c(c)c`
is the interior of a benzene ring, present in 4430 molecules. **Only eleven are radius-two
fragments and only four or five of those are clean enough to name.** The median bit collects eleven
distinct environments.

Two further limits, measured rather than assumed. The ranking agrees with the fold-wise ranking on
**36 of 50** bits, so the list itself is about seventy per cent stable. And the enzymes barely
share it: pairwise overlap is 7 to 11 bits and **three** are common to all four, which is what
forces the correction to item 150 above.

**The chemistry, controlled.** A raw difference in mean pIC50 between carriers and non-carriers
would confound the fragment with the size and lipophilicity of the molecules that carry it, so each
is also measured on the residual after linearly removing MolWt, MolLogP, NumAromaticRings, TPSA,
NumHAcceptors and FractionCSP3. Raw / residual:

    фрагмент                       молекул      1A2          2C9          2D6          3A4
    c(cn)nc   пиримидин, C2 между
              двумя кольцевыми N       150   +.78/+.67    +.57/+.65    +.66/+.68    +.80/+.66
    c(nc)n(-c)c  N-арил-азол            74   +.64/+.61    +.50/+.62    +.77/+.84    +.88/+.97
    n(c)c     пиридиновый N            624   +.57/+.42    +.42/+.38    +.33/+.37    +.75/+.72
    c(cc)nc   углерод азина            222   +.32/+.13    +.52/+.46    +.19/+.23    +.60/+.62
    c(OC)(cc)cc  метоксиарил           147   +.09/+.00    +.19/+.19    +.51/+.53    +.50/+.34
    C(CC)N(C)C  трет. алифат. амин     328   -.85/-.58    +.10/+.12    +.00/-.12    -.09/-.06
    C(C)N     алифатический амин      2104   -.51/-.23    -.23/-.05    +.08/-.09    -.20/-.09
    C(C)O     алифатический гидроксил  588   -.35/-.11    -.26/-.08    +.16/+.05    -.30/-.20
    c(cc)c(-c)c  биарил                344   +.43/+.21    +.24/+.02    -.15/-.03    +.23/-.15

**Every fragment that survives the control is a ring nitrogen with a lone pair available to the
haem iron.** Pyrimidine, N-aryl azole, pyridine, azine -- the classical coordinating pharmacophore,
and the azole's effect *grows* under the control on CYP2D6 and CYP3A4 (+0.77 to +0.84, +0.88 to
+0.97). The most used bit in the whole fingerprint, 126 splits against the runner-up's 37, is the
pyridine nitrogen. The single largest effect belongs to the pyrimidine C2 at +0.66 to +0.68 on all
four enzymes at once.

**Everything that collapses under the control was a size marker.** The biaryl bit goes from +0.43
to +0.21 on CYP1A2 and from +0.23 to -0.15 on CYP3A4 -- it was measuring how large the molecule is.
The aliphatic hydroxyl and the generic aliphatic amine lose half to two thirds. Had the raw column
been read alone, three of nine fragments would have been over-read as chemistry.

**One enzyme-specific statement survives, and it is CYP1A2's.** The tertiary aliphatic amine costs
**-0.58 of pIC50 after the control**, on CYP1A2 alone -- the other three enzymes are flat. That is
the sp3, basic, non-planar character being rejected by a narrow planar site, which is the shape
argument item 134 tried to test on the residue and could not settle. It is settled here in the
opposite place: not among the compounds where the screen beats the model, but across the labels
themselves.

**What this is worth.** Not a feature: the model already has these bits and is already using them,
so nothing here raises a score. What it is worth is that a 2048-column block whose contribution was
a number is now a mechanism -- the fingerprint's job is to carry haem-coordinating ring nitrogen
and sp3 basicity, two things no descriptor in the 217-column block encodes directly. That also
explains item 138's shape: the trees consult the fingerprint rarely because only a handful of its
columns say anything, and they cannot drop it because nothing else says that.

**157. The NCGC panel's harm was its weight, its gain is CYP1A2 and nothing else, and the macro is
back to zero rather than positive.** `src/ablncgc.py` with a row weight, which item 153 named as the
design error without fixing it. The base reproduces exactly (0.5462 in both runs), so the ladder is
comparable throughout.

    вес панели   доля таблицы   MACRO rho   против базы      1A2      2C9      2D6      3A4
    нет                    0 %     0.5462       +0.0000   +0.000   +0.000   +0.000   +0.000
    w = 1.0               74 %     0.5286       -0.0176   +0.010   -0.012   -0.048   -0.020
    w = 0.3               46 %     0.5392       -0.0070   +0.013   -0.011   -0.026   -0.005
    w = 0.1               22 %     0.5510       +0.0048   +0.023   -0.001   +0.000   -0.003

**Monotone in the weight on the macro and on three of the four enzymes.** The pre-registered branch
fired: an arm above 0.5462 says the weight was a cause, and w = 0.1 is above it. Item 153's reading
-- that a foreign protocol at 74 per cent of the table dominates the objective while the source
indicator has to undo it from a minority position -- is confirmed by the ladder rather than argued.

What may not be claimed, and it is the larger half. **+0.0048 macro is below the 0.007 floor.**
Correcting the weight moved the panel from clearly harmful to indistinguishable from not using it.
That is a repair, not a result.

The structure underneath is sharper than the macro. **CYP1A2 gains at every weight and is the only
enzyme that never loses** (+0.010, +0.013, +0.023), while **CYP2D6 is the largest casualty** at
-0.048 and its damage disappears exactly as the weight falls. At w = 0.1 three enzymes sit within
0.003 of the base and one is +0.023 -- so whatever the panel is worth, it is worth it on CYP1A2
alone.

Two hypotheses for that, both post-hoc and neither tested here. CYP1A2 has the thinnest own labels
relative to what NCGC supplies (1412 against 5169, a ratio of 3.7 against 2.4 to 2.9 for the
others) and the largest untapped gap in item 139 at 0.502. CYP2D6 is the enzyme whose signal is the
most mechanistically specific -- the Glu216 salt bridge and the basic amines -- which is exactly the
kind of structure a different protocol would scramble. Written down as candidates, not as findings.

**The macro floor is not the right floor for a single-enzyme claim.** 0.007 was measured on the
macro, which averages four enzymes and therefore has less variance than any one of them, so +0.023
on CYP1A2 is not automatically above its own noise. `--seeds 1,2,3` on the base and the w = 0.1 arm
is queued and is what settles it. Until it returns the honest summary is: the weight explained the
damage, and the only surviving candidate for a gain is one enzyme awaiting replication.

**158. Pooling's reversal is not column subsampling, and the screen works better with no pooling at
all.** `src/ablaux.py` with two arms the first design should have had. Both answers are large and
one of them is the best data-side result in this file.

    рука                      MACRO пара   MACRO rho     1A2     2C9     2D6     3A4
    независимо                    0.7125      0.5684  0.4966  0.6052  0.4116  0.7601
    независимо+скрининг           0.6865      0.5929  0.5087  0.6634  0.4177  0.7818
    пул mf1.0                     0.7379      0.5416  0.4807  0.5736  0.3821  0.7301
    пул+скрининг mf1.0            0.7150      0.5644  0.4803  0.6013  0.4226  0.7532

**Item 151's hypothesis is refuted.** The enzyme indicator being starved by `max_features=0.3` was
mechanical and testable, and at `max_features=1.0`, where the indicator is available at every split,
pooling gives **0.5416 against 0.5420** at 0.3 -- identical, and still 0.027 below the per-enzyme
model. Column subsampling has nothing to do with it.

What is left is the uncomfortable reading. **Pooling gains +0.0141 on HistGradientBoostingRegressor
and loses 0.027 on plain depth-5 trees, at either subsampling.** The submission is built around an
effect that survives one learner and reverses on the other, and its mechanism (item 132, contrast
through the enzyme indicator) does not explain why a different tree ensemble cannot use the same
indicator. That is now the largest open question in this file and it outranks anything queued.

**The screen, isolated from pooling, is worth +0.0245 of rank.** Screening rows added to a
per-enzyme model -- no indicator, no pooled table, just that enzyme's own curves plus its own
screening readings on compounds it has no curve for -- take 0.5684 to **0.5929**, and ST-RAE from
0.7125 to **0.6865**. Three and a half times the noise floor. Per enzyme: **CYP2C9 +0.058**, CYP3A4
+0.022, CYP1A2 +0.012, CYP2D6 +0.006, so CYP2C9 dominates exactly as item 152 pre-registered it
would (2.4 times more screening than curves, the highest ratio of the four).

Compare the four ways the screen has now been given to a model:

    пул+скрининг mf0.3           0.5727
    пул+скрининг mf1.0           0.5644
    независимо+скрининг          0.5929

**Every pooled variant is worse than the per-enzyme one.** Item 152 tied the screen to the pooled
table because item 132's contrast mechanism made that look natural, and that decision cost most of
the effect. The screen does not need pooling; it needed only to be a target instead of a feature.
Seeds 1-3 are queued.

**159. Early stopping was armed on `пул+TDI` and never fired, so item 125 is confounded only
mildly.** `verify/k41_earlystop.py`, one fold:

    рука       early_stopping   обучающих   do_early_stopping_   n_iter_        rho
    пул                  auto        5154                False   300/300     0.5730
    пул                 False        5154                False   300/300     0.5730
    пул+TDI              auto       10318                 True   300/300     0.5781
    пул+TDI             False       10318                False   300/300     0.5756

The threshold crossing of item 142 is confirmed -- `do_early_stopping_` is True on `пул+TDI` and
False on `пул`, exactly as the row counts predicted. But **`n_iter_` is 300 of 300**: the stopping
criterion was armed and never triggered, so the arm was not truncated.

The difference that remains is not truncation but the holdout: with early stopping armed,
scikit-learn reserves ten per cent of the rows as a validation set and trains on the other ninety.
That is worth **0.0025 of rank on this fold** -- real, in the direction that flatters the TDI arm,
and far too small to overturn item 125's null. Item 142 should be read as "the comparison had a
defect worth about 0.0025" rather than "the comparison was between two learners".

Setting `early_stopping=False` explicitly in every new arm stays right for a different reason: the
defect is small here and there is no guarantee it stays small on a table of 14424 or 60702 rows,
where the ten per cent held out is a much larger absolute number of molecules.

A defect in the check itself, recorded rather than quietly fixed: `k41` computes its ST-RAE column
by calling `fit_apply` with every row assigned to fold 0, which leaves the affine pair no training
folds and returns NaN. The rank column is unaffected and is what the item rests on.

**160. The screening level channel is chemistry on three enzymes and plate systematics on CYP3A4.**
`verify/k35_plate.py` re-runs item 13's channel contribution under a split blocked by assay plate
instead of by Butina cluster.

    фермент     протокол      без уровня   с уровнем     вклад
    CYP1A2      кластерный        0.4957      0.8626   +0.3669
    CYP1A2    по планшетам        0.5070      0.8555   +0.3485
    CYP2C9      кластерный        0.5972      0.8925   +0.2952
    CYP2C9    по планшетам        0.5649      0.8799   +0.3150
    CYP2D6      кластерный        0.4027      0.8283   +0.4256
    CYP2D6    по планшетам        0.3806      0.8210   +0.4404
    CYP3A4      кластерный        0.7646      0.9083   +0.1437
    CYP3A4    по планшетам        0.7482      0.6012   **-0.1470**

On CYP1A2, CYP2C9 and CYP2D6 the contribution holds within 0.02 of itself, so the channel carries
chemistry and item 13 survives a harder test than it was originally given. **On CYP3A4 it reverses
by 0.29**: block the plates and the level channel becomes actively harmful.

The explanation is item 129 and it is specific. CYP3A4 is the only enzyme whose training set
contains a separate analog campaign -- 530 compounds, run apart from the diversity screen -- and
compounds run together sit on the same plates. A Butina split leaks plate identity through
structural similarity, so the level channel on CYP3A4 has been partly reading which plate a
compound came from. That is the first measured instance of the plate confound this file has
suspected since item 35, and it is confined to one enzyme.

**161. The dead-zone-by-pairwise run is retracted: it lacked the control that would make it
readable, and its L1 was the wrong estimator.** `src/abldeadpair.py` returned the dead-zone pass
*costing* 0.027 to 0.031 of rank on the own booster, where `abldead` measured it gaining 0.027 to
0.036 on HistGradientBoostingRegressor. A sign reversal between learners is a finding only if the
alternative explanations are excluded, and two were not.

**No `L1 по метке` arm.** `abldead` has one and it matters: without it there is no way to tell "the
dead zone fails on this learner" from "this script's L1 is bad", and the whole comparison rests on
the L1 step.

**The L1 was bad.** `abldead` uses the pinned scikit-learn's `loss="absolute_error"`, which does the
line search that absolute-loss boosting requires. This script fitted trees to `sign(y - s)` under a
squared criterion, which is a different estimator: the leaf then holds the mean of a set of plus and
minus ones, the magnitude of the residual is discarded entirely, and the step size is arbitrary.
Friedman's LAD boosting fits the tree to the residuals and replaces each leaf value with the
**median** of the residuals in it; that is now implemented, the control arm is added, and the run is
requeued on two seeds. Nothing from the first run is carried forward, including its additivity
arithmetic.

**162. The dead zone in every ensemble member holds on four seeds, and it is worth most to the
members the submission weights least.** `src/abldzens.py`, seeds 1-3 against seed 0:

    рука                       MACRO пара   MACRO rho
    базовый ансамбль               0.6828      0.6002
    мёртвая зона в бустингах       0.6804      0.6007
    мёртвая зона везде             0.6615      0.6196

Seed 0 gave 0.6029 and 0.6203, so **+0.019 on all four seeds**, and applying the pass to the
boosters alone remains worth nothing. Per member, averaged over seeds 1-3:

    поферментно   0.5623 -> 0.6059   (+0.0435)
    GP            0.5535 -> 0.5997   (+0.0462)
    гребневая     0.5601 -> 0.5789   (+0.0188)
    пул           0.5758 -> 0.5932   (+0.0174)

**The two members that gain most are the per-enzyme booster and the Gaussian process**, and the two
that gain least are the pooled booster and the ridge. Read next to item 158 that is a pattern rather
than a coincidence: the per-enzyme model gains +0.0435 from the dead zone and +0.0245 from the
screening rows, both independently, while the pooled model gains little from either.

**163. The NCGC censored rows survive the re-weighting.** Fourth arm of item 157's ladder: all
54177 rows including the censored ones at w = 0.3 gives **0.5483**, against 0.5392 for the fitted-AC50
rows alone at the same weight and 0.5462 for the base. The censored rows are worth +0.0091 at w = 0.3
after being worth +0.0079 at w = 1.0, so "not active up to the top concentration" is information at
both weights. Per enzyme it buys CYP1A2 +0.025 and CYP2C9 +0.015 and costs CYP2D6 -0.028, the same
CYP2D6 sensitivity item 157 found. Still below the base on the macro; still the arm to keep if the
panel is used at all.

**164. The dead zone had been measured on four members while five are submitted, and on five it is
worth +0.0167.** `verify/k46_five.py`. Caught by reading rather than by running: `src/abldzens.py`
reads `oof_pool_all`, `oof_gp` and `oof_weak`, so its ensemble is the per-enzyme boosting, the
pooled boosting, the Gaussian process and the ridge -- **four**. Item 121 wired the joint-likelihood
trunk in as a fifth behind `--mode ансамбль5`, and that is what `src/submit.py` builds. So item
162's 0.6196, the best number in this file, belonged to a configuration that is not submitted.

Both `submit.py` (line 230) and `abldzens.py` combine members by an unweighted mean, so the five-
member ensemble is `(4 * четыре + ствол) / 5` exactly and can be composed from saved predictions in
minutes rather than recomputed in hours. The composition check passes: the plain five-member arm
reproduces item 120's rank of **0.6063** to the fourth decimal.

    состав                        пара      ранг
    четыре, базовый ансамбль    0.6819    0.6009
    пять, обычный (подаётся)    0.6824    0.6063
    четыре, мёртвая зона везде  0.6611    0.6198
    пять, МЗ в четырёх членах   0.6653    0.6230

**+0.0167 of rank in the configuration that is actually submitted**, on four seeds. Two things
follow. The trunk still adds +0.0032 on top of the dead zone, against +0.0054 without it, so the
pass takes a little of what the trunk was contributing but not most of it. And the number is a
**lower bound**: the trunk is not reprojected here, and every other member gained from the pass --
+0.0435, +0.0462, +0.0188, +0.0174 -- so a fully reprojected five-member ensemble should score
better than 0.6230. Doing that means re-running `src/trunk.py` under torch against the projected
target, which is the obvious next build.

**165. The per-enzyme noise floor is larger than the macro floor, and it had never been written
down.** Item 70 measured 0.007 macro at fixed seed, and that number has been used ever since to
judge per-enzyme claims. It cannot be: the macro averages four enzymes and therefore has less
variance than any one of them. Measured as the spread across four seeds at a fixed arm, on the
saved own-boosting sweep:

    фермент   sd по сидам
    CYP1A2         0.0061
    CYP2C9         0.0071
    CYP2D6         0.0049
    CYP3A4         0.0033
    МАКРО          0.0036

**The macro is quieter than three of the four enzymes**, and CYP2C9 is twice as noisy as CYP3A4.
Two consequences, both retroactive. A single-enzyme gain has to clear roughly 0.006 to 0.007 for
CYP1A2 and CYP2C9 and only about 0.003 for CYP3A4 -- so the same number means different things on
different enzymes, and this file has been treating them alike. And item 157's +0.023 for the NCGC
panel on CYP1A2 is about **four times that enzyme's own spread**, which is the check it was waiting
for: the claim survives, and the queued seed replication now tests reproducibility rather than
significance.

Differences between arms at the same seed share the fold structure, so their noise is smaller than
the spread of a single arm. The figures above are therefore conservative when used on paired
comparisons, which is the usual case here.

**166. "New information survives, rearrangement dies" is nearly the rule and breaks in two places
worth naming.** An outside reading sorted the log by what an arm *does* and proposed that
interventions supplying new information or a new function class survive while those permuting what
the model already has do not. Sorted against the file it is close, and the two exceptions are the
informative part.

**It fails on new observations.** The NCGC panel is new information by any definition -- 13126
structures, a different laboratory, measurements we did not have -- and it costs 0.0176 of rank at
weight 1.0 and clears zero only at weight 0.1 (item 157). ChEMBL was new information too and all
three handlings die under the affine pair (item 77). Meanwhile the screening table, which is also
new observations, is worth +0.0245 (item 158). The discriminator is not novelty, it is
**protocol**: the screen is the same laboratory, the same assay and the same compounds, the panel
is a different protocol with an offset of +0.44 to +0.87, and ChEMBL is selected by publishability.

**It fails on new columns.** The mechanistic block is not new information in the observational
sense -- every one of its thirty columns is computed from the same SMILES the descriptors are
computed from -- and it is worth +0.0163, up to +0.0313 in the test regime. Meanwhile FCFP,
chemprop, MoLFormer-XL and ChemBERTa are also computed from the same SMILES and all four die. The
discriminator here is not novelty either: it is whether the column encodes a quantity **not
derivable from the block already present**. Protonation state at pH 7.4 and salt-bridge geometry
are not in the 217 RDKit descriptors; a different fingerprint or a learned embedding is the same
structural information in another basis, and a tree ensemble that already has Morgan counts plus
descriptors gains nothing from a change of basis. Four re-encodings tried, zero survivors.

So the rule that fits the file is one step narrower, and it is a single idea in three costumes:

    новые НАБЛЮДЕНИЯ   выживают, если тот же протокол
    новые КОЛОНКИ      выживают, если величина не выводится из уже имеющегося блока
    новые ЦЕЛИ         выживают, если несут ПОСОЕДИНЕНИЕВУЮ структуру, которой нет у
                       аффинной пары --- ширина полосы, различимость пары, тождество фермента

Everything in the dead column is a monotone re-expression of what the model already outputs
(heteroscedastic layer 93, band edges 114, exact Bayes action 126, rank averaging 106) or a change
of basis (101, 117, 154), and item 77 explains why: the affine pair is a shrinkage fitted per fold
directly to the metric, so it already performs every monotone repair. **The prior for the remaining
five weeks is therefore not "prefer new data" but "prefer quantities the affine pair cannot
manufacture and the existing block cannot derive"** -- which is why the band width was worth more
than four external datasets.

**167. Every feature intervention that has ever worked here works on CYP2D6, and CYP3A4 has never
been helped by anything.** `verify/k47_perenzyme.py`. Item 165 measured the per-enzyme floor and
observed that this file had been judging per-enzyme claims by the macro figure. This is the audit
that follows from it: every ablation whose predictions are on disk, re-read per enzyme against that
enzyme's own floor. It computes nothing new.

Only the four-seed rows are trustworthy -- a single-seed per-enzyme difference has a spread of
0.005 to 0.007, so anything under about 0.015 on one seed says nothing.

    вмешательство            сидов      1A2      2C9      2D6      3A4
    механистический блок (81)    4  -0.0019  +0.0148  +0.0454  +0.0008
    пулирование                  4  +0.0157  +0.0071  +0.0374  -0.0057
    блок формы                   4  +0.0060  -0.0003  +0.0087  -0.0008
    блок кислот                  4  +0.0012  -0.0023  -0.0025  -0.0006
    пол этого фермента (165)        0.0061   0.0071   0.0049   0.0033

**CYP3A4 clears its floor on nothing and is pushed below it by pooling.** CYP1A2 clears it on
pooling alone. CYP2D6 clears it on three interventions out of four and by margins seven times its
floor. This is not a new measurement; it is the same measurements read against the right yardstick,
and the pattern was invisible while every one of them was reported as a macro.

Three specific consequences.

**The shape block is a CYP2D6 effect and item 119 closed it as a macro null.** +0.0087 on CYP2D6
against a floor of 0.0049, four seeds, while the macro is +0.0034 and reads as nothing. The
mechanism is in the block itself and was never used to read it: `shape3d.npz` carries
`bN_arom_ang_min` and `bN_arom_ang_mean`, the **angles of the basic nitrogen to the aromatic
system** -- that is the Glu216 salt-bridge geometry in three dimensions, and CYP2D6 is the only
enzyme it describes. The block was thrown in globally and scored globally.

**The acid block is null even per enzyme, and that refutes a hypothesis before it was built.** The
obvious next mechanistic block was CYP2C9's, whose site has Arg108 binding anionic ligands, so acid
strength and anionic fraction at pH 7.4 should pay there. Measured: **-0.0023 on CYP2C9 over four
seeds**, inside the floor, and nothing on the other three. The hypothesis is not refuted in general
-- the mechanistic block already carries `n_acid`, `n_tetrazole`, `n_sulfonamide` and `ph74_n_anion`,
so what is measured is the *marginal* value of six more acid columns -- but the cheap version of
"build CYP2C9's site block" is already done and already negative.

**The attention is inverted with respect to the headroom.** Item 139 measured the distance from our
model to a single calibrated screening column: CYP1A2 0.502, CYP2C9 0.340, CYP2D6 0.306, CYP3A4
0.289. So the enzyme that receives every working intervention is the one with the **second smallest**
gap, and the enzyme with the largest gap has received one. Item 134 said this in prose and was half
retracted; here it is arithmetic over four seeds.

**168. Proteochemometrics: the coordinate form is closed by arithmetic, the interaction form is
three-quarters unbuilt, and item 153 conflated them.** Item 153 wrote that the CYP2C19 arm was "the
only testable form of a proteochemometric enzyme coordinate available here". That was too strong
and the two halves need separating.

**The coordinate form is closed, and not by measurement.** With four enzymes, all four seen in
training, a one-hot indicator is a *sufficient statistic* for enzyme identity: any per-enzyme
function is exactly representable, so a vector of protein descriptors is a four-row lookup table --
a change of basis on the one-hot. Item 166 puts changes of basis in the column with four attempts
and zero survivors, and `src/ablcoord.py` measured this exact case: a two-column coordinate from the
measured screening correlations against the four-column one-hot, **+0.0007**. Protein descriptors
pay when the number of targets is large enough that sharing beats fitting each, or when a target
must be predicted that was never trained on. We have four targets and no zero-shot requirement.
Sequence, structure or cavity descriptors *of the enzyme* cannot add information the indicator does
not already carry.

**The interaction form is the thing that works, and it exists for one enzyme.** `src/feats.py`
opens with a docstring reading `Mechanistic feature block for CYP2D6`, and its geometric columns --
`topo_bN_to_arom_min`, `topo_bN_to_arom_mean`, `n_bN_geom`, `pharm_2d6`, `pharm_2d6_x_prot` -- are
Glu216 pharmacophore descriptors. They are ligand columns conditioned on a known active site, which
is what proteochemometrics is for when the protein set is small: not "which enzyme is this" but
"which ligand property matters for this enzyme". That block is worth +0.0454 on the enzyme it was
built for, and item 167 shows nothing else in this repository comes close.

CYP2C9's +0.0148 from the same block is an **accident**: `n_acid`, `n_tetrazole` and `ph74_n_anion`
are there for generic reasons and CYP2C9's site has Arg108. Nobody designed it and it clears the
floor twice over. That is the strongest available evidence that the deliberate version, for the two
enzymes that have nothing, is worth building -- and item 156 already supplies CYP1A2's recipe from
the data rather than from the literature: the tertiary aliphatic amine costs **-0.58 of pIC50 after
controlling for size, lipophilicity and aromaticity, on CYP1A2 and on no other enzyme.** A narrow
planar slot rejecting sp3 basicity is a testable feature specification, not an analogy.

What is genuinely untried and passes item 166's gate: a **docking score into the four crystal
structures**. A score for the pair (ligand, enzyme) is not derivable from the ligand block, because
it depends on the cavity -- unlike a fingerprint, an embedding or a protein coordinate, all of which
are functions of things we already have. There is not one mention of docking, SMARTCyp or a PDB
identifier in 167 items. It is also the largest compute commitment yet proposed: 4905 molecules by
four structures is about twenty thousand runs.

**169. One (E, h) does not describe CYP3A4 at all -- and the 530 unscreened compounds, where the
defect was expected, are not where it lives.** `verify/k48_calpop.py`.

The concern arrived from outside and was well aimed: `verify/g1_calib.py` fits the instrument
constants on the mask "has a curve AND has a screening reading", and item 129 established that
CYP3A4's training set is 1805 screened library compounds plus a 530-compound analog campaign that
was never screened. So on CYP3A4 the calibration is fitted on 1805 rows and applied to 2335, which
is item 83's error class. The other three enzymes have **zero** labelled compounds without a
reading, so this is a CYP3A4 question only.

    группа        n   медиана y     IQR   медиана полосы   доля y < pC0
    в скрине   1805       4.203   1.635            0.412         53.6 %
    кампания    530       4.454   1.189            0.315         44.0 %

**The proposed test cannot be run, and the reason is the concern itself.** The calibration residual
needs a screening reading to exist; the 530 have none. There is no population on which to measure
the misfit, because the missing reading is what defines the population.

So the question was asked one step back: how population-dependent is (E, h) at all? Split the
screened set at the median of the axes the campaign actually differs on, refit each half, and
transport the constants across. The bench reproduces the published constants first --- 0.728/1.260,
0.621/1.112, 0.867/1.242, 0.931/1.967, identical to `CAL_E` and `CAL_H` in `src/trunk.py`.

    фермент   раскол по полосе    раскол по уровню метки
    CYP1A2              +0.030                    +0.036
    CYP2C9              +0.035                    +0.019
    CYP2D6              +0.205                    +0.067
    CYP3A4              +0.141                    **+1.365**

**On CYP3A4 the Hill slope fitted on the lower half of the labels is 0.503 and on the upper half
2.516 --- a factor of five inside one enzyme's own screened population.** Transporting the low
half's constants to the high half gives a residual sd of 2.100 against that half's own 0.736, a
2.9-fold blow-up. That is not an instrument constant being applied to the wrong population. **It is
a single pair of constants failing to describe the population it was fitted on**, which is a worse
and more general defect than the one suspected, and it was never recorded.

**And the campaign is not where it bites.** Measured in units of the split just performed, the 530
sit 0.15 of a split away on label level and 0.05 on band width --- an order of magnitude closer to
the screened set than the halves are to each other. So "fitted on 1805, applied to 2335" is true
and is the *least* important instance of the problem. Extrapolating to the campaign is a small step
inside a map that is already wrong across its own domain.

**The biochemistry the outside reading offered turns out to be the explanation, in a different
place than proposed.** CYP3A4 is the textbook homotropic cooperative CYP: a large cavity that binds
two ligands, with a Hill slope that varies with occupancy rather than sitting at one value. That is
visible in our own published constant --- `CAL_H` is 1.968 on CYP3A4 against 1.261, 1.112 and 1.243
elsewhere, the outlier of the four --- and the 0.503-to-2.516 swing measured here is that
cooperativity resolved across the label range. CYP2D6's +0.205 on the band axis says it is not
purely a CYP3A4 phenomenon, but CYP3A4's is an order of magnitude larger.

**What this promotes.** A tree whose leaves carry local (E, h) instead of a constant --- proposed
from outside as the one unclosed form of informed boosting, at a medium prior --- is no longer a
speculative form. It is the direct fix for a defect now measured: the split by structure defines a
chemotype and the instrument physics is fitted inside it, which is exactly what a calibration whose
parameters swing five-fold across a population requires. The prior on it should be raised
accordingly.

**What this puts at risk.** Three things on CYP3A4 pass through these constants: the trunk's
screening channel (item 79, -0.0264 of pair), the band model the dead zone reprojects onto (items
148, 164), and the estimate of delta. None is invalidated by this --- the constants were fitted to
minimise residual over the whole screened set and remain the best single pair --- but every one of
them inherits a map with a 0.612 residual sd on CYP3A4, the worst of the four enzymes by a factor
of at least 1.2, and now with a named reason.