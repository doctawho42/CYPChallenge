# Verification

Fourteen scripts. The `f*` run was a sweep over everything that had been computed and
written by that point; `g*` answers four questions raised against the document. There are
no dependencies between scripts, so the order is free. All of them expect the data in
`data/` and the organisers' metric code (`git submodule update --init`).

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

Four claims were established for the first time: the MNAR to MAR transition was confirmed
by simulation; the Spearman gain from the mechanistic block was confirmed by bootstrap;
the Bayes optimum under ST-RAE was verified numerically; and our cross-validation was
measured to be 0.153 harder than the test set.

No leakage was found in the pipeline. There is no contamination between training and test:
the overlap by name, by canonical SMILES and by SMILES with stereochemistry stripped is
zero across all three training files.

## Three more, found while setting up the environment

These are documentation errors rather than analysis errors — the code was right and the
prose was wrong — and they are corrected in the English documentation but **not yet in
the PDF**:

9. The mechanistic block has **30** features, not 28. `data/mech_names.csv`, written by
   `feats.py`, has always listed 30 and matches the copy committed in `results/`.
10. There are **fourteen** verification scripts, not twenty-three, in both this file and
    the top-level README.
11. `f12_cvhard.py` could never have run to completion as committed: line 34 referenced an
    undefined `fRES`. Fixed, and the script now reproduces the documented -0.153.
