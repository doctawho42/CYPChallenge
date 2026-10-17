# What first place costs, and where the gap lives

`verify/k128_boardgap.py`, the nineteen committed board snapshots, latest pull 11 October. Every
number below comes out of that script; nothing is re-fitted and nothing needs the data files.

The project has used its board snapshots for two things — reading our own row's movement (items 313,
319, 344) and reading three entrants' published code (item 316). The other 272 rows are a sample,
and they answer a question no out-of-fold measurement can.

## Where we stand

    макро, 11 октября       ST-RAE     ро     место
    мы                      0.5357  0.7062   64 / 273
    первое (preheat-to-450) 0.3706  0.7844
    отрыв                  +0.1651 +0.0782

Item 344's two values still hold: ordering moved once, and the only thing that moved it was the
frozen-encoder probe. We are 64th now rather than 53rd, because the field grew and the teams above
us kept uploading.

## The decomposition, and it reverses the project's priorities

The board publishes a ST-RAE *and* a Spearman ρ on the same predictions. ρ sees only the ordering;
ST-RAE sees the ordering and the placement, because a prediction inside a compound's published band
scores exactly zero whatever its rank. So regress one on the other across the competitive band and
read the residuals.

    ST-RAE = 2.4721 − 2.6022·ро      n = 198, sd остатка 0.0975

    остаток       мы  −0.0987   (−1.01 sd)
                 1-е  −0.0603

> **WITHDRAWN 17 October, by item 354.** The conclusion below — that our placement is better than
> first place's and the axis is finished — is wrong, and the error is the comparison class. This
> regression is fitted across 198 entrants including a long tail with no post-processing at all
> (ST-RAE 0.69 to 7.03 at respectable rho), and that tail sets the slope, so beating the line says
> little. Matching on ORDERING instead — every entrant within 0.03 of our rho on an enzyme, best
> ST-RAE among them — says **+0.0872 of macro ST-RAE is unclaimed at the ordering we already have**,
> 53 per cent of the gap to first place: CYP1A2 +0.0701, CYP2C9 +0.0638, CYP2D6 +0.1429, CYP3A4
> +0.0721. The mechanism is in item 354: the board's MAE and R2 columns bound the test's label
> spread at 1.02 to 1.48 times our training spread, and `fit_shrinkage` can only express a shift of
> the mean, so the shipped lambda keeps a compression the wider test punishes. The arithmetic below
> is left standing as the record of what it said.

**Our placement is better than first place's, by 0.038 of ST-RAE.** We are a full standard deviation
better than the field at converting ordering into score, and the team in first is not. Of the 0.1651
gap, the cross-sectional slope attributes **−0.2035 to ordering** and gives **+0.0384 back on
placement**.

That number is the finding. The dead zone, the affine pair, the delta tilt, the monotone placement
budget, the threshold rules — items 128, 148, 256, 277, 282–285, 294, 313, 322 — are all placement,
and they are the bulk of this project's effort. They worked: that is what the −0.0987 is. And they
are now **finished as a source of places**, because the teams ahead of us are ahead on ordering while
being *worse* at placement than we are. Every further idea in that family is competing for a
resource we have already cornered.

A control that the regression is not an artefact: the local price of a board place around our row
comes out at **0.0031 of ST-RAE**, which is exactly what item 346 measured independently.

## What ordering has to deliver

To reach first place's ST-RAE while keeping our placement advantage:

    поперечным наклоном (−2.60)   +0.0634 ро
    нашим наклоном (пункт 322)    +0.1019 ро

The honest figure is somewhere between: the cross-sectional slope is contaminated by teams that are
good at both axes, and item 322's −1.62 is a within-team transform measured on one before-and-after.

Set that against the largest lever this project has ever *realised on the board*: the frozen probe,
**+0.0127 of ρ** (item 327, a 0.92 transfer from +0.0138 out of fold). First place is **5 to 8 such
levers.** For scale, here is everything we have that is positive and measured, in out-of-fold macro
rank:

    зонд, ядровая RBF вместо гребневой       +0.0166   пункт 360, 2 сида — ЗАМЕНА головы, не член
    зонд surrogate_adme, шестым членом       +0.0148   пункт 355, 4 сида (chemprop_medium недостижим)
    зонд chemprop_medium, шестым членом      +0.0138   пункт 324, знак 8/8, ОТГРУЖЕН
    панель NCGC, поферментно w0.1            +0.0112   пункт 356: ПРОЧИТАНО, умирает над составом
    minimol шестым                           +0.0048   пункт 361, но седьмым −0.0020 при 0/4
    стек: chemeleon + gbm-голова             +0.0047   пункт 346, знак 4/4
    surrogate_adme седьмым                   +0.0038   пункт 349, лучший одиночный кандидат
    chemeleon седьмым                        +0.0033   пункт 345

**The list does not sum, and since items 355–361 it is clearer than ever why.** Its top four entries
are the SAME slot — one sixth member built on a frozen encoder — and item 360's kernel head replaces
the head on the member that ships rather than adding one, so it is worth +0.0012 over item 355's
линейная and not +0.0166 on top of it. The error correlations between the candidates for that slot
are 0.907 to 0.939 (items 349, 356, 361), item 356 read the NCGC panel arm and watched it die over
the ensemble (−0.0092 as a seventh, sign 0/2) exactly as its pre-registered prediction said, and item
361 measured the same for `minimol` (−0.0020 at 0/4). The best *buildable* combination on this list is
the probe with a kernel head plus item 346's stack, at most +0.021 of out-of-fold rank if they add,
which items 346 and 349 say they mostly do not. **No combination of the project's historical lever
sizes reaches first place.** It needs a different class of model, not one more member.

The one thing on record priced at the right order of magnitude is item 341's: the frozen pretrained
encoder reaches 0.6261 of macro rank where our own from-scratch net reaches 0.5778, so **pretraining
is worth at least +0.048 on our own folds** — three to four times the best lever we have shipped.
And the configuration that converts it is published in outline: item 343 records that the best
readable entrant's model is a fine-tune warm-started from that checkpoint **with seventeen heads**,
and item 341 records what it is trained on — "four nets averaged and challenge+ChEMBL+qHTS behind
them". That is the only lever of the right size anywhere in this record, and this project has never
built it. `docs/transfer_axis.md` is the audit of why, and of what is and is not closed around it.

## Where the gap lives, by enzyme

    фермент   нас      ST-RAE     ро   |   1-е  ST-RAE     ро    лучшее ро   отрыв ST-RAE   доля
    CYP2D6    # 72     0.6786  0.4359      0.4135  0.5428     0.5745        +0.2651       38.5 %
    CYP1A2    # 65     0.5656  0.7513      0.3846  0.8561     0.8561        +0.1810       26.3 %
    CYP3A4    # 75     0.4729  0.8146      0.3515  0.8741     0.8741        +0.1214       17.7 %
    CYP2C9    # 68     0.4258  0.8231      0.3056  0.8808     0.8836        +0.1202       17.5 %

**CYP2D6 is 38.5 per cent of the gap on its own**, and the field's best ρ there is only 0.5745 — so
it is not a label-resolution ceiling that everyone has hit, it is 0.1386 of ordering that someone
else has found and we have not. This is the enzyme the journal already says has an address: item 134's
potent basic amines and the Glu216 salt bridge, and item 239's measurement that permuting the ten
salt-bridge columns costs 0.1545 of rank on CYP2D6 above a size-matched random block — the largest
single-block dependence in the feature matrix. `docs/teammate_asks.md` reached the same conclusion
from out-of-fold numbers alone ("2D6 — самый слабый и у него адресный хемотип"); the board confirms it
with 38.5 per cent of the money.

**Two things this does NOT license, because both have already been settled.** First,
`results/perenzyme.csv` is not a shortcut. Item 167 built it by re-reading every ablation per enzyme
against item 165's per-enzyme floors, and it carries a column literally called *скрыто макро* —
hidden by the macro. Four-seed rows with a CYP2D6 gain above that enzyme's 0.0049 floor are there,
and the biggest ones are already shipped (the mechanistic block +0.0454, pooling +0.0374). The one
that looks most like a free lunch — the shape block, CYP2D6 +0.0087 at sign 4/4, closed by item 119
on a macro of +0.0034 against a 0.0036 floor — was reopened on exactly that argument and then closed
properly over the ensemble: *"the enzyme the shape block helps most standalone is the enzyme that
dilutes it most in the ensemble"*, with CYP3A4's +0.0000 serving as the harness control. So a
member-level per-enzyme gain still has to survive the composition, and on CYP2D6 specifically the
composition is where they die.

Second, item 167's closing conclusion needs revising against the board rather than repeating. It
reads: *"The attention is inverted with respect to the headroom … the enzyme that receives every
working intervention is the one with the second smallest gap"*, from item 139's distance to a
calibrated screening column (CYP1A2 0.502, CYP2C9 0.340, CYP2D6 0.306, CYP3A4 0.289). That ranks
distance to an **oracle**. The board ranks distance to **competitors**, and it puts CYP2D6 first by a
wide margin — 38.5 per cent of the gap, against CYP1A2's 26.3. For a competition the second is the
operative quantity, so the attention CYP2D6 has received is not misallocated; what item 167 actually
identified is that CYP1A2 has the most room nobody can reach, which is a different claim.

Note also that our board ordering and our out-of-fold ordering disagree about which enzyme we are
good at. Out of fold: 3A4 0.796 ≫ 2C9 0.665 > 1A2 0.547 > 2D6 0.471. On the board: 2C9 0.8231 > 3A4
0.8146 > 1A2 0.7513 > 2D6 0.4359. Different populations, 375 compounds against our training rows, so
it is not a contradiction — but any per-enzyme decision taken on the out-of-fold ordering alone is
being taken on the wrong ranking of our own strengths.

## Is first place a moving target? Barely

    снимок               n      1-е      2-е     10-е       мы
    2026-09-15         164   0.3944   0.3952   0.4428        —
    2026-09-25T0150Z   239   0.3680   0.3776   0.4361   0.5760
    2026-09-26T0727Z   227   0.3680   0.3776   0.4315   0.5381
    2026-10-02T1733Z   248   0.3772   0.3776   0.4367   0.5794
    2026-10-06T0934Z   246   0.3706   0.3776   0.4301   0.5357
    2026-10-11T1926Z   273   0.3706   0.3776   0.4301   0.5357

**First place moved −0.0238 of ST-RAE in twenty-six days and second place has not moved at all since
24 September** (0.3776, eight snapshots). Tenth moved −0.0127. The target is close to static, which is
the one piece of good news in this file: we are not chasing a receding leader, we are 0.165 behind a
nearly stationary one.

Our own row moved −0.0403 over the same window — and that includes one upload on 2 October that
*lost* 0.0437 (0.5357 → 0.5794) before being reverted by the 6 October pull. The re-upload rule is
one valid submission per twelve hours with the **last valid one counting**, so a bad upload left
standing at the close is a real risk, not a recoverable experiment.

**Submission closes 3 November 2026** (README, `docs/md_task_spec.md` §10). That is seventeen days
from this snapshot, against a required +0.165 of ST-RAE and a historical rate of −0.040 per three
weeks.

## One adopted gain has never reached the board

Item 344 tabulates nine uploads, the last submitted **3 October**. Items 345 to 351 all come after it,
and one of them adopted a change: item 346, *"the three sub-floor members ADD, and the sum clears both
floors: +0.0047 of macro rank and −0.0085 of pair, sign 4/4 on each. It is the first adoptable gain
since item 324."* No upload is recorded anywhere after item 344's ninth, and METHOD.md still documents
the composition as five members plus the probe (line 403) rather than item 346's eight. Item 346
itself prices it: *"At the measured 0.0031 of ST-RAE per board place this is about three places."*

So before any new modelling, there is a measured, controlled, already-adopted gain of about three
places sitting unuploaded with seventeen days left. That is the cheapest board movement available and
it needs no new ideas at all.

**And it could not have been uploaded, which is why it wasn't.** `src/submit.py` held exactly one
extra member — `_combine(parts, e, extra=None)`, with `PROBE_EMB` hard-wired to
`emb_chemprop_medium.npz` — so item 346's composition had no route to a submission file. Two thirds
of it is a parameterisation and is now done: `_combine` takes an `extras` list beside `extra`,
`_probe_load` takes `emb/emb_te/tag` instead of loading one fixed pair, and `oof_predictions` threads
it. Added as a second parameter rather than by widening the first, because telling a list of four
arrays from a list of lists of four arrays means guessing about shapes. Inert by default, and
`tests/test_extra_members_inert.py` pins both halves — nothing moves on the shipped paths, and the
new members reach **every** enzyme rather than being eaten by `_keep`/SOLO, which is the hazard item
324's own test was written for. The pre-existing `tests/test_probe_member_inert.py`, which compares
against the pre-probe module recovered from git, still passes.

What is deliberately **not** done, because it cannot be executed without the data and this is the
file that builds the actual submission: `main()` is still wired for one extra member, and item 335's
gbm head has no device in `submit.py` at all. So the remaining work for action 0 is (a) thread a list
through `main()`'s five probe sites — load, test prediction, out-of-fold, the twin control, the
`probe_parts.json` dump — and (b) either add the gbm head's device or ship `chemeleon` alone, which
is +0.0033 of rank and −0.0054 of pair of item 346's +0.0047/−0.0085, about 1.7 places of the 2.7.
A bad upload left standing at the close is not recoverable, so this part wants a machine that can run
`make submit` and read its validators before anything goes up.

## What this says to do

0. **Upload item 346's stack.** It is adopted, controlled at sign 4/4, worth about three places, and
   has not reached the board. Nothing else on this list is free.
1. ~~**Stop spending on placement.**~~ **Withdrawn — the opposite is true.** Item 354: at our own
   ordering the board says 0.0872 of macro ST-RAE is unclaimed, which is 53 per cent of the gap and
   the single largest actionable quantity in this project. Fit the per-enzyme SCALE against the
   test's label spread (bounded off the board itself) rather than against our training spread. It
   cannot cost ordering: expansion is strictly monotone, so Spearman is invariant exactly.
2. ~~**Read what is already measured**: `results/preds/oof_ncgc_per.json`, +0.0112 of macro rank at
   2/2.~~ **Done — item 356.** It reproduces at the member level (+0.0113 at 2/2, shuffle control
   firing) and **dies over the ensemble**: −0.0092 as a seventh member and −0.0044 as a substitution,
   sign 0/2 on both, error correlation 0.935 with the existing composition. The pre-registered
   prediction P0 in `docs/transfer_axis.md` held on both halves. Nothing to upload from it.
3. **Put everything else into ordering, and weight it toward CYP2D6**, which is 38.5 per cent of the
   gap and where the field proves there is signal we are missing.
4. **The only lever sized for first place is representational**: a warm-started multi-head fine-tune
   over the release plus ChEMBL plus the qHTS panel, which is the published shape of the best readable
   entrant's model and the one configuration items 319 and 343 did *not* measure — item 319 added one
   head to a trunk trained from scratch, item 343 fine-tuned one head from a warm start. Nobody here
   has built both halves at once.

## Is the gap a label ceiling? No, and that is now a number

`verify/k138_noiseceiling.py`, item 362. The organisers publish a per-compound assay sd for every
label, and the band is exactly 2·1.96 times it. Converted into a ceiling on **Spearman** — not on
R², which `src/ceiling.py` has had all along and which attenuates differently — the labels permit:

    фермент   потолок ро   наш ро   лучший на доске   наш %   лучш %
    CYP1A2       0.9892    0.7513        0.8561        76 %    87 %
    CYP2C9       0.9665    0.8231        0.8836        85 %    91 %
    CYP2D6       0.9786    0.4359        0.5745        45 %    59 %
    CYP3A4       0.9783    0.8146        0.8741        83 %    89 %
    МАКРО        0.9782    0.7062        0.7971        72 %    81 %

**The +0.1019 of ρ the table above demands is 37 per cent of the +0.2719 the labels leave
unclaimed**, and three quarters of what is unclaimed sits on CYP1A2 and CYP2D6 — the same two
enzymes this file already singles out. The band's scale would have to be understated by a factor of
3.2 for the board's best ordering to already be at the ceiling, on a quantity item 220 measured as
3.8 times too *wide*. Item 354's bounds on the test's label spread raise the ceiling further, to
0.9861. And item 362's sharper half: the published error sits almost entirely on the weak compounds
(95 per cent of the noise variance), and we take **44 per cent** of the ceiling on the potent half
against 58 per cent on the weak one — the deficit is worst where the labels are cleanest.

So nothing in this file can be explained away as label resolution. The constraint is modelling.

And the arithmetic should be said plainly rather than implied: on this evidence first place is not
reachable by the kind of step this project has been taking. It needs the 0.048-class intervention to
work, and then it needs most of the smaller ones to add on top of it without the redundancy that ate
items 176, 182, 191 and 194. Item 362 removes the one excuse that would have made that shape
acceptable — the labels permit four times the ordering first place needs from us — so the target is
honestly shaped as a modelling gap, not a measurement one.
