# Transfer, multi-task and somebody else's data: what is closed, what is open, and the result that was already on disk

Written in answer to a teammate's question — *were the external-data ideas closed too early? We have
other people's panels and structures; we could pretrain an encoder on them and transfer it, or make
one trunk solve two tasks at once.* The instinct is right, and the evidence for it turned out not to
be a new idea at all.

Nothing below is a new run. Everything is either a quotation from the journal, a number recomputed
from a committed cache, or a prediction written down before the measurement that would test it.

## The headline: two panel arms are sitting in committed caches and the journal reports neither

`src/ablncgc.py` wrote four caches into `results/preds/`. The journal reports two of them —
`oof_ncgc.json` at item 153 and `oof_ncgc_w.json` at item 157 — and reports neither
`oof_ncgc_w123.json` nor `oof_ncgc_per.json`. Both files are named nowhere in `verify/README.md`,
`METHOD.md` or `docs/`, and neither arm's result appears anywhere in any form. (Git says nothing about
when they were produced: this repository's history begins with one squashed import of 463 files, so
every cache carries the same commit date.)

The only trace of the second file is item 187, which records that the per-enzyme panel arm *"died
instantly three separate times"* on a parser defect in `src/ablncgc.py`, *"twice inside a queue whose
wrapper reported success"*, and that the parser was replaced with a regular expression. The repaired
run's output is what that file holds. Nobody read it.

**One: the replication item 157 queued did return.** Item 157 said `--seeds 1,2,3` on the base and
the w = 0.1 arm *"is queued and is what settles it"*; item 165 said the +0.023 on CYP1A2 is about
four times that enzyme's own spread and that *"the queued seed replication now tests reproducibility
rather than significance."* Four seeds, pooled arm, read out of the caches:

    ранг, пул + NCGC w0.1     1A2       2C9       2D6       3A4     MACRO
    сид 0                  +0.0228   -0.0007   +0.0001   -0.0030   +0.0048
    сид 1                  +0.0177   -0.0118   -0.0146   -0.0041   -0.0032
    сид 2                  +0.0016   +0.0015   -0.0108   +0.0093   +0.0004
    сид 3                  +0.0142   -0.0134   +0.0073   -0.0023   +0.0014
    среднее                +0.0141   -0.0061   -0.0045   -0.0000   +0.0008
    знак                       4/4       1/4       2/4       1/4       3/4
    пол                     0.0061    0.0071    0.0049    0.0033    0.0036

So item 157's surviving claim survives: **CYP1A2 +0.0141 at sign 4/4, 2.3 times its own floor** — but
at 60 per cent of the headline figure, and the macro stays at zero (+0.0008, below 0.0036). That is
the honest four-seed version of item 157, and it has been computable for a month.

**Two, and this is the larger half: the per-enzyme arm inverts item 157's conclusion and is the best
external-data number this project has.** Item 157 concluded *"whatever the panel is worth, it is
worth it on CYP1A2 alone"* and *"+0.0048 macro is below the 0.007 floor … a repair, not a result."*
Both sentences are about the POOLED model. The per-enzyme arm, against the per-enzyme base, with the
shuffled-panel control run beside it on both seeds:

    поферментно + NCGC w0.1        1A2       2C9       2D6       3A4     MACRO   перемешанный
    ранг, сид 0                 +0.0126   +0.0194   +0.0108   +0.0123   +0.0137       -0.0009
    ранг, сид 1                 -0.0006   +0.0291   -0.0054   +0.0124   +0.0088       -0.0036
    ранг, среднее               +0.0060   +0.0243   +0.0027   +0.0123   +0.0112       -0.0022
    ST-RAE после пары, сид 0    -0.0027   -0.0263   +0.0024   -0.0196   -0.0115       +0.0032
    ST-RAE после пары, сид 1    +0.0001   -0.0214   +0.0130   -0.0151   -0.0059       +0.0041

The macro agrees in both currencies on both seeds, and the shuffled panel goes the wrong way in both
currencies on both seeds. The macro rank gain is **+0.0112**, three times the four-seed macro floor;
the post-pair macro ST-RAE gain is **−0.0087** mean, which clears item 70's 0.007 fixed-seed floor on
seed 0 and not on seed 1. It is carried by **CYP2C9 (+0.0243 of rank, −0.024 of ST-RAE) and CYP3A4
(+0.0123, −0.017)** — the two enzymes the pooled ladder had the panel hurting or doing nothing for —
while CYP1A2, item 157's whole story, contributes +0.0060 at sign 1/2 and CYP2D6 is the one cell that
disagrees between the currencies.

The mechanism is not mysterious and the project measured it twice already. Item 157's damage was *"a
foreign protocol at 74 per cent of the table dominates the objective while the source indicator has
to undo it from a minority position"*; item 188 measured that pooling costs 0.020 of rank because a
pooled model spends depth isolating the enzyme before it can model anything conditional on it. A
per-enzyme fit has no indicator to starve and no enzyme to isolate, so the panel's rows arrive in the
one place where a foreign protocol cannot be confused with an enzyme contrast. The pooled ladder was
measuring two things at once.

**What may not be claimed.** Two seeds, where the convention is four. The gain is over the
**per-enzyme member**, not over the shipped composition — and `docs/teammate_asks.md` states the
ceiling that kills exactly this kind of number: a member-level δ sinks over the ensemble both by
averaging (~δ/N) and by redundancy, 3 to 22× on the fused enzymes. Items 290 and 319 are the worked
example of a real channel dying on that ceiling. So this is a candidate, not a result: the decisive
measurement is the same one item 290 applied and item 324's probe passed — nested over the shipped
composition, four seeds.

What does **not** need redoing is the currency. `src/ablncgc.py:269-278` scores exactly as
`verify/k106_probe.py:score()` does — per enzyme on that enzyme's masked vector, `spearmanr(y, p)`
for rank and `strae` after `fit_apply(p, lo, hi, fold, u)` for the pair, macro as the mean over the
four — so these numbers are in the same units as item 324's +0.0138 and item 346's +0.0047 and are
directly comparable to them. Both halves of the decisive measurement are therefore cheap: two more
seeds, and the nesting, which needs no retraining at all because the arm's out-of-fold predictions
are in the cache. That is the same route item 117 used to re-read item 68 and item 290 used to close
the external head.

## The grid the question is really about: source × role

The project has never tabulated its external-data results by ROLE, which is why "we closed external
data" and "we closed pretraining" sound like the same sentence. METHOD.md already makes the
distinction once, for the chemprop family — concatenation into the same learner is closed (items 61,
117, 154), a ridge head on the frozen encoder ships (items 324, 325, 327) — and calls the difference
between those two "the finding". Extended to every source:

| Source | rows, pooled | rows, per-enzyme | columns | own head on the trunk | pretrain → frozen member | pretrain → fine-tuned |
|---|---|---|---|---|---|---|
| in-release screen, 17504 | 158: +0.0245 | 158 | — | **ships**: 79, 80, 121 | — | — |
| in-release pre-incubation arm | 240: inside floor | — | — | 243: AUC yes, MCC no | — | — |
| ChEMBL 8068 (CheMeleon baseline set) | 61, 77: harmful after the pair | — | 61, 117, 154: −0.0075 | **319: −0.0005, 1/4** | — | — |
| **NCGC qHTS panel, 54177 rows** | 157: +0.0048, floor 0.007 | **+0.0112, 2/2, unread** | **never** | **never** | **never** | — |
| organisers' second CYP release | 201: impossible, 14/1340 overlap | — | — | — | — | — |
| published checkpoints (15 + 2 + 3) | — | — | 61, 117, 154 | — | 324–350: 18 cases, 5 helped | **343: −0.0648** |
| structures / cocrystal | — | — | 85, 111, 112 | — | 350: 0/4 on three 3D checkpoints | — |

Four cells are empty, and all four are the same source.

## Half one of the question: "one trunk, two tasks, the second external" is closed, and closed by name

This is not a near neighbour of the proposal. Item 290 names it:

> `src/trunkext.py` (item 66) — the shared trunk with 8004 external ChEMBL CYP pIC50 rows in the
> auxiliary head instead of the screen — **is the artefact that most directly instantiates the "one
> shared trunk, two tasks, the second external" proposal**

It followed the project's own pattern: built at item 66, left "still not settled" because it had only
ever been scored in raw ST-RAE, which item 77 established is not evidence. Item 290 closed it over
the ensemble by re-reading committed predictions. Item 318 then built the real vehicle — a third head
beside `head_pic` and `head_scr`, proven inert bit-for-bit on the discriminating device, with the
generator wound back around the new head because a head that merely comes last still moves every
prediction through the dropout masks (removing the wind-back moves predictions by 0.92 pIC50,
measured by planting the defect). Ship gate pre-registered at **≥ +0.0105 macro Spearman** over the
ensemble, nested, sign ≥ 3/4, deliberately above item 290's own +0.0085.

Item 319 ran it. Four seeds, both arms:

    C − A   -0.0005   знак 1/4   трёхголовый против штатного — в 4.5 раза НИЖЕ пола 0.0036
    C − B   +0.0043   знак 4/4   сам канал: половина от +0.0085 пункта 290
    B − A   -0.0048   знак 0/4   цена вставки 8004 строк

The closing argument does not depend on the vehicle, which is why it closes the proposal and not just
that build: **even at zero insertion cost the channel gives +0.0043, still below the gate, so no
better vehicle rescues this line.** And the channel halved exactly as the vehicle improved — a weak
trunk has more room to be helped, so part of item 290's +0.0085 was the weakness it repaired rather
than the information it added. Two further facts make a rematch unattractive: the three-head trunk is
0.0129 *below* the member it displaces even with the channel on (0.5837 against 0.5966), and the
board price differs by a factor of 8.6 between currencies — zero places on rank, at least one place
lost on post-pair ST-RAE, same arm, same run.

What item 319 leaves open, in its own words: `lam_ext = 3.0` was carried in, so it measures λ = 3.0
and not the best λ. A grid could only have been scored on the channel-on arm, which is item 245's
defect by construction, so the carried-in value is the honest choice rather than a gap.

**And the multi-task family is already exploited where it pays.** Item 188 put the four enzymes in one
tree with four values per leaf: +0.0069 of macro rank over the per-enzyme reference and +0.0273 over
pooling — CYP1A2 +0.0208, CYP2C9 +0.0075, CYP2D6 +0.0166, and CYP3A4 **−0.0171**, the enzyme with the
most labels. That is ordinary negative transfer on the data-rich task, and it is the shape to expect
from any further multi-task arm here: sharing helps the tasks with least data and can cost the task
with most.

## Half two: pretrain-then-transfer is the project's best external result, and the record names the one form left

The frozen-encoder device is the only thing this project shipped that moved the board's *ordering*.
Item 324: a RidgeCV head on 600 frozen dimensions of somebody else's encoder scores standalone macro
rank 0.6261 — above every member we have — and adds **+0.0138 at sign 8/8**. The control that makes
it a result about the weights rather than about adding a member: the same RidgeCV on our own
DESC+MECH block is **−0.0151 at sign 0/4**. Item 327 measured the transfer to the board at 0.92
against 0.52 and 0.437 for post-processing arms, and item 344 records that across nine uploads our
board Spearman took exactly two values — ordering moved once, and this is what moved it.

Three measured constraints shape anything new here.

**Fine-tuning is the wrong instrument — in the configuration that was measured.** Item 343: the best
fine-tuned arm reaches 0.5613 against **0.6261** for the same checkpoint left frozen — 0.0648,
eighteen times the rank floor — with discriminating same-architecture random-weight controls at
0.4942 (`rand600`) and 0.4833 (`rand2048`). But item 343 states its own scope in a sentence that is
easy to read past, and it is the most important sentence on this axis:

> his headline "freeze, don't fine-tune" is about naive fine-tuning on the challenge's labels, while
> his best-ranking chemprop model IS a fine-tune warm-started from this checkpoint with **seventeen
> heads**. The configuration matters and the naive one loses.

So what is closed is the naive instrument: one head, our labels, from a pretrained start. The
configuration that actually outranks us is a warm-started **multi-head** fine-tune — and item 341
records what it is trained on: `briford`, 16th at ρ 0.7516 against our 0.7062, uses "four nets
averaged and **challenge+ChEMBL+qHTS** behind them". That is the teammate's two proposals in one
object: transfer the encoder *and* make one trunk carry many tasks, over the union of the release,
ChEMBL and the NCGC panel. This project has measured the two halves separately and never together —
item 319's head was one extra head on a trunk trained from scratch, and item 343's fine-tune was one
head on a warm start. The cell with both is empty, it is the published configuration of the best
readable entrant, and item 341 already measured the size of the prize on our own folds: "pretraining
is worth at least +0.048 here".

**The admission screen is standalone rank, and it is now pinned.** Item 348 read a threshold of 0.60
off thirteen candidates and said the next candidate would be its test rather than its confirmation.
Five arrived; the band is empty and 0.0089 wide — `minimol` at 0.5990 below, item 342's graph net at
0.6079 above — over eighteen cases with no exception, Spearman(standalone, Δ-as-member) = +0.8184.

**The class of the pretraining target decides it, and same-class encoders do not add.** Item 349's
table, extended by item 350:

    эндпойнт, измеренный          n=2   помогли 2
    эндпойнт, дистиллированный    n=1   помогли 1
    дескриптор                    n=5   помогли 0
    представление чужой модели    n=1   помогли 0
    языковая модель над SMILES    n=2   помогли 0

And the ceiling that prices the proposal honestly: `surrogate_adme` is the best single candidate on
record (+0.0038, sign 4/4, sd 0.0007) and adds **+0.0007** to item 346's stack — a fifth of the floor
— because `chemeleon ~ surrogate_adme` sits at **0.914**, the most correlated pair measured here,
recovering 60 per cent of its summed parts. The mechanism that explains why an endpoint-class encoder
works is the mechanism that stops a second one adding.

**And there is a cheaper candidate than building an encoder, declined once for a reason that applied
to a different experiment.** `src/embed.py:24` records it: *"The CYP-fine-tuned baseline in the same
organisation was NOT used, because a model already trained on CYP pIC50 data would make the
comparison circular."* That was correct for item 61's question — does a *learned* representation of
descriptor information beat handing the information over — where a CYP-trained encoder would answer a
different question. It is not the question item 324's device asks, which is only whether a frozen
prediction adds over our members. The organisers' own CheMeleon CYP baseline is an endpoint-trained
checkpoint in the organisation whose training rows `data/fetch.sh` already downloads, with **zero of
the 750 blinded test molecules** in that set and 64 of our 4905 training molecules (item 318) — so the
circularity is bounded at 1.3 per cent of our rows and is not a rules question at all. It is absent
from item 350's enumerated stock (fifteen organiser checkpoints, two third-entrant, three other public
encoders), and it has never been scored. If that repository publishes the weights beside the training
data it publishes, this is one download and one `k121` run — cheaper than every remaining descriptor
checkpoint and in the one class that has never failed.

Item 350 then states the remaining cell itself, and never scopes it:

> the next encoder worth building would have to be pretrained on **measured activity we are allowed
> to use**, which is a different and much larger piece of work, and `surrogate_adme` shows it would
> also have to be **decorrelated from the three ADME encoders we already hold**

That is the teammate's proposal, named in the record as the one live option on this axis, with both
of its conditions attached and no measurement against it. It is not closed. It was deferred on cost.

Independently corroborated, recorded at item 316 because someone else paid for it: entrant adlvdl's
stated rule is *"sharing a learned encoder helps; pooling rows of a frozen one does not"* — the same
shape as our items 61/117/154 against 290, reached by a team that never read this journal.

## Why the NCGC panel is the source that fits, in every role

`src/fetch1851.py` refetches it from PubChem PUG REST (AIDs 410, 883, 891, 899, 884, plus the 411
luciferase counter-screen) and `src/ncgcmerge.py` folds it into `data/ncgc/panel.csv`. It is the only
corpus in the ledger that is simultaneously refetchable from a committed script and measured activity
on our own scale: one laboratory, one protocol, one campaign, the whole library at sixteen
concentrations, no publication filter, medians 5.30 (CYP1A2) and 5.10 (CYP3A4) against our labels'
5.13. 54177 rows over 13126 structures, 18711 with a fitted AC50, inactives kept as censored
measurements, a fifth isoform (CYP2C19) we do not have, and **zero of the 750 blinded test molecules**
(item 144).

And the reason it died as pooled rows is a statement about its scale meeting our objective, nothing
else. Item 144: the per-enzyme offsets are +0.442 to +0.868 and only **11 to 48** molecules bridge
the two protocols per enzyme, giving the offset a standard error of 0.21 against a quantity of 0.45,
so the arm that subtracts a fitted shift *was never built* and a source indicator replaced it. Item
203 measured 0.0075 for a source that is calibrable, on 2334 shared compounds. Item 201 states the
rule: overlap is a property of the library's kind, not of its size, so more NCGC rows never repair it.

In the frozen-encoder role **none of that can act**. The encoder is trained on the panel out of
process, the weights are frozen, and the head that reads them is refit in fold on *our* labels, which
is where the scale comes from. There is no offset to calibrate, so the eleven-compound overlap — the
binding constraint on every rows route and on the head route — is irrelevant. That is the teammate's
point, and it is correct.

## Pre-registered predictions

Written before anything is built, against the thresholds the record already fixes. Honest priors in
brackets; two of them disagree with wanting this to work.

**P0 — the unread arm, and it goes first because it is nearly free.** The per-enzyme panel arm
reproduces on seeds 2 and 3 at sign 4/4 on macro rank, and **does not** survive nesting over the
shipped composition above the 0.0036 floor. *[the first likely — two seeds, a firing control and two
agreeing currencies; the second likely too, because `docs/teammate_asks.md`'s second ceiling has
eaten four member-level gains of this size, items 176, 182, 191 and 194.]* If both halves hold, the
deliverable is a per-enzyme feature-side gain that must then be fought for over the ensemble, and
CYP2C9 is where to fight — it is the arm's own largest cell, and item 324's probe gained most there
too (+0.0283).

**P1 — the encoder's screen.** An encoder pretrained on NCGC pAC50 over the five isoforms, frozen,
read by `probe_oof`, scores standalone macro rank **≥ 0.60** and passes item 348's threshold.
*[likely: the class is measured activity, which is 2 for 2, and the panel is larger than anything the
three successes were trained on.]*

**P2 — the stack, and this is the prediction that matters.** Added to item 346's stack its
contribution lands **below the 0.0036 macro floor**, because it is the fourth endpoint-class encoder
and item 349's mechanism is monotone in class identity: its error correlation with `chemprop_medium`,
`chemeleon` and `surrogate_adme` will be **≥ 0.88**, and same-class pairs recover 57 to 60 per cent of
their summed single contributions. *[likely, and not in tension with P1 — that is exactly the pattern
`surrogate_adme` set: best single candidate on record, +0.0007 in the stack.]*

**P3 — the one place the encoder can pay.** On a single enzyme it clears that enzyme's own floor at
sign ≥ 3/4 on fresh seeds, and **it is CYP2C9 rather than CYP1A2**, because the unread per-enzyme arm
puts the panel's information there (+0.0243 against CYP1A2's +0.0060) and because CYP2C9 has the
fewest own labels of the four (1285). *[a coin-flip, and the most informative number in the design
either way. Note this contradicts what item 157 would have predicted, which is the point of writing
it down.]*

**P4 — the head route, for completeness.** If the panel is put through `head_ext` instead, the channel
C − B exceeds item 319's +0.0043 (more rows, one protocol) and C − A still misses the +0.0105 gate,
because the insertion cost scales with the appended rows and item 319 already paid −0.0048 for 8004
of them. *[likely. This is why the head route is not the one to run first.]*

## The checkpoint the record refused, and P5–P8 written before it is scored

**Added 17 October, before any score exists.** Item 350 named the remaining form of this axis as a
warm-started multi-head fine-tune over challenge+ChEMBL+qHTS and deferred it on GPU cost. The
organisers have already trained half of it and published it, and this project refused that
checkpoint once for a reason that does not cover the competition.

`openadmet/cyp1a2-cyp2d6-cyp3a4-cyp2c9-chemeleon-v1`: a multitask chemprop D-MPNN warm-started from
CheMeleon, 2048-wide message passing into 1024×3, **four outputs — CYP3A4, CYP2D6, CYP2C9, CYP1A2 —
trained on 8068 ChEMBL 37 structures of pIC50 on exactly our endpoint**, no validation split, to
convergence. `src/embed.py` records the refusal verbatim: *"The CYP-fine-tuned baseline in the same
organisation was NOT used, because a model already trained on CYP pIC50 data would make the
comparison circular."* That is right about the comparison `embed.py` makes — learned representation
against hand descriptors — and says nothing about whether the checkpoint is a member. It is external
public prior information of the exact class item 349 scores 2-for-2, and the training has already
been paid for by someone else.

Two gates passed before the model was run at all. **Overlap with the blinded test set is 0 of 750**
by canonical SMILES of the largest fragment, with the matcher verified against the test set itself;
overlap with our training set is **64, which reproduces item 290's recorded count exactly**. And the
checkpoint's unspecified hyper-parameters are *identified* rather than guessed: a bare state dict
loads at any depth because W_h is shared, so the no-split model's own training labels were used as
the key — depth 6 / relu / mean gives ρ 0.9987 and R² 0.9977 against 0.9329 and R² −25.5 for the
runner-up, a gap of 0.0658 of rank. `src/embed_cypmodel.py` re-checks that fit on every run and
stops below 0.99.

Why this is not items 60–63 again, and it is the whole argument: **the ChEMBL offset is a statement
about LEVEL** (+0.22 to +0.60 against our labels, which is what killed the row merge, because
selection is not correctable by a shift). Spearman is invariant to level exactly, the affine pair
absorbs what is left for ST-RAE, and item 353 says the gap is ordering. The defect that closed the
external-rows axis does not touch the quantity we are behind on.

**P5 — the screen.** The 2048-d frozen embedding, read by `probe_oof` unchanged, scores standalone
macro rank **≥ 0.60** and passes item 348's threshold. *[likely: endpoint-measured class is 2 for 2
over eighteen cases, and this is the only candidate ever screened whose pretraining target IS the
scored endpoint.]*

**P6 — it beats the best member we can build.** Its standalone rank exceeds `surrogate_adme`'s
0.6287 (items 355/360). *[likely but not safe — `surrogate_adme` carries 25 ADME endpoints against
this one's four, and item 324's shipped probe reaches 0.6261 on a checkpoint trained on none of
them, so breadth has been worth as much as aim before.]*

**P7 — the sixth slot.** As a sixth member on the five it is positive at sign 2/2 above the 0.0036
floor. *[likely, and the least informative of the four: item 361 measured that an empty sixth slot is
worth about 0.006 more than a seventh, so almost anything that passes P5 passes this.]*

**P8 — the seventh slot, and this is the prediction that decides whether the axis is open.** As a
**seventh** member beside five + `surrogate_adme` it is positive and clears the 0.0036 floor, and its
error correlation with the composition is **below 0.90** — lower than any candidate on record.
*[genuinely uncertain, and the honest prior is against it. Three pairs have now been measured at
0.907, 0.914 and 0.939 (items 349, 356, 361) and every seventh-slot candidate died. The case for an
exception is that all three were generic encoders — descriptors, 25 ADME endpoints, a bio-trained
representation — while this one is fitted to the four numbers we are scored on, which is a different
information channel rather than a different view of the same one. If P8 fails, the seventh slot is
closed for good and the frozen-encoder axis is finished at one member; if it holds, it is the first
new channel since item 324.]*

**And the part that has no precedent in the record at all:** this checkpoint also emits **four direct
estimates of our four targets**, which no previous candidate did — every one of the twenty-seven was
an encoder with no head on our endpoint. Those four columns are measured in the same device, as a
member and as features, with no pre-registered threshold because nothing comparable exists to set one
from. That half is exploratory and will be reported as such.

## The order to run it in, and what each step costs

1. **Read what is already measured.** `uv run python verify/k129_ncgcread.py --seeds 0,1` nests the
   cached arm over the composition and runs the four controls; it needs no retraining, because the
   arm's predictions are in the cache. Then `src/ablncgc.py --seeds 2,3` for the two missing seeds.
   Hours, not days, and it decides P0.
2. **The arm item 316 costed at one line and nobody ran.** The luciferase counter-screen flag is
   already a column in `data/ncgc/panel.csv`; `grep -n luc src/ablncgc.py` returns nothing. If the
   panel's gain is partly a reporter artefact, dropping AID 411's actives says so, and P0's
   interpretation depends on it.
3. **The censored rows as bounds rather than as a point target.** `src/ablncgc.py` sets
   `panel["target"] = panel.pAC50.fillna(panel.enzyme.map(floor))` — the per-enzyme minimum, a point
   target under squared error. The project owns a one-sided loss already (the dead zone, item 148,
   +0.033 of rank) and on our four enzymes **25170 of the panel's 43881 rows are censored** — 18711
   carry a fitted AC50, the rest are "not active up to the top concentration", which is a bound and
   is currently being fitted as a point. This is the cheapest unexploited information in the source.
4. **Only then the encoder.** Pretrain out of process in the environment items 343 and 349
   established — chemprop's pins never touch the pipeline's, the contract is a file. Start at
   `chemprop_medium`'s shape (d_h 600, depth 4, 0.82 M encoder parameters), not `rdkit2d`'s 2048/6:
   it is the best member we hold, it is ten times cheaper, and item 343's `rand600` control at 0.4942
   already exists for that shape. Embed through `src/embed_medium.py`'s contract, score with
   `verify/k121_chemeleon.py --emb/--tag/--out` (item 345 generalised it for exactly this) and
   `verify/k122_stack.py` (item 349 generalised it to any frozen embedding). Four seeds.

**A baseline hazard for steps 1 and 4.** Item 346's adopted gain — six members + `chemeleon` + the
gbm head, +0.0047 of macro rank and −0.0085 of pair at sign 4/4, "the first adoptable gain since item
324" — is not in METHOD.md, which still documents the composition as five members plus the probe as a
sixth (line 403). Item 348 had to carry that distinction explicitly once already, because item 345's
candidate was a seventh on six, item 324's a sixth on five, and item 347's three a ninth on eight.
Anything measured here must say which stack it was added to, or it is not comparable with the
eighteen-case record.

**The controls are not optional and all of them already exist.** Item 345's alignment control, since
a self-built embedding carries no `molecule_name` either: fit on the array rolled by one row, which
must collapse to ~0. The same-architecture random-weights control (item 343), which separates
pretraining from architecture. The DESC+MECH ridge control (item 324), which separates the encoder
from the member shape. Error correlation against the five members and pairwise against the three
endpoint-class encoders, because P2 is a prediction about exactly that number. The blind-overlap
guard on **largest-fragment InChIKey connectivity blocks**, not whole-record canonical SMILES — item
318 found the 750 test compounds give 748 connectivity blocks, so an enantiomer or salt form of a
test compound passes the committed guard undetected. And for any arm that appends foreign rows to the
feature matrix, item 319's `Ipc` overflow: the guard `xs[xs < 1e-6] = 1.0` cannot catch `inf`, and the
column silently becomes identically zero, which is a third difference hiding inside a two-arm
comparison.

**What runs in a cloud session and what does not, measured rather than assumed.** `uv sync` builds
the pinned environment (pypi is reachable), `git submodule update --init` clones the tutorial, and
`uv run pytest` passes — so the apparatus is fine. What is denied at the proxy is
`huggingface.co`, `pubchem.ncbi.nlm.nih.gov` and `arxiv.org`, which is where the challenge CSVs, the
panel and every checkpoint live. Two consequences, in order of how much they cost:

- **Step 1 needs the four challenge label columns and nothing else.** The arm's own predictions are
  committed; what is missing is `y`, `conf_low` and `conf_high` to score against. No committed
  artefact carries them (`band.json`'s `lo`/`hi` are a bootstrap band, not the compounds' bounds).
- **The six-member nesting additionally needs `data/emb_chemprop_medium*.npz`**, which is not
  committed either and is built from a checkpoint on Hugging Face. Without it the panel arm can
  still be nested over the five committed members in `results/preds/members_seed*.json`, which is a
  real ensemble-level test and the one item 290 applied, but it is not the shipped six.

`verify/k129_ncgcread.py` is written against this and exits 3 with the missing path named. Allowing
`huggingface.co` under the environment's Network access setting is all step 1 needs; the pretraining
in step 4 additionally wants `pubchem.ncbi.nlm.nih.gov` and a GPU or an overnight run on more than
four cores.

## One control that was never run, and the criterion rests on it

`chemprop_chemeleon.pt`'s pretraining target is not recorded anywhere in this tree. Item 349's table
classes it as *измеренная активность*; the provenance on record is that it is the second, larger
checkpoint in the same third entrant's repository as `chemprop_medium.pt` — `jeremycheminf/
openadmet_scripts`, the `adme_pretrain` D-MPNNs (item 324) — which is a reasonable basis for the
label and is not a measurement of it.

The reason it is worth one call anyway: CheMeleon itself is a chemprop D-MPNN pretrained to predict
**Mordred descriptors** for about a million PubChem molecules, and its published shape —
2048-dimensional message passing at depth 6, ~8.7 M encoder parameters — is exactly the shape items
343 and 349 report for this checkpoint. A file named `chemprop_chemeleon.pt` inside an
`adme_pretrain` family is therefore ambiguous on its face: it could be that entrant's activity-trained
model, or the CheMeleon backbone they warm-started from. If it is the latter, one of item 349's three
endpoint successes is a **descriptor** case that helped at +0.0033 sign 4/4, item 350's class table is
wrong where it reads *дескриптор n=5 помогли 0*, and the criterion the whole frozen-checkpoint line
rests on has a counterexample inside its own success column. Settled by reading the checkpoint's own
hyper-parameters — the same one-call control that items 345 and 349 were each caught by, neither of
them by argument.

## The short answer

Not closed prematurely — closed on one half, deferred on the other, and already measured on a third
nobody looked at.

*One trunk, two tasks, the second external* is genuinely finished: built twice, pre-registered, and
closed at item 319 on the channel rather than on the vehicle, which is the closure that does not
invite a rematch. *Pretrain somebody else's panel and transfer the encoder frozen* is the project's
single best external-data result and the only intervention that ever moved the board's ordering; item
350 names doing it on activity data we may use as the one remaining form, deferred on cost, with one
corpus in reach that meets both of its conditions. And the thing to do first is neither: it is to
read `results/preds/oof_ncgc_per.json`, where the panel is worth **+0.0112 of macro rank and −0.0087
of post-pair ST-RAE on two seeds with its shuffle control firing**, carried by the two enzymes item
157 said it hurt. The external-data axis was not closed too early. It was closed on the pooled model
and never reopened on the per-enzyme one, and the number that says so has been in the repository for
a month.
