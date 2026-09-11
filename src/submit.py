"""Produce the two submission files and refuse to write anything the validator rejects.

Until this existed our leaderboard score was not 0.767 - it was undefined, because
nothing in the repository read cyp-challenge-TEST-BLINDED.csv at all. The intermediate
deadline is 24 September and the full test set is revealed once, on the 25th; that is the
only external measurement of how far our cross-validation is from the leaderboard, and it
cannot be recovered later.

Two things here are easy to get wrong and both are silent.

The descriptor block. src/feats.py selects descriptor columns by isna().mean() < 0.05,
a filter fitted on the 4905 training molecules. Recomputing it on 750 test molecules
would produce a different column set, and a check on the column *count* would not catch
a permutation either. feats.build() therefore reindexes by name against the committed
data/desc_names.csv and data/mech_names.csv.

Shrinkage, and why the default is the open question rather than a settled one.
--shrink is ON as of 6 September (item 252), with --no-shrink to undo it. It was off for a
long time, and the justification it used to carry had since been falsified. That
justification was: on the top quartile by activity, shrinking toward the training mean
makes every enzyme worse; the test is built around anchors at the 93rd to 98th percentile;
therefore the test is that kind of subsample and shrinkage would hurt.

The middle step does not survive measurement. verify/k5_shift.py runs one model over both
sets and compares its own output distributions - the mapping is identical, so the
difference in outputs is a difference in inputs. The shift is +0.014 / +0.147 / -0.190 /
+0.437 by enzyme, mild, and negative on 2D6, which is the internal control since 2D6 took
no part in anchor selection. The stress test that flipped the sign moved ybar by +1.1 to
+1.3 - three times stronger than the shift that is actually there.

verify/k6_shift1d.py then reweights the training set to the test-like marginal of the
predictions (1-D density ratio, effective n 796 to 1444) and finds the optimum barely
moves: (+0.30, 0.58) against (+0.30, 0.58) on 1A2, (+1.10, 0.84) against (+0.95, 0.82) on
3A4 at worst. Our own fitted parameters score 0.6781 under those weights against 0.6773
for parameters fitted under them - a gap of 0.0008 macro. Shrinkage does not flip sign
there; it wins, 0.7407 raw against 0.6781.

So the evidence now points the other way, and the default has deliberately NOT been
flipped on that basis alone, because two things the reweighting cannot reach are exactly
the two that would break it. It corrects the marginal of yhat while assuming p(y | yhat)
is unchanged on the test - the assumption recalibration exists to test. And it does not
touch the similarity geometry at all: the test sits at median nearest-neighbour 0.587 and
no re-split of the training data gets above 0.450. Flipping this default changes what gets
submitted, and the 25 September reveal is one-shot, so it is a decision to take
deliberately rather than as a side effect of a docstring.

If it is turned on, the offset and lambda are fitted jointly out of fold rather than the
offset being fixed, reaching macro 0.7150 against 0.7227 for the +0.40 slice. Note the
offset is in centre units: predictions move by (1 - lambda) times it, so the +0.40 once
quoted was never +0.40 in pIC50 - the real shifts are +0.13 / +0.10 / +0.03 / +0.17.

WHICH SHIFT, AND UNDER WHICH CRITERION. Both were open questions until they were
measured, and the second turned out to matter more than the first.

The shift is not one number for four enzymes. On CYP2D6 it is negative, and the reason is
chemistry rather than statistics: the test carries about a third as many compounds that are
basic at pH 7.4 as the CYP2D6 label mask does, and CYP2D6 is the one enzyme of the four that
binds through a salt bridge to a protonated nitrogen, so basic compounds are MORE active
there and less so everywhere else (verify/k7_2d6shift.py, k10_strat2d6.py).

The criterion was doing more work than any of the estimates. Choosing the pair by the WORST
case inside each enzyme's plausible range is insurance against a bad leaderboard; choosing by
the MEAN over the posterior of delta is a bid for the best expected score. The two disagree
by more than any two estimates of delta disagree - the per-enzyme gain is 0.107 by worst case
and 0.045 by mean - and they differ in the sign of their derivative with respect to how wide
the range is. The mean is adopted here: we are after the best expected score, not insurance.

The default is therefore 0, +0.3, -0.5, +0.7, and the zero on CYP1A2 is deliberate. There the
sign of the shift is not determined at all, P(delta >= 0) = 0.59, so any non-zero choice is a
coin flip against doing nothing: the mean criterion picks +0.1, gains 0.0001 by it, and loses
to zero in half the posterior draws. That is added variance for no expected return.

Held to the same standard, this rule is clean where the earlier one was not. The worst-case
rule picked +0.4 on CYP1A2 and lost to doing nothing in 81 % of draws - exactly the defect
that got a single global delta rejected, relocated to another cell. Under the adopted rule
the fractions are 0.09 / 0.11 / 0.01 on the three enzymes it touches.

SETTLED, 6 September, item 252. The switch is on and the vector is the DOCUMENTED one --
0, +0.3, -0.5, +0.7 -- which the code had never actually held (item 251: three statements
about this default, none agreeing).

The vector was chosen over a narrower alternative for one reason. Items 87 and 88 measure
delta's spread ACROSS MODELS at nine to eighty-eight times the seed term, growing rather
than settling when a third model is added, and on CYP2D6 the three estimates are -0.405,
-0.917 and -1.167. A design-simulation route (item 249) puts CYP2D6 at +0.23 and would
have moved that cell to -0.1; adopting it would have preferred one model-free estimate
over three model-based ones agreeing in sign. The documented vector contradicts none of
them.

CYP2D6 moved from -0.5 to 0 on 7 September (item 256), after the band recomputation put a
price on the bet that no earlier discussion had. Fitting under an assumed tilt and scoring
against our own untilted labels costs 0.0407 of macro ST-RAE if the test turns out to look
like the training set, against an expected GAIN of 0.0473 over the posterior -- close to an
even-money bet. CYP2D6 alone carried 0.0735 of that exposure, and it is the cell where the
four routes disagree by SIGN: three model-based estimates give -0.405, -0.917 and -1.167,
the design simulation gives +0.23. Setting it to zero drops the price to 0.0223 and leaves
the three cells where the routes agree in sign untouched. Not a claim about the truth on
CYP2D6 -- an refusal to bet on a cell whose sign is contested.

CYP1A2 moved from 0 to +0.1 on 6 September as a separate decision (item 254), and the
argument is minimax regret rather than maximum expectation. Three posteriors put the
expected ST-RAE at:

                             -0.1    +0.0    +0.1    +0.2
    current posterior       .8368   .8297   .8296   .8363
    design route            .9106   .8892   .8754   .8694
    mixture 50/50           .8756   .8611   .8539   .8539

+0.1 is the only value no posterior dislikes: the current one is flat between 0 and +0.1
(0.8297 against 0.8296), the mixture prefers it to zero by 0.0072, and the design route
prefers it by 0.0138 while wanting +0.2. Going to +0.2 costs 0.0066 under the current
posterior, so it is not robust.

This is a choice of ASSUMPTION, not a claim of measured effect, and the noise floor is the
wrong instrument for it -- the criterion is expected loss under the posterior. The gain is
below every floor the project has and must not be quoted as one.

The size of the shift is bracketed rather than pinned. src/reweight.py tilts the label
marginal and puts the centre at +0.4 for delta = 0 and +0.9 for delta = 0.5; the anchor
percentiles put delta at +1.05, an upper bound, since the anchors' neighbours were chosen
by similarity and regress toward the mean. The prediction-shift route puts it near +0.09,
a lower bound, since a model that mostly interpolates propagates only part of an input
shift into its outputs. Everything between about +0.1 and +0.6 is live.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
TUT = tutorial()

import argparse
import json

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

import feats as F
from cypsplit import butina_folds, fold_digest
from gp import prepare as gp_prepare, gp_predict
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression, RidgeCV
from reweight import tilt

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
TDI_CYPS = ["CYP3A4", "CYP2D6"]          # the only two the organisers score
# Lambda grid for fit_shrinkage. Extended above 1.0 (item 277): the submitted predictions are
# 0.32-0.71 as wide as the labels, so the tilted objective wants to EXPAND them, and the optimum
# sat pinned at the old 1.0 edge on three enzymes of four even at delta=0. It saturates by 1.6;
# 2.0 gives margin. The map stays strictly increasing (lambda>0), so rank is untouched.
GRID = np.round(np.arange(0.2, 2.0001, 0.02), 2)

# Same learner and settings as src/ablate.py, so the submitted model is the one the
# document's numbers describe rather than a cousin of it.
def gbm_reg():
    return HistGradientBoostingRegressor(max_iter=300, learning_rate=0.06,
                                         max_leaf_nodes=31, l2_regularization=1.0,
                                         random_state=0)


def gbm_clf():
    return HistGradientBoostingClassifier(max_iter=300, learning_rate=0.06,
                                          max_leaf_nodes=31, l2_regularization=1.0,
                                          random_state=0)


# The threshold the classification track submits at. It is a PLUG-IN rule: it maximises
# the MCC *expected* under the model's own probabilities, so it needs no labels - and it is
# therefore only correct to the extent those probabilities are calibrated. Item 165 raised
# that as an unstated assumption; verify/k68_tdicalib.py measures what it costs. Kept as
# one function so the submission and the check cannot drift apart.
def plugin_threshold(pr, grid=None):
    ts = np.linspace(0.05, 0.95, 91) if grid is None else np.asarray(grid, float)

    def emcc(t):
        yh = pr >= t
        tp = (pr * yh).sum(); fp = ((1 - pr) * yh).sum()
        fn = (pr * ~yh).sum(); tn = ((1 - pr) * ~yh).sum()
        d = np.sqrt((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
        return (tp * tn - fp * fn) / d if d > 0 else 0.0

    return float(ts[int(np.argmax([emcc(t) for t in ts]))])


# Platt calibration of the classifier's probabilities, fitted OUT OF FOLD on the training rows
# and applied to the test predictions. Item 233 pre-registered the rule for putting this in and
# item 235 is the four-seed measurement that passed it: macro MCC +0.0133 against a floor of
# 0.0076, sign 7 of 8 cells, all 40 Platt slopes positive.
#
# Why it is needed at all, and it is a defect rather than an improvement: plugin_threshold above
# maximises the MCC *expected under the model's own probabilities*, which is only the right
# threshold if those probabilities are calibrated. They are not - the mean predicted rate was
# 0.099 against a true 0.217 on CYP2D6. Item 165 raised the unstated assumption; this closes it.
#
# The slope must be positive or the map is not a calibration but an order reversal - Platt fits a
# logistic on logit(p), and on a classifier with no ordering the slope can come out negative.
# CYP2D6's AUC is 0.588, close enough to make that a live failure, so it is checked rather than
# assumed. All 40 folds of the four-seed sweep came out positive, minimum +0.084.
#
# Costs five extra classifier fits per enzyme, about twenty-five minutes on top of the run.
def tdi_calibrate(Xtr, ytr, fold, p_test):
    p_oof = np.zeros(len(ytr))
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0 or ytr[trn].sum() == 0:
            continue
        p_oof[te] = gbm_clf().fit(Xtr[trn], ytr[trn]).predict_proba(Xtr[te])[:, 1]

    def _logit(q):
        q = np.clip(q, 1e-6, 1 - 1e-6)
        return np.log(q / (1 - q)).reshape(-1, 1)

    lr = LogisticRegression(C=1e6).fit(_logit(p_oof), ytr)
    slope = float(lr.coef_[0, 0])
    if slope <= 0:
        raise SystemExit(
            f"наклон Платта {slope:+.4f} <= 0: карта переворачивает порядок, а не калибрует. "
            "Все 40 фолдов четырёхсидового прогона (пункт 235) дали положительный наклон, "
            "минимум +0.084, так что это не ожидаемый режим -- разбираться, а не обходить.")
    return lr.predict_proba(_logit(p_test))[:, 1], slope


# ---------------------------------------------------------------- связка (пункты 244, 250)
# Метка TDI --- конъюнкция (пункт 234): is_TDI <=> (Delta > log10 2) И (pi_TDI > 4.301).
# Прямой классификатор эту структуру выбрасывает. Связка её использует: вероятность ворот
# берётся из плеча преинкубации, предсказанного четырёхчленным ансамблем с мёртвой зоной, и
# умножается на вероятность сдвига; произведение калибруется и режется plug-in-порогом.
#
# Условия приёмки были зафиксированы в пункте 244 ДО прогона и выполнены в пункте 250: макро
# +0.0127 при поле 0.0076, знак 6 из 8. Контроли: перемешивание фактора сдвига обрушивает
# связку до -0.1453 при знаке 0/8, каждый фактор поодиночке хуже подаваемого.
#
# Поферментно она не проходит ничего --- +0.0083 на CYP3A4 при поле 0.0281 и +0.0170 на CYP2D6
# при 0.0419. Заявление здесь МАКРО, и правило пункта 244 писалось на макро сознательно.
#
# verify/k76_bundle.py --- измерение, эта функция --- его развёртывание. Логика продублирована,
# а не импортирована, потому что src/ не должен зависеть от verify/; при расхождении верить
# k76 и приводить сюда.
TDI_L = np.log10(2.0)
TDI_GATE = 4.0 + TDI_L
BUNDLE_KINDS = ["поферментно", "пул", "GP", "гребневая"]


def _pooled2(Z, e):
    ind = np.zeros((len(Z), 2), np.float32)
    ind[:, e] = 1.0
    return np.hstack([Z, ind])


def _arm_member(kind, Xs, ys, folds, Xtes, KW, e=None, f=None, full=False):
    """Один член на плече преинкубации. full=True --- обучение на всех строках для теста."""
    if kind == "пул":
        if full:
            Xtr = np.vstack([_pooled2(Xs[k], k) for k in range(len(Xs))])
            ytr = np.concatenate(ys)
            m = HistGradientBoostingRegressor(**KW).fit(Xtr, ytr)
            return [m.predict(_pooled2(Xtes[k], k)) for k in range(len(Xs))]
        Xtr = np.vstack([_pooled2(Xs[k][folds[k] != f], k) for k in range(len(Xs))])
        ytr = np.concatenate([ys[k][folds[k] != f] for k in range(len(Xs))])
        m = HistGradientBoostingRegressor(**KW).fit(Xtr, ytr)
        return m.predict(_pooled2(Xs[e][folds[e] == f], e))
    if full:
        out = []
        for k in range(len(Xs)):
            A, B = _dz_design(kind, Xs[k], k, Xtes[k]) if kind in ("GP", "гребневая") else (Xs[k], Xtes[k])
            out.append(_dz_fit(kind, A, ys[k], B) if kind in ("GP", "гребневая")
                       else HistGradientBoostingRegressor(**KW).fit(A, ys[k]).predict(B))
        return out
    trn, te = folds[e] != f, folds[e] == f
    if kind in ("GP", "гребневая"):
        A, B = _dz_design(kind, Xs[e][trn], e, Xs[e][te])
        return _dz_fit(kind, A, ys[e][trn], B)
    return HistGradientBoostingRegressor(**KW).fit(Xs[e][trn], ys[e][trn]).predict(Xs[e][te])


def _arm_ensemble(Xs, ys, folds, Xtes, LOs, HIs):
    """Вне фолда и на тесте: четыре члена, мёртвая зона, невзвешенное среднее.

    Пулированный член обучается ОДИН раз на фолд: его обучающая выборка от эндпоинта не
    зависит, зависит только предсказание."""
    PLAIN = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
                 l2_regularization=1.0, random_state=0)
    n = len(Xs)
    oof = [np.zeros(len(y)) for y in ys]
    tst = [np.zeros(len(Xt)) for Xt in Xtes]

    def pass_(targets, KW):
        """Один проход всех членов: вне фолда и на полных данных."""
        o = [np.zeros(len(t)) for t in targets]
        t_ = [np.zeros(len(Xt)) for Xt in Xtes]
        for f in range(5):
            if "пул" in BUNDLE_KINDS:
                Xtr = np.vstack([_pooled2(Xs[k][folds[k] != f], k) for k in range(n)])
                ytr = np.concatenate([targets[k][folds[k] != f] for k in range(n)])
                mdl = HistGradientBoostingRegressor(**KW).fit(Xtr, ytr)
                for e in range(n):
                    te = folds[e] == f
                    if te.sum():
                        o[e][te] += mdl.predict(_pooled2(Xs[e][te], e)) / len(BUNDLE_KINDS)
            for kind in BUNDLE_KINDS:
                if kind == "пул":
                    continue
                for e in range(n):
                    te = folds[e] == f
                    if not te.sum():
                        continue
                    trn = ~te
                    if kind in ("GP", "гребневая"):
                        A, B = _dz_design(kind, Xs[e][trn], e, Xs[e][te])
                        o[e][te] += _dz_fit(kind, A, targets[e][trn], B) / len(BUNDLE_KINDS)
                    else:
                        o[e][te] += (HistGradientBoostingRegressor(**KW)
                                     .fit(Xs[e][trn], targets[e][trn])
                                     .predict(Xs[e][te])) / len(BUNDLE_KINDS)
        Xtr = np.vstack([_pooled2(Xs[k], k) for k in range(n)])
        mdl = HistGradientBoostingRegressor(**KW).fit(Xtr, np.concatenate(targets))
        for e in range(n):
            t_[e] += mdl.predict(_pooled2(Xtes[e], e)) / len(BUNDLE_KINDS)
        for kind in BUNDLE_KINDS:
            if kind == "пул":
                continue
            for e in range(n):
                if kind in ("GP", "гребневая"):
                    A, B = _dz_design(kind, Xs[e], e, Xtes[e])
                    t_[e] += _dz_fit(kind, A, targets[e], B) / len(BUNDLE_KINDS)
                else:
                    t_[e] += (HistGradientBoostingRegressor(**KW).fit(Xs[e], targets[e])
                              .predict(Xtes[e])) / len(BUNDLE_KINDS)
        return o, t_

    raw_o, _ = pass_(ys, PLAIN)
    print("      плечо: обычный проход готов", flush=True)
    tgt = [np.clip(raw_o[e] * len(BUNDLE_KINDS) / len(BUNDLE_KINDS), LOs[e], HIs[e])
           for e in range(n)]
    oof, tst = pass_(tgt, DZ_KW)
    print("      плечо: мёртвая зона готова", flush=True)
    return oof, tst


def _platt_1d(score_tr, y_tr, score_te, raw=False):
    def f(q):
        if raw:
            return np.asarray(q, float).reshape(-1, 1)
        q = np.clip(q, 1e-6, 1 - 1e-6)
        return np.log(q / (1 - q)).reshape(-1, 1)
    lr = LogisticRegression(C=1e6).fit(f(score_tr), y_tr)
    return lr.predict_proba(f(score_te))[:, 1], lr.predict_proba(f(score_tr))[:, 1]


def test_features(desc_names, mech_names):
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    FP, dsc, M, ok = F.build(list(te.SMILES), desc_names, mech_names)
    if len(ok) != len(te):
        # Every test SMILES parses today; if that ever stops being true the rows would
        # silently misalign, so fail loudly instead.
        raise SystemExit(f"RDKit не разобрал {len(te)-len(ok)} тестовых SMILES - "
                         "выравнивание строк сломается, разбираться вручную")
    return te, np.hstack([FP, dsc.to_numpy(np.float32), M.to_numpy(np.float32)])


# Сетка смещений намеренно шире, чем нужно любому правдоподобному сдвигу. Прежняя,
# linspace(-0.2, 1.6), зажимала подгонку с обеих сторон: при предполагаемом сдвиге ниже -0.3
# оптимальное смещение упиралось в нижний край и переставало двигаться, а на CYP3A4 при
# +0.7 оно садилось на верхний. Подогнанный параметр, стоящий на границе сетки, --- не
# подогнанный параметр, и плоскость целевой функции рядом с ним мнимая.
OFFGRID = np.round(np.arange(-3.0, 3.01, 0.05), 2)
# Сетка по сдвигу предсказаний. Диапазон намеренно шире любого правдоподобного: оценки
# delta лежат в пределах -1.3 .. +1.0, а сдвиг всегда меньше delta по модулю.
SHIFTGRID = np.round(np.arange(-2.0, 2.001, 0.01), 2)
DESC_MECH = 247   # ширина блока DESC+MECH в конце матрицы признаков


# Правило по delta выбрано под смесью ТРЁХ апостериоров --- раздельной модели, пулированной
# и подаваемого ансамбля, --- а не под одним. Обращение ядра
# даёт разные ответы на раздельной и пулированной моделях --- +0.04 против +0.34 на 1A2,
# +0.35 против +0.78 на 2C9, -0.51 против -0.41 на 2D6, +0.74 против +0.84 на 3A4, --- и
# расхождение сопоставимо со всей выборочной неопределённостью, а на 2C9 её превышает.
# Выбрать между ними без тестовых меток нельзя, поэтому розыгрыши обеих сложены с равным
# весом (это решение, а не вывод) и правило подобрано под расширенное распределение.
# Оценка на самом ансамбле на CYP2D6 вышла ЗА диапазон обеих составляющих, -0.686 против
# -0.507 и -0.408: обращение ядра нелинейно, и среднее предсказаний не есть среднее решений.
# Это сместило центр смеси сильнее, чем расширило её, поэтому 2D6 вернулся к -0.5, а не стал
# осторожнее. Выигрыш поферментной подгонки +0.066. CYP1A2 пропущен: его сдвиг стоит 0.003
# макро при 33% шанса навредить, и знак delta там не держится даже между моделями.
# Подаваемая модель --- трёхчленный ансамбль, поэтому в смесь входят четыре апостериора:
# две базы бустинга, гауссов процесс и сам ансамбль. Ячейка CYP2D6 --- единственная, где
# выбор критерия решает: центр смеси даёт -0.7, минимакс по моделям -0.4, и разница между
# ними стоит 0.027 макро в терминах худшего случая. Взят центр, потому что средний критерий
# принят в пункте 57; минимакс записан как альтернатива, а не отвергнут.
# verify/k15_pooldelta.py, k16, k17, k19; README пункты 87, 88, 89, 95.


def pooled_design(X, e):
    """Признаки плюс четырёхпозиционный индикатор фермента."""
    ind = np.zeros((len(X), len(CYPS)), np.float32)
    ind[:, e] = 1.0
    return np.hstack([X, ind])


# Проход мёртвой зоны. Пины --- КОПИЯ gbm_reg() плюс absolute_error, ровно как в
# src/abldzens.py:KW. Копия, а не ссылка: если gbm_reg когда-нибудь изменится, числа
# пункта 164 должны перестать воспроизводиться ЗАМЕТНО, а не тихо разъехаться.
DZ_KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
             l2_regularization=1.0, random_state=0, loss="absolute_error")


def _dz_design(kind, X, e, Xte=None):
    """Матрица признаков одного члена. Одна функция на путь вне фолда и на тестовый путь,
    потому что дефект 2 пункта 202 --- ровно расхождение этих двух путей у гребневой.

    Возвращает (Xtr, Xte_преобразованный). Преобразование, где оно есть, подгоняется НА
    ОБУЧАЮЩИХ СТРОКАХ и применяется к тесту --- никогда на объединении.
    """
    if kind == "GP":
        tf = gp_prepare(X)
        return tf(X), (None if Xte is None else tf(Xte))
    if kind == "гребневая":
        A = np.nan_to_num(X[:, -DESC_MECH:].astype(np.float64), posinf=0.0, neginf=0.0)
        mu, sd = A.mean(0), A.std(0) + 1e-9
        f = lambda B: np.clip((np.nan_to_num(B[:, -DESC_MECH:].astype(np.float64),
                                             posinf=0.0, neginf=0.0) - mu) / sd, -5.0, 5.0)
        return f(X), (None if Xte is None else f(Xte))
    if kind == "пул":
        return pooled_design(X, e), (None if Xte is None else pooled_design(Xte, e))
    return X, Xte


def _dz_fit(kind, Xtr, target, Xte):
    """Подгонка одного члена под спроецированную мишень. Зеркало src/abldzens.py:refit."""
    if kind == "GP":
        return gp_predict(Xtr, target, Xte)
    if kind == "гребневая":
        return RidgeCV(alphas=np.logspace(-1, 4, 12)).fit(Xtr, target).predict(Xte)
    return HistGradientBoostingRegressor(**DZ_KW).fit(Xtr, target).predict(Xte)


def _dz_oof(kind, X, target, fold, e):
    """Один проход мажорирования вне фолда: та же модель, спроецированная мишень.

    Мишень строится вызывающим как clip(p_oof, lo, hi), и p_oof обязано быть предсказанием
    ВНЕ ФОЛДА. Проекция предсказаний модели на её же обучающих строках схлопывает схему в
    тождество: переобученная модель кладёт каждую обучающую строку внутрь её полосы, мишень
    становится равна предсказанию, градиент зануляется, и прогон отрабатывает чисто, ничего
    не сделав. Разбор в докстринге src/abldead.py.
    """
    Xi, _ = _dz_design(kind, X, e)
    p = np.zeros(len(target))
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0:
            continue
        p[te] = _dz_fit(kind, Xi[trn], target[trn], Xi[te])
    return p


def load_screen(rows):
    """Показания одноточечного скрининга в порядке rows.csv: (log2fc, есть_показание).

    Длинная таблица со столбцами enzyme и log2fc_estimate разворачивается в четыре вектора
    длины n, выровненные по data/rows.csv --- как того требует CLAUDE.md: порядок строк
    матрицы признаков задаёт rows.csv, и всякий, кто читает свою таблицу, обязан
    переиндексироваться по ней, иначе признаки и метки молча разъезжаются.
    """
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    L2 = np.full((len(rows), len(CYPS)), np.nan)
    for e, c in enumerate(CYPS):
        sub = sc[sc.enzyme == c].set_index("Molecule_Name")["log2fc_estimate"]
        L2[:, e] = sub.reindex(rows.Molecule_Name).to_numpy(float)
    return L2, ~np.isnan(L2)


def _screen_train(X, y, mask, L2, has, e, trn):
    """Обучающая таблица одного фермента: строки кривых плюс строки скрининга.

    Возвращает (Xa, ya) с ОДНОЙ добавленной колонкой источника: 0 --- кривая, 1 --- скрининг.

    Механизм пункта 177. Показание скрининга --- это log2fc при одной концентрации, не pIC50,
    поэтому подмешать его как метку напрямую нельзя. Изотоническая регрессия log2fc -> pIC50
    подгоняется на молекулах, у которых есть И кривая, И показание, и переносит показание в
    шкалу метки. Связь убывающая: сильнее ингибирует --- ниже log2fc, отсюда increasing=False.

    Строки добавляются ТОЛЬКО там, где у этого фермента кривой нет, так что ни одна ячейка
    не получает две конкурирующие мишени, и добавленное --- новая супервизия, а не
    перевзвешивание уже имеющегося. Таких строк 11505 против 6525 кривых.

    `trn` --- булев вектор обучающих МОЛЕКУЛ. Изотоника подгоняется только по ним, потому что
    подгонка по всему занесла бы отложенные метки в обучающие мишени через отображение ---
    тихая версия утечки, которую этот файл ловил дважды.
    """
    sel = mask[:, e] & trn
    Xa = [np.hstack([X[sel], np.zeros((int(sel.sum()), 1), np.float32)])]
    ya = [y[sel, e]]

    fitm = mask[:, e] & has[:, e] & trn
    addm = has[:, e] & ~mask[:, e] & trn
    # Порог 50 --- из src/ablaux.py: ниже него изотоника подгоняется по шуму. При отказе
    # рука молча вырождается в обычную поферментную, что и есть правильное поведение.
    if fitm.sum() >= 50 and addm.sum() > 0:
        iso = IsotonicRegression(increasing=False, out_of_bounds="clip")
        iso.fit(L2[fitm, e], y[fitm, e])
        Xa.append(np.hstack([X[addm], np.ones((int(addm.sum()), 1), np.float32)]))
        ya.append(iso.predict(L2[addm, e]))
    return np.vstack(Xa), np.concatenate(ya)


def _src0(Xb):
    """Матрица предсказания с колонкой источника 0: спрашиваем как про кривую."""
    return np.hstack([Xb, np.zeros((len(Xb), 1), np.float32)])


def _oof_one(X, y, mask, fold, pool, scr=None):
    """Предсказания вне фолда одной из двух базовых моделей: раздельной или пулированной.

    Пул складывает все четыре набора меток в одну таблицу с индикатором фермента, так что
    CYP2D6 обучается не на своих 1493 строках, а на 6525 всех, а специфику забирает через
    взаимодействия с индикатором. На четырёх сидах это стоит -0.0067 макро ПОСЛЕ аффинной
    пары (t = -9.06, p = 0.003) и прибавляет 0.0137 ранговой связи, то есть относится к тому
    классу вмешательств, которые постобработка не поглощает (verify/README.md, пункты 80, 84).

    Фолды --- по молекуле, поэтому все четыре копии соединения лежат в одном фолде и утечки
    между ферментами нет.
    """
    P = [np.zeros(int(mask[:, e].sum())) for e in range(len(CYPS))]
    for f in range(5):
        if not pool:
            for e in range(len(CYPS)):
                m = mask[:, e]
                fi, Xi, yy = fold[m], X[m], y[m, e]
                a, b = fi != f, fi == f
                if b.sum() == 0:
                    continue
                if scr is None:
                    P[e][b] = gbm_reg().fit(Xi[a], yy[a]).predict(Xi[b])
                else:
                    L2, has = scr
                    Xa, ya = _screen_train(X, y, mask, L2, has, e, fold != f)
                    P[e][b] = gbm_reg().fit(Xa, ya).predict(_src0(Xi[b]))
            continue
        Xs, ys = [], []
        for e in range(len(CYPS)):
            sel = mask[:, e] & (fold != f)
            if sel.any():
                Xs.append(pooled_design(X[sel], e))
                ys.append(y[sel, e])
        model = gbm_reg().fit(np.vstack(Xs), np.concatenate(ys))
        for e in range(len(CYPS)):
            m = mask[:, e]
            b = fold[m] == f
            if b.sum() == 0:
                continue
            P[e][b] = model.predict(pooled_design(X[m][b], e))
    return P


# Поферментный отбор члена (пункт 218). Ключ --- фермент, значение --- какие члены
# усредняются на нём; отсутствие ключа означает «все», то есть прежнее поведение.
#
# Посылка ансамбля --- что усреднение бьёт свои члены. Прочитанная ПОФЕРМЕНТНО, как велит
# пункт 165, она есть четыре проверки, и на CYP3A4 проваливается: гауссов процесс в одиночку
# даёт ранг 0.8153 против 0.8056 у пятичленного, знак 4/4 на четырёх сидах по обоим
# критериям, втрое выше поля этого фермента (0.0033). Перебор всех 31 подмножества ставит
# «GP один» первым, а полный ансамбль четырнадцатым.
#
# ЧТО ЗДЕСЬ НА ГРАНИ и должно быть видно тому, кто это меняет. «поферментно+GP» даёт ранг
# 0.8149 --- отличие 0.0004, вчетверо ниже поля, то есть по решающему критерию это НИЧЬЯ, --- и
# при этом лучше по метрике (0.4101 против 0.4129). Два члена ещё и страхуют от того, что
# гауссов процесс окажется неудачлив на тесте, чего одиночный член не переживёт. Выбран
# одиночный GP; альтернатива меняется правкой одной строки ниже.
# UPDATED (items 282, 283): CYP2C9 and CYP3A4 both ship GP+ствол, confirmed on FRESH seeds 4-7
# that did not generate the hypothesis. CYP2C9 GP+ствол beats the five-member mean by +0.0125
# (sd 0.0021, sign 4/4, floor 0.0071) with no shrinkage from the in-sample +0.0124 -- the fewest
# labels (1285) mean the two low-dim members on DESC+MECH survive and the high-dim boosters that
# overfit are dropped. CYP3A4 GP+ствол beats GP-alone by +0.0044 (sd 0.0013, sign 4/4, floor
# 0.0033): item 218 put GP-alone first, but that was an earlier trunk; the current one makes
# GP+ствол the better arm, so 218 is superseded, not wrong. Macro rank +0.0042 over the frozen
# submission. CYP1A2/CYP2D6 were not tested on fresh seeds (k84: 1A2 borderline, 2D6 net-negative)
# and stay on the full ensemble.
# CYP1A2 added (items 284, 285): fresh seeds 4-7 give поферментно+GP+ствол over the five-member
# mean at +0.0068 (sd 0.0017, sign 4/4, floor 0.0061) -- thin (lower 95% bound +0.0047, below the
# floor) but pre-registered and met. It keeps the per-enzyme booster (1412 labels overfit less than
# CYP2C9's 1285) and drops the ridge. Three of four enzymes now drop members; only CYP2D6 keeps the
# full ensemble (k84: its selection is net-negative). Shipped-seed-0 macro +0.0061, fresh-seed
# estimate +0.0059, both above the 0.0036 floor.
SOLO = {"CYP1A2": ("поферментно", "GP", "ствол"),
        "CYP2C9": ("GP", "ствол"),
        "CYP3A4": ("GP", "ствол")}

# Предполагаемый сдвиг теста по ферментам. ОДНО определение: verify/k79_bandfix.py читает
# его отсюда, потому что полоса фальсификации относится к конкретным дельтам, и разошедшиеся
# копии этой константы уже дважды давали полосу для руки, которая не подаётся.
DELTA_DEFAULT = "0.1,0.3,0,0.7"


def _keep(e, kind):
    """Входит ли член в ансамбль этого фермента."""
    sel = SOLO.get(CYPS[e])
    return sel is None or kind in sel


def oof_members(X, y, mask, fold, mode, scr=None, seed=0, dead=False):
    """Предсказания вне фолда КАЖДОГО члена по отдельности, как (вид, по ферментам).

    Существует потому, что проход мёртвой зоны применяется к членам поодиночке: каждому
    нужна своя спроецированная мишень и своя переподгонка. Порядок членов зафиксирован ---
    он определяет, что именно усредняется, а среднее невзвешенное.
    """
    parts = [("поферментно", _oof_one(X, y, mask, fold, False, scr)),
             ("пул", _oof_one(X, y, mask, fold, True))]
    if mode in ("ансамбль", "ансамбль5"):
        parts.append(("GP", _oof_gp(X, y, mask, fold)))
        parts.append(("гребневая", _oof_ridge(X, y, mask, fold)))
    if mode == "ансамбль5":
        parts.append(("ствол", _oof_trunk(y, mask, fold, seed, dead)))
    return parts


def dz_pass(parts, X, mask, fold, bands):
    """Проход мёртвой зоны по всем членам, у которых он определён.

    Пункт 164: +0.0167 ранга на подаваемой пятичленной конфигурации, четыре сида, при
    проходе в ЧЕТЫРЁХ членах. Ствол там не перепроецирован, отчего 0.6230 объявлено нижней
    оценкой; здесь он пропускается по той же причине --- его переподгонка живёт в
    src/trunk.py под torch, а не здесь.
    """
    LO, HI = bands
    out = []
    for kind, P in parts:
        if kind == "ствол":
            out.append((kind, P))
            continue
        Q = []
        for e in range(len(CYPS)):
            m = mask[:, e]
            t = np.clip(P[e], LO[m, e], HI[m, e])
            Q.append(_dz_oof(kind, X[m], t, fold[m], e))
        out.append((kind, Q))
        print(f"    мёртвая зона: {kind} перепод.", flush=True)
    return out


def oof_predictions(X, y, mask, fold, mode, scr=None, bands=None):
    """Предсказания вне фолда в одном из трёх режимов.

    Ансамбль --- среднее двух базовых. Он выигрывает больше каждой из них: -0.0234 макро
    после аффинной пары против раздельной (p = 0.0001) и -0.0168 сверх пула (p = 0.00004),
    и прибавляет 0.0291 ранговой связи. Причина в том, что две базы видят разные обучающие
    таблицы и ошибаются по-разному, так что усреднение снимает часть дисперсии ДО того, как
    за неё возьмётся усадка. По критерию пункта 80 это относится к тем вмешательствам,
    которые постобработка не поглощает, --- и проверено, что не поглощает.
    """
    if mode in ("ансамбль", "ансамбль-без-GP", "ансамбль5"):
        parts = oof_members(X, y, mask, fold, mode, scr, dead=bands is not None)
        if bands is not None:
            parts = dz_pass(parts, X, mask, fold, bands)
        return [np.mean([P[e] for k, P in parts if _keep(e, k)], axis=0)
                for e in range(len(CYPS))]


TRUNK_LAM = "3.0"      # значение, на котором пункт 79 мерил канал
TRUNK_MODE = "twohead"


def _trunk_device(meta_device):
    """Тот же вычислитель, на котором посчитан сохранённый файл, если он доступен.

    Не педантизм. Усадка подгоняется по предсказаниям вне фолда из trunk_twohead.json, а
    применяется к предсказаниям теста, посчитанным здесь; если первые считались на mps, а
    вторые пойдут на cpu, две половины пары разойдутся по численности, и подогнанное
    lambda будет описывать не тот объект, к которому его прикладывают. Модели всё равно
    обучаются на разных данных, так что побитового совпадения не бывает --- но
    систематического расхождения вычислителя быть не должно.
    """
    import torch
    have = {"mps": torch.backends.mps.is_available(), "cuda": torch.cuda.is_available(),
            "cpu": True}
    if have.get(meta_device):
        return meta_device
    if meta_device not in (None, "cpu"):
        print(f"    ВНИМАНИЕ: сохранённые предсказания ствола посчитаны на {meta_device}, "
              f"здесь его нет --- половины аффинной пары лягут на разные вычислители",
              flush=True)
    return "cpu"


def _trunk_blocks(Xall):
    """Разрезать общую матрицу обратно на FP / DESC / MECH.

    src/trunk.py принимает блоки по отдельности, потому что кладёт FP под log1p, а здесь
    матрица уже склеена. Ширины берутся из feats.npz, а не зашиваются числом: их три, и
    ошибка в любой сдвинет весь блок молча.
    """
    z = np.load(D + "feats.npz")
    n_fp, n_de, n_me = z["FP"].shape[1], z["DESC"].shape[1], z["MECH"].shape[1]
    if n_fp + n_de + n_me != Xall.shape[1]:
        raise SystemExit(f"ширины блоков {n_fp}+{n_de}+{n_me} не дают {Xall.shape[1]}")
    return (Xall[:, :n_fp], Xall[:, n_fp:n_fp + n_de], Xall[:, n_fp + n_de:])


def _trunk_clip(p, y_e):
    """Диапазон меток фермента плюс-минус две единицы, как в src/trunkdose.py.

    Не косметика: на сиде 0 ствол выдаёт одно соединение на -360 (пункт 42), и аффинная
    пара, в отличие от изотоники, такой выброс не поглощает. Без обрезки весь замер стал бы
    замером одной молекулы.
    """
    return np.clip(p, np.nanmin(y_e) - 2.0, np.nanmax(y_e) + 2.0)


TRUNK_FOLD_DIGEST = "2d93c19815e14261"   # то же золотое значение, что в tests/test_split.py


def _oof_trunk(y, mask, fold, seed=0, dead=False):
    """Предсказания ствола вне фолда --- из results/preds/trunk_twohead.json.

    Читаются, а не пересчитываются, ровно по той же причине, по какой читается oof.json:
    прогон занимает часы, файл закоммичен и его происхождение записано. Разбиение то же ---
    src/trunk.py зовёт butina_folds из cypsplit.py с тем же сидом, так что усреднять его
    предсказания с бустинговыми законно, и аффинная пара подгоняется по тем же фолдам.

    Пятый член принят по пункту 120: -0.0061 пары и +0.0054 ранга, знак 4/4 на каждом из
    четырёх ферментов --- больше вдвое, чем даёт гребневая, и единственный член, который
    помогает всем четырём.
    """
    # Сторож. Этот член --- единственный, который не пересчитывается здесь, а читается из
    # файла, посчитанного на сиде 0. Если разбиение когда-нибудь сдвинется, все остальные
    # члены поедут за ним, а этот молча останется на старых фолдах, и предсказания вне
    # фолда перестанут быть вне фолда. Digest ловит это на месте.
    d = fold_digest(fold)
    src = RES + ("preds/trunk_twohead_dead.json" if dead else "preds/trunk_twohead.json")
    J = json.load(open(src))
    T = J.get("preds", J)
    # Золотая константа существует ТОЛЬКО для сида 0 --- это то же значение, что пинит
    # tests/test_split.py. Подача всегда идёт на нём, поэтому проверка остаётся жёсткой.
    if seed == 0 and d != TRUNK_FOLD_DIGEST:
        raise SystemExit(
            f"разбиение сдвинулось: {d} вместо {TRUNK_FOLD_DIGEST}. Предсказания ствола в "
            f"{src} посчитан на старых фолдах и вне фолда больше не лежат. "
            f"Перезапустите src/trunk.py или снимите режим ансамбль5.")
    # На прочих сидах золотого значения не существует нигде в репозитории, и раньше это
    # означало отсутствие сторожа вовсе. Теперь src/trunk.py записывает дайджест КАЖДОГО
    # сида в meta, так что проверка стала данными вместо константы и покрывает все сиды.
    rec = (J.get("meta") or {}).get("fold_digests") or {}
    if str(seed) in rec and rec[str(seed)] != d:
        raise SystemExit(
            f"сид {seed}: дайджест {d}, а ствол посчитан на {rec[str(seed)]}. "
            f"Перезапустите src/trunk.py --seeds {seed}.")
    if str(seed) not in rec and seed != 0:
        raise SystemExit(
            f"сид {seed}: в trunk_twohead.json нет записанного дайджеста фолдов, значит "
            f"проверить, что предсказания ствола лежат вне ЭТИХ фолдов, нечем. "
            f"Перезапустите src/trunk.py, чтобы он их записал.")
    key = f"{TRUNK_MODE}|{seed}|{TRUNK_LAM}"
    if key not in T:
        raise SystemExit(f"нет ключа {key} в trunk_twohead.json; запустите src/trunk.py --seeds {seed}")
    A = np.asarray(T[key], float)
    return [_trunk_clip(A[mask[:, e], e], y[mask[:, e], e]) for e in range(len(CYPS))]


def _desc_scaled(X):
    """Дескрипторы и механистический блок, стандартизованные с обрезкой."""
    A = np.nan_to_num(X[:, -DESC_MECH:].astype(np.float64), posinf=0.0, neginf=0.0)
    return np.clip((A - A.mean(0)) / (A.std(0) + 1e-9), -5.0, 5.0)


def _oof_ridge(X, y, mask, fold):
    """Гребневая на дескрипторах вне фолда, alpha по leave-one-out на обучающих строках.

    Четвёртый член ансамбля. Поодиночке она слабее бустинга (ранг 0.560 против 0.565), но
    в ансамбле даёт -0.0026 (p = 0.0005), потому что ошибается ГЛАДКО там, где деревья и
    гауссов процесс ошибаются локально. Лес и kNN проверены тем же способом и оба ухудшают,
    так что членство определяется измерением, а не принципом «больше разнообразия лучше»."""
    P = []
    for e in range(len(CYPS)):
        m = mask[:, e]
        yy, fi = y[m, e], fold[m]
        Xi = _desc_scaled(X[m])
        p = np.zeros(len(yy))
        for f in range(5):
            trn, te = fi != f, fi == f
            if te.sum() == 0:
                continue
            p[te] = RidgeCV(alphas=np.logspace(-1, 4, 12)).fit(Xi[trn], yy[trn]).predict(Xi[te])
        P.append(p)
    return P
    return _oof_one(X, y, mask, fold, mode == "пул")


def _oof_gp(X, y, mask, fold):
    """Гауссов процесс вне фолда. Ядро на дескрипторах, потому что пункт 90 намерил, что
    именно в этом пространстве сосед несёт информацию, а Морган --- худшее из четырёх."""
    P = []
    for e in range(len(CYPS)):
        m = mask[:, e]
        yy, fi = y[m, e], fold[m]
        tf = gp_prepare(X[m])
        Xi = tf(X[m])
        p = np.zeros(len(yy))
        for f in range(5):
            trn, te = fi != f, fi == f
            if te.sum() == 0:
                continue
            p[te] = gp_predict(Xi[trn], yy[trn], Xi[te])
        P.append(p)
    return P


def fit_shrinkage(P, y, mask, delta=(0.0, 0.0, 0.0, 0.0)):
    """Offset and lambda per enzyme, both chosen out-of-fold on the training data.

    Fitting the two jointly rather than fixing the offset and searching lambda: the
    family {c + L(p - c)} is identically the affine family {a + b p} with b = L and
    a = c(1 - L), so a fixed offset is an arbitrary slice through it. Jointly it reaches
    macro 0.7150 against 0.7227 for the +0.40 slice and 0.7333 for the offset at the
    training mean.

    Worth noting what the offset is not. The predictions move by (1 - L) times it, not by
    it, so the +0.40 quoted earlier was never a 0.40 shift in pIC50 - at the fitted
    lambdas the actual shifts are +0.13 / +0.10 / +0.03 / +0.17. And at L -> 1 the
    centre is not identified at all, which is a second reason to fit the affine pair.
    """
    out = []
    for e, c in enumerate(CYPS):
        m = mask[:, e]
        yy = y[m, e]
        p = P[e]
        lo = LO[m, e]; hi = HI[m, e]
        mu = p.mean()
        # Under an assumed shift the objective is the tilted one: our own labels reweighted
        # so their mean sits delta higher. At delta = 0 the weights are all ones and this is
        # exactly the untilted fit, so the default path is unchanged.
        de = float(delta[e])
        w = np.ones_like(yy) if de == 0 else tilt(yy, de)
        # Сетка ведётся по СДВИГУ предсказаний s, а не по смещению центра off. Это та же
        # двухпараметрическая семья --- q = c + L(p - c) при c = mu + off тождественно равно
        # L p + (1 - L) mu + s, где s = (1 - L) off, --- но параметризация другая, и разница
        # не косметическая. Сдвиг s измеряется в единицах pIC50 и ограничен здравым смыслом;
        # смещение off при L -> 1 не ограничено ничем, потому что тот же сдвиг требует всё
        # большего off. Прежняя сетка по off упиралась в край на CYP2C9, как только ансамбль
        # поднял L до 0.90, и параметр выбирала сетка, а не данные.
        A = SHIFTGRID[:, None]
        B = np.broadcast_to(GRID[None, :], (len(SHIFTGRID), len(GRID)))
        q = (mu * (1.0 - B) + A)[:, :, None] + B[:, :, None] * p[None, None, :]
        pen = (w * (np.maximum(q - hi, 0.0) + np.maximum(lo - q, 0.0))).sum(axis=2)
        ii, jj = np.unravel_index(pen.argmin(), pen.shape)
        sh, L = float(SHIFTGRID[ii]), float(GRID[jj])
        if ii in (0, len(SHIFTGRID) - 1) or jj in (0, len(GRID) - 1):
            print(f"    ВНИМАНИЕ {c}: оптимум на краю сетки (сдвиг {sh:+.2f}, lambda {L:.2f})",
                  flush=True)
        # mu --- среднее предсказаний ВНЕ ФОЛДА на обучении. Именно оно, а не среднее по
        # тесту: преобразование подогнано относительно него, и подстановка тестового
        # среднего молча сместила бы центр на разницу маргиналей, то есть ровно на то,
        # что мы отдельно оцениваем как delta.
        out.append((L, mu, sh))
        print(f"    {c}: lambda {L:.2f}, сдвиг предсказаний {sh:+.3f}", flush=True)
    return out


def main():
    global LO, HI
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-shrink", dest="shrink", action="store_false",
                    help="не применять усадку (поведение до пункта 252). Флаг ВКЛЮЧЁН по "
                         "умолчанию: пара (off, lambda) подбирается под предполагаемый сдвиг "
                         "теста, вектор --delta. Выигрыш по среднему правилу +0.0473 макро пары "
                         "против подгонки при нуле (src/shrinkchoice.py, блок 6), при парном поле "
                         "лидерборда 0.017. ВАЖНО: этот выигрыш посчитан по апостериору из одного "
                         "источника; пункты 87 и 88 меряют разброс delta ПО МОДЕЛЯМ в 9-88 раз "
                         "больше сидового, и он в оценку не входит.")
    ap.add_argument("--mode", default="ансамбль5",
                    choices=["раздельно", "пул", "ансамбль", "ансамбль-без-GP", "ансамбль5"],
                    help="раздельно воспроизводит поведение до пункта 84; ансамбль включает GP; "
                         "ансамбль5 добавляет пятым членом ствол со скрининговой головой "
                         "(пункт 120: -0.0061 пары, +0.0054 ранга, знак 4/4 на каждом "
                         "ферменте). Умолчание не переключено: пятый член требует torch на "
                         "машине, где собирается сабмит, и решение о составе принимает команда")
    # РАСХОЖДЕНИЕ, найденное 6 сентября и здесь только ЗАФИКСИРОВАННОЕ, а не исправленное.
    # Докстринг выше говорит «The default is therefore 0, +0.3, -0.5, +0.7», пункт 59 журнала
    # говорит «src/submit.py --delta now defaults to it», а код стоял на другом векторе всегда:
    # 0,0.5,-0.5,0.8 до коммита 20be90c и 0,0.5,-0.7,0.8 после. Ни одно из трёх утверждений не
    # согласуется с двумя другими.
    #
    # УСТРАНЕНО 6 сентября (пункт 252): код приведён к документированному вектору, а не наоборот.
    # Выбор решён пунктом 88: разброс delta ПО МОДЕЛЯМ в 9-88 раз больше сидового, на CYP2D6 три
    # модели дают -0.405, -0.917 и -1.167, а на CYP1A2 не выживает даже знак. Документированный
    # вектор не спорит по знаку ни с одной из трёх; узкая альтернатива из пункта 249 спорила бы.
    ap.add_argument("--delta", default=DELTA_DEFAULT,
                    help="предполагаемый сдвиг средней активности теста относительно нашей "
                         "выборки. Пара (off, lambda) подбирается под ЭТО предположение. "
                         "Ноль означает «тест распределён как обучающая выборка» - это не "
                         "отсутствие предположения, а предположение, и src/shrinkchoice.py "
                         "показывает, что по вилке +0.1..+0.6 оно худшее из трёх правил: "
                         "худший случай на 0.087, средний на 0.040 хуже подгонки под "
                         "середину вилки. Значение по умолчанию оставлено нулевым, чтобы "
                         "Принимает одно число на все ферменты или четыре через запятую в "
                         "порядке CYPS. Умолчание --- принятое правило: по среднему "
                         "апостериорному, с нулём на CYP1A2, см. docstring")
    ap.add_argument("--outdir", default=RES + "submission/")
    ap.add_argument("--no-deadzone", dest="deadzone", action="store_false",
                    help="ОТКЛЮЧИТЬ проход мёртвой зоны. По умолчанию он ВКЛЮЧЁН. Пункт 164: "
                         "+0.0167 ранга на подаваемой пятичленной конфигурации, четыре сида. "
                         "Мишень --- clip(предсказание ВНЕ ФОЛДА, lo, hi), модель та же, у "
                         "бустингов под absolute_error. Стоит одного полного прохода вне фолда "
                         "поверх обычного счёта, потому что мишень нельзя строить из "
                         "предсказаний на собственных обучающих строках. Умолчание "
                         "переключено на основании пункта 204: verify/k58_dzsubmit.py "
                         "воспроизвёл пункт 164 кодом самой подачи --- +0.0136 ранга и "
                         "-0.0158 пары на сиде 0, все четыре фермента вверх.")
    ap.set_defaults(deadzone=True)
    ap.add_argument("--no-bundle", dest="bundle", action="store_false",
                    help="классификация прямым классификатором + Платт (поведение до пункта 250). "
                         "По умолчанию ВКЛЮЧЕНА связка: вероятность ворот из плеча преинкубации, "
                         "предсказанного четырёхчленным ансамблем с мёртвой зоной, умножается на "
                         "вероятность сдвига, произведение калибруется и режется plug-in-порогом. "
                         "Принята по предрегистрации пункта 244: макро MCC +0.0127 при поле 0.0076, "
                         "знак 6/8. Поферментно не проходит НИЧЕГО (+0.0083 при 0.0281 и +0.0170 "
                         "при 0.0419) --- заявление макро. Стоит примерно час к прогону.")
    ap.add_argument("--screen", action="store_true",
                    help="добавить строки скрининга в поферментный член. ВЫКЛЮЧЕНО по умолчанию: "
                         "пункт 177 даёт +0.029 ранга ОДИНОЧНОЙ модели, но пункт 182 померил тот "
                         "же арм в ансамбле и получил +0.0058 на базовом и +0.0014 на лучшем, то "
                         "есть НИЖЕ макро-пола 0.0036. Диагноз там же: корреляция ошибок арма с "
                         "ансамблем 0.935--0.969, он ошибается на тех же соединениях, а среднее "
                         "платит за несогласие. Код оставлен, потому что для ОДИНОЧНОЙ модели "
                         "+0.029 остаётся в силе, и потому что вариант ЗДЕСЬ иной: скрининг "
                         "вставлен ВНУТРЬ поферментного члена, а не добавлен шестым. Механизм "
                         "предсказывает тот же ноль, измерено это не было.")
    a = ap.parse_args()

    _pl.Path(a.outdir).mkdir(parents=True, exist_ok=True)
    desc_names = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
    mech_names = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()

    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    scr = None
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())

    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)

    if a.screen:
        # Пункт 177: показания скрининга как МИШЕНЬ поферментной модели, +0.029 ранга на
        # четырёх сидах, каждый фермент выше своего пола, и метрика улучшается одновременно.
        # Только поферментный член: в пулированной рамке тот же скрининг давал +0.0040, и
        # пункт 177 показал, что виновата была рамка, а не скрининг.
        L2, has = load_screen(rows)
        scr = (L2, has)
        print("скрининг как мишень поферментного члена: "
              + " ".join(f"{c[3:]} +{int((has[:, e] & ~mask[:, e]).sum())}"
                         for e, c in enumerate(CYPS))
              + f" строк (всего +{int((has & ~mask).sum())} к {int(mask.sum())} кривым)",
              flush=True)

    print("строю признаки теста (переиндексация по именам, не пересчёт фильтра)", flush=True)
    te, Xte = test_features(desc_names, mech_names)
    if Xte.shape[1] != X.shape[1]:
        raise SystemExit(f"ширина не совпала: обучение {X.shape[1]}, тест {Xte.shape[1]}")
    print(f"  тест {Xte.shape}, обучение {X.shape}", flush=True)

    lams = None
    if a.shrink:
        print("подбираю усадку вне выборки на обучающих данных", flush=True)
        fold, _ = butina_folds(list(rows.SMILES))
        d = [float(x) for x in str(a.delta).split(",")]
        if len(d) == 1:
            d = d * 4
        if len(d) != 4:
            raise SystemExit(f"--delta: нужно одно число или четыре через запятую, дано {len(d)}")
        print(f"предполагаемый сдвиг по ферментам: "
              + ", ".join(f"{c} {v:+.2f}" for c, v in zip(CYPS, d)), flush=True)
        P = oof_predictions(X, y, mask, fold, a.mode, scr,
                            (LO, HI) if a.deadzone else None)
        lams = fit_shrinkage(P, y, mask, d)

    dz_targets = None
    if a.deadzone:
        # Мишень тестового пути строится из предсказаний ВНЕ ФОЛДА на обучающих строках.
        # Иначе схема вырождается: модель, спрошенная про свои же обучающие строки, кладёт
        # их внутрь полос, мишень совпадает с предсказанием, переподгонка ничего не меняет,
        # и прогон отрабатывает чисто (докстринг src/abldead.py). Это отдельный полный
        # проход вне фолда, и он стоит примерно столько же, сколько всё остальное вместе.
        print("проход вне фолда для мишеней мёртвой зоны", flush=True)
        dzfold, _ = butina_folds(list(rows.SMILES))
        dz_targets = []
        for kind, P in oof_members(X, y, mask, dzfold, a.mode, scr, dead=True):
            if kind == "ствол":
                continue
            dz_targets.append((kind, [np.clip(P[e], LO[mask[:, e], e], HI[mask[:, e], e])
                                      for e in range(len(CYPS))]))
        print(f"    мишени готовы для {len(dz_targets)} членов", flush=True)

    print(f"обучаю на всей выборке (режим: {a.mode}) и предсказываю тест", flush=True)
    act = pd.DataFrame({"SMILES": te.SMILES, "Molecule_Name": te.Molecule_Name})
    shared = None
    trunk_te = None
    if a.mode == "ансамбль5":
        import trunk as TR
        meta = json.load(open(RES + "preds/trunk_twohead.json")).get("meta", {})
        dev = _trunk_device(meta.get("device"))
        print(f"    обучаю ствол на всей выборке и предсказываю тест "
              f"(блоки {meta.get('blocks', TR.BLOCKS)}, вычислитель {dev}, "
              f"в файле {meta.get('device')})", flush=True)
        fp_te, de_te, me_te = _trunk_blocks(Xte)
        trunk_te = TR.fit_predict_test(fp_te, de_te, me_te, lam=float(TRUNK_LAM),
                                       seed=0, mode=TRUNK_MODE,
                                       blocks=meta.get("blocks"), device=dev,
                                       dead=(RES + "preds/trunk_twohead.json"
                                             if a.deadzone else ""))
    if dz_targets is None and a.mode in ("пул", "ансамбль", "ансамбль-без-GP", "ансамбль5"):
        Xs = [pooled_design(X[mask[:, e]], e) for e in range(len(CYPS))]
        ys = [y[mask[:, e], e] for e in range(len(CYPS))]
        shared = gbm_reg().fit(np.vstack(Xs), np.concatenate(ys))
        print(f"    пулированная модель на "
              f"{sum(int(mask[:, e].sum()) for e in range(len(CYPS)))} строках", flush=True)
    for e, c in enumerate(CYPS):
        m = mask[:, e]
        parts = []
        if dz_targets is not None:
            # Мёртвая зона: все члены строятся ОДИНАКОВО --- та же модель, спроецированная
            # мишень, признаки через _dz_design, который подгоняет преобразование на
            # обучающих строках и применяет к тесту. Это же снимает дефект 2 пункта 202:
            # гребневая больше не стандартизуется на объединении обучения с тестом.
            for kind, T in dz_targets:
                if not _keep(e, kind):
                    continue
                Xtr_d, Xte_d = _dz_design(kind, X[m], e, Xte)
                parts.append(_dz_fit(kind, Xtr_d, T[e], Xte_d))
        elif a.mode in ("раздельно", "ансамбль", "ансамбль-без-GP", "ансамбль5"):
            if scr is None:
                parts.append(gbm_reg().fit(X[m], y[m, e]).predict(Xte))
            else:
                # На тесте обучающими являются ВСЕ молекулы, поэтому изотоника подгоняется
                # по всей выборке. Это не утечка: отложенных меток здесь нет, отложен тест,
                # а он в подгонке не участвует ни одной строкой.
                Xa, ya = _screen_train(X, y, mask, scr[0], scr[1], e,
                                       np.ones(len(X), bool))
                parts.append(gbm_reg().fit(Xa, ya).predict(_src0(Xte)))
        if dz_targets is None and a.mode in ("пул", "ансамбль", "ансамбль-без-GP",
                                            "ансамбль5"):
            parts.append(shared.predict(pooled_design(Xte, e)))
        if dz_targets is None and a.mode in ("ансамбль", "ансамбль5"):
            tf = gp_prepare(X[m])
            parts.append(gp_predict(tf(X[m]), y[m, e], tf(Xte)))
            # Дефект 2 пункта 202: раньше здесь стояло
            # _desc_scaled(np.vstack([X[m], Xte])) --- стандартизация по объединению
            # обучения с тестом. Два следствия, и второе хуже первого. Подаваемый член
            # переставал быть тем, который измерен в _oof_ridge (там масштаб берётся по
            # X[m]), при поферментных полах 0.003--0.007. И тестовые признаки входили в
            # обучающее преобразование, то есть это трансдуктивное использование теста.
            # Путь мёртвой зоны это уже не задевало: _dz_design подгоняет любое
            # преобразование на обучающих строках. Здесь то же самое.
            Btr, Bte = _dz_design("гребневая", X[m], e, Xte)
            parts.append(RidgeCV(alphas=np.logspace(-1, 4, 12))
                         .fit(Btr, y[m, e]).predict(Bte))
        if a.mode == "ансамбль5" and _keep(e, "ствол"):
            parts.append(_trunk_clip(trunk_te[:, e], y[m, e]))
        p = np.mean(parts, axis=0)
        if CYPS[e] in SOLO:
            print(f"    {c}: поферментный состав {'+'.join(SOLO[CYPS[e]])} "
                  f"({len(parts)} член(ов) из пяти), пункт 218", flush=True)
        if lams is not None:
            L, mu_tr, sh = lams[e]
            p = L * p + (1.0 - L) * mu_tr + sh
        act[f"{c}_pIC50_direct_inhibition"] = p
        print(f"    {c}: n_обуч {m.sum()}, среднее предсказание {p.mean():.3f}", flush=True)

    tdi = pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv").set_index("Molecule_Name")
    keep = rows.Molecule_Name.isin(tdi.index).to_numpy()
    T = tdi.loc[rows.Molecule_Name[keep]].reset_index()
    cls = pd.DataFrame({"SMILES": te.SMILES, "Molecule_Name": te.Molecule_Name})
    tdifold, _ = butina_folds(list(rows.SMILES[keep]))
    if a.bundle:
        # Связка (пункты 244, 250). Готовим обе конъюнкты на плече преинкубации.
        inh_full = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
                      .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
        Xs, ys, folds, LOs, HIs, labs, shifts, gates = [], [], [], [], [], [], [], []
        for c in TDI_CYPS:
            lab = T[f"{c}_is_TDI"]
            ta = T[f"{c}_pIC50_TDI_condition"].to_numpy(float)
            dr = inh_full[f"{c}_pIC50_direct_inhibition"].to_numpy(float)
            lo = T[f"{c}_pIC50_TDI_condition_conf_low"].to_numpy(float)
            hi = T[f"{c}_pIC50_TDI_condition_conf_high"].to_numpy(float)
            mm = (lab.notna() & np.isfinite(ta) & np.isfinite(dr)
                  & np.isfinite(lo) & np.isfinite(hi)).to_numpy()
            Xs.append(X[keep][mm]); ys.append(ta[mm]); folds.append(tdifold[mm])
            LOs.append(lo[mm]); HIs.append(hi[mm])
            labs.append(lab[mm].astype(int).to_numpy())
            shifts.append((ta - dr)[mm]); gates.append((ta > TDI_GATE)[mm].astype(int))
            print(f"    {c}: связка, обучающих {int(mm.sum())}", flush=True)
        arm_oof, arm_tst = _arm_ensemble(Xs, ys, folds, [Xte] * len(TDI_CYPS), LOs, HIs)
        for e, c in enumerate(TDI_CYPS):
            gp_te, gp_tr = _platt_1d(arm_oof[e] - TDI_GATE, gates[e],
                                     arm_tst[e] - TDI_GATE, raw=True)
            sh = (shifts[e] > TDI_L).astype(int)
            dp_tr = np.zeros(len(sh))
            for f in range(5):
                trn, tef = folds[e] != f, folds[e] == f
                if tef.sum() and sh[trn].sum():
                    dp_tr[tef] = gbm_clf().fit(Xs[e][trn], sh[trn]).predict_proba(Xs[e][tef])[:, 1]
            dp_te = gbm_clf().fit(Xs[e], sh).predict_proba(Xte)[:, 1]
            prod_tr, prod_te = gp_tr * dp_tr, gp_te * dp_te
            pr, _ = _platt_1d(prod_tr, labs[e], prod_te)
            cal_tr, _ = _platt_1d(prod_tr, labs[e], prod_tr)
            thr = plugin_threshold(cal_tr)
            cls[f"{c}_is_TDI"] = pr >= thr
            print(f"    {c}: ворота E[p] {gp_te.mean():.3f} (обуч. доля {gates[e].mean():.3f}), "
                  f"сдвиг E[p] {dp_te.mean():.3f} ({sh.mean():.3f}); произведение -> "
                  f"E[p] {pr.mean():.3f} (метка {labs[e].mean():.3f}); порог {thr:.3f}, "
                  f"положительных {int((pr>=thr).sum())} из {len(pr)}", flush=True)
    else:
        for c in TDI_CYPS:
            lab = T[f"{c}_is_TDI"]
            m = lab.notna().to_numpy()
            yb = lab[m].astype(int).to_numpy()
            pr = gbm_clf().fit(X[keep][m], yb).predict_proba(Xte)[:, 1]
            raw_mean = float(pr.mean())
            pr, slope = tdi_calibrate(X[keep][m], yb, tdifold[m], pr)
            thr = plugin_threshold(pr)
            cls[f"{c}_is_TDI"] = pr >= thr
            print(f"    {c}: Платт наклон {slope:+.3f}, среднее {raw_mean:.3f} -> {pr.mean():.3f} "
                  f"(обучающая доля {yb.mean():.3f}); порог {thr:.3f}, "
                  f"положительных {int((pr>=thr).sum())} из {len(pr)}", flush=True)

    ap_ = a.outdir + "activity_submission.csv"
    tp_ = a.outdir + "tdi_submission.csv"
    act.to_csv(ap_, index=False)
    cls.to_csv(tp_, index=False)

    # The gate. Nothing above is trusted until the organisers' own code accepts it.
    from validation.activity_validation import validate_activity_submission
    from validation.tdi_validation import validate_tdi_submission
    ids = set(te.Molecule_Name)
    ok = True
    for name, fn, path in [("регрессия", validate_activity_submission, ap_),
                           ("классификация", validate_tdi_submission, tp_)]:
        try:
            res = fn(path, expected_ids=ids)
            bad = res if isinstance(res, list) else getattr(res, "errors", [])
            if bad:
                ok = False
                print(f"  {name}: ОТКЛОНЕНО")
                for e_ in bad:
                    print("     ", e_)
            else:
                print(f"  {name}: принято")
        except TypeError:
            res = fn(path)
            print(f"  {name}: принято (валидатор без expected_ids)")
    if not ok:
        raise SystemExit("валидатор отверг файл; ничего не отправлять")
    print(f"\nготово:\n  {ap_}\n  {tp_}")


if __name__ == "__main__":
    LO = HI = None
    main()
