"""What leaderboard score the out-of-fold estimate predicts, and what would falsify it.

Written before 24 September rather than after it, because a threshold chosen once the score is
known is not a threshold. Item 145 argued the intermediate leaderboard is the only sample from
the test distribution that will exist before the close; this turns that into a number.

What has to be modelled, and why it is not just "our OOF score". ST-RAE is a ratio and the
organisers' denominator is a constant predictor at `mean(y_true)` put through the same soft
threshold:

    ST-RAE = sum_i max(0, lo_i - p_i, p_i - hi_i) / sum_i max(0, lo_i - ybar, ybar - hi_i)

Both halves are computed on the **test** labels and the **test** bands, neither of which we have.
Item 130 measured which half is the problem: prediction quality contaminates the shift estimate by
only 2 to 12 per cent, while the denominator eats everything, and on CYP2D6 it eats it entirely.
So a leaderboard prediction is not one number with a confidence interval on the model; it is a
ratio of two quantities that both move under the shift.

Three sources of movement, each given its own interval so they can be argued about separately.

  **выборка**   the test is 750 molecules. Even with no shift at all, a 750-molecule draw moves
                the score. Bootstrapped over Butina clusters, because clusters are the unit the
                split assigns and molecules within one are not independent.

  **сид**       f3 measured that the split seed alone moves macro ST-RAE by 0.016. That is
                variation in our own estimate, not in the test, and it belongs in the interval.

  **сдвиг**     item 123 measured the chi-square divergence between the test and out-of-fold
                similarity distributions at 2.838. For any statistic s with variance sigma^2 under
                our distribution P, and any Q with chi2(Q||P) <= rho,

                    |E_Q[s] - E_P[s]| <= sigma * sqrt(rho)

                by Cauchy-Schwarz -- the same bound `src/ablrobust.py` uses to build its weights.
                With rho = 2.838 the multiplier is 1.685.

**And applied to this metric that bound is vacuous, which is the finding.** ST-RAE is a ratio, and
the denominator's per-compound values are mostly zero -- a constant predictor at the mean lands
inside the band for a large share of compounds -- so its standard deviation exceeds its mean and
the bound permits a reweighting that drives the denominator to zero. The worst-case score is then
unbounded above, and no leaderboard value whatsoever could fall outside it. A bound that excludes
nothing excludes nothing, however rigorously it was derived.

So the script reports the two halves separately: the numerator band, which is what "our model
under the shift" means and is finite, and the denominator band, which is where the vacuity comes
from. The usable pre-registration is stated under a **named** assumption -- that the test's label
spread and band widths resemble the training set's -- rather than under a hidden one. Naming it is
the point: it is exactly the assumption items 123 and 130 say cannot be checked from here.

Reads results/preds/oof.json and whatever ensemble file is named. Writes nothing. ~3 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json

import numpy as np
import pandas as pd

from cypsplit import butina_folds, cluster_ids
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
CHI2 = 2.838          # пункт 123: расхождение тест / вне-фолда по сходству
SEED_SPREAD = 0.016   # f3: сид расщепления двигает macro ST-RAE на столько
NTEST = 750           # молекул в ослеплённом тесте
NBOOT = 2000


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--preds", default="oof.json")
    ap.add_argument("--arm", default="FP+DESC+MECH")
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    blob = json.load(open(RES + "preds/" + a.preds))
    P = blob.get("preds", blob)
    fold, ncl = butina_folds(list(rows.SMILES))
    # Кластеры Бутины --- единица, которой оперирует расщепление; бутстрап идёт по ним,
    # иначе соседи внутри кластера считались бы независимыми и интервал вышел бы узким.
    # cluster_ids отдаёт КОРТЕЖ (cid, n_clusters); молчаливого отката здесь нет намеренно ---
    # первая версия имела try/except и тихо бутстрапила по молекулам, давая узкий интервал.
    clus, ncl2 = cluster_ids(list(rows.SMILES))
    clus = np.asarray(clus)
    assert len(np.unique(clus)) == ncl2 == ncl, "кластеризация разошлась с расщеплением"
    rng = np.random.default_rng(0)

    print(f"кластеров Бутины {len(np.unique(clus))}, молекул {len(rows)}, "
          f"тест {NTEST} молекул\n")
    print(f"{'фермент':9s} {'ST-RAE OOF':>11s} {'выборка (750)':>15s} "
          f"{'числитель chi2':>17s} {'знаменатель chi2':>19s}")

    per, boots, bands, dead = {}, {}, {}, {}
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy()
        lo = tr.loc[m, col + "_conf_low"].to_numpy()
        hi = tr.loc[m, col + "_conf_high"].to_numpy()
        fi, cl = fold[m], clus[m]
        key = f"{a.arm}|{c}" if f"{a.arm}|{c}" in P else f"0|{a.arm}|{c}"
        p = fit_apply(np.asarray(P[key], float), lo, hi, fi, np.ones(len(y)) / len(y))

        num = np.maximum(0.0, np.maximum(lo - p, p - hi))
        den0 = np.maximum(0.0, np.maximum(lo - y.mean(), y.mean() - hi))
        per[c] = num.sum() / den0.sum()

        # --- выборка: бутстрап по кластерам до размера теста ---
        uc = np.unique(cl)
        frac = NTEST / len(rows)
        k = max(1, int(round(frac * len(uc))))
        vals = []
        for _ in range(NBOOT):
            pick = rng.choice(uc, size=k, replace=True)
            idx = np.concatenate([np.where(cl == u)[0] for u in pick])
            if len(idx) < 20:
                continue
            yb = y[idx]
            db = np.maximum(0.0, np.maximum(lo[idx] - yb.mean(), yb.mean() - hi[idx]))
            if db.sum() <= 0:
                continue
            vals.append(num[idx].sum() / db.sum())
        boots[c] = np.percentile(vals, [2.5, 97.5]) if vals else (np.nan, np.nan)

        # --- сдвиг: граница Коши-Буняковского, порознь на числитель и знаменатель ---
        mult = np.sqrt(CHI2)
        n_lo = max(num.mean() - mult * num.std(), 0.0)
        n_hi = num.mean() + mult * num.std()
        d_lo = den0.mean() - mult * den0.std()
        d_hi = den0.mean() + mult * den0.std()
        # Числитель при НАБЛЮДЁННОМ знаменателе: это «наша модель под сдвигом».
        bands[c] = (n_lo / den0.mean(), n_hi / den0.mean())
        dead[c] = d_lo <= 0            # знаменатель может уйти в ноль -> отношение без границы

        print(f"{c:9s} {per[c]:11.4f} {boots[c][0]:7.3f}-{boots[c][1]:<7.3f} "
              f"{bands[c][0]:7.3f}-{bands[c][1]:<9.3f} "
              f"{d_lo:8.4f}..{d_hi:<8.4f}{'  ВЫРОЖДЕНА' if dead[c] else ''}")

    macro = float(np.mean([per[c] for c in CYPS]))
    mb = (float(np.mean([boots[c][0] for c in CYPS])),
          float(np.mean([boots[c][1] for c in CYPS])))
    ms = (float(np.mean([bands[c][0] for c in CYPS])),
          float(np.mean([bands[c][1] for c in CYPS])))
    print(f"\n{'МАКРО':9s} {macro:11.4f} {mb[0]:7.3f}-{mb[1]:<7.3f} "
          f"{ms[0]:7.3f}-{ms[1]:<9.3f}")

    nvac = sum(dead.values())
    lo_all = min(mb[0], ms[0]) - SEED_SPREAD
    hi_all = max(mb[1], ms[1]) + SEED_SPREAD
    print(f"""
ПРЕДРЕГИСТРАЦИЯ, записана до 24 сентября.

  точечная оценка вне фолда                    {macro:.4f}
  выборочный разброс на 750 молекулах          {mb[0]:.3f} .. {mb[1]:.3f}
  сид расщепления (f3)                         +/- {SEED_SPREAD}
  числитель под сдвигом chi2={CHI2}             {ms[0]:.3f} .. {ms[1]:.3f}
  ВСЁ ВМЕСТЕ, при названном допущении          {lo_all:.3f} .. {hi_all:.3f}

Допущение названо, а не спрятано: разброс тестовых меток и ширины их полос похожи на
обучающие. Без него интервала нет вовсе, и вот почему.

  знаменатель ушёл в ноль в пределах хи-квадрат-шара: {nvac} фермента(ов) из 4.

Знаменатель ST-RAE --- константный предсказатель на среднем, и у большой доли соединений он
попадает ВНУТРЬ полосы, давая ровно ноль. Поэтому его стандартное отклонение больше среднего,
и граница Коши-Буняковского при радиусе {CHI2} допускает перевзвешивание, обнуляющее знаменатель.
Отношение тогда не ограничено сверху ничем.

**Строгая граница на эту метрику пуста.** Никакой счёт лидерборда не может выйти за неё, а
значит и опровергнуть ею ничего нельзя. Это не изъян вывода --- это свойство отношения, у
которого знаменатель почти всюду нулевой, и оно приходит к тому же, к чему пункты 123 и 129
пришли с двух других сторон.

Что из этого следует ДО 24 сентября, а не после:

  счёт вне интервала при названном допущении опровергает само допущение --- то есть говорит,
  что тестовые полосы устроены не как обучающие. Это единственное, что лидерборд здесь может
  опровергнуть, и это стоит знать;

  счёт внутри не подтверждает ничего;

  поэтому подачу выбирать надо не под «подтвердить оценку», а под РАЗНОСТЬ двух подач
  (пункт 145): разность делит один и тот же неизвестный знаменатель и потому от него не
  зависит. Это единственная величина, которую лидерборд измеряет чисто.""")


if __name__ == "__main__":
    main()
