"""The oracle gate: how much can ANY per-compound correction add over a monotone one.

Why this exists as a procedure rather than a result. Item 126 closed a proposal that argued it
escaped item 77's ceiling *by construction* -- the correction depended on the posterior spread,
and the spread is not a function of the point prediction. The premise was true and the conclusion
was false: the resulting action came out 97 to 99.7 per cent monotone in the point anyway. The
error has a general form, and it is worth stating as such:

    item 77's ceiling is not a claim about a class of functions. It is a claim about
    **magnitudes**. Showing that a correction lies outside the affine family says nothing;
    what decides is whether its non-affine part survives the curvature of the loss.

So any future proposal of the form "a per-compound correction on top of the predictions" needs an
a-priori bound on the size of its non-affine part, not a proof that one exists. This file computes
that bound once, for all such proposals at once, and it needs nothing built.

The oracle correction. For each compound take the correction that a perfect per-compound layer
would apply -- move the prediction to the nearest point of its own true band:

    c_i = clip(p_i, lo_i, hi_i) - p_i

This drives the ST-RAE numerator to exactly zero, so it is the unimprovable per-compound
correction. It is unachievable by construction: it uses the true label's band. That is the point.
Whatever any real layer does is a subset of what the oracle does.

What the oracle correction actually is. Up to the clipping, c is the negative of the model's
residual: where the prediction misses the band, the correction is exactly the distance back to it.
That identity is checked below and it is what makes the gate general, because it says the target
of any per-compound post-hoc layer **is the model's own error**. A layer that could compute it
from information available at prediction time would be information that belonged in the model.

Two numbers per enzyme therefore, and the second is the one that binds.

  R^2 of isotonic(c | p)   how much of the oracle correction is a monotone function of the
                           prediction. Near zero is expected and is not a licence: it means the
                           correction is not a function of the prediction *at all*, which rules
                           out monotone post-processing and every other function of p alone.

  монотонный потолок       the best monotone map of p, in the metric. Computed by
                           majorise-minimise: isotonic regression onto clip(m, lo, hi), iterated.
                           That is the same reduction src/abldead.py uses -- the soft-threshold
                           loss is L1 against a target reprojected onto the band -- applied here
                           to the monotone class instead of to a learner. The first attempt at
                           this file used p + isotonic(c | p) instead and got a "ceiling" that
                           scored *worse* than the affine pair, which is impossible for a ceiling
                           and is how the error was caught.

The fits are deliberately in-sample. An out-of-fold version would be weaker and make the gate look
more permissive than it is; in-sample is generous to the monotone class and therefore conservative
in the direction a gate needs.

Reads results/preds/oof_pool_all.json, oof_gp.json, oof_weak.json. Writes nothing. ~1 minute.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
SEEDS = (0, 1, 2, 3)
MEMBERS = [("oof_pool_all", "независимо"), ("oof_pool_all", "пул"),
           ("oof_gp", "GP"), ("oof_weak", "гребневая")]


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    cache = {f: json.load(open(RES + f"preds/{f}.json"))["preds"]
             for f in {m[0] for m in MEMBERS}}

    acc = {}
    for seed in SEEDS:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for c in CYPS:
            col = f"{c}_pIC50_direct_inhibition"
            m = tr[col].notna().to_numpy()
            y = tr.loc[m, col].to_numpy()
            lo = tr.loc[m, col + "_conf_low"].to_numpy()
            hi = tr.loc[m, col + "_conf_high"].to_numpy()
            fi = fold[m]
            P = np.column_stack([np.asarray(cache[fn][f"{seed}|{arm}|{c}"], float)
                                 for fn, arm in MEMBERS])
            p = P.mean(1)

            oracle = np.clip(p, lo, hi) - p
            iso = IsotonicRegression(out_of_bounds="clip").fit(p, oracle)
            mono_add = iso.predict(p)
            r2 = 1.0 - ((oracle - mono_add) ** 2).sum() / max(((oracle - oracle.mean()) ** 2).sum(), 1e-12)

            # Тождество: оракульная поправка --- это обрезанный остаток. Когда предсказание
            # ниже полосы, обе величины положительны, поэтому корреляция должна быть БЛИЗКА
            # К +1; знак здесь --- часть проверки, а не оформление.
            resid = y - p
            ident = float(np.corrcoef(oracle, resid)[0, 1])

            # Настоящий монотонный потолок: мажорирование. Потеря с мёртвой зоной есть L1
            # против мишени, спроецированной на полосу, поэтому итерируем изотонику по
            # clip(m, lo, hi) --- та же редукция, что в src/abldead.py.
            m_cur = p.copy()
            for _ in range(30):
                tgt = np.clip(m_cur, lo, hi)
                m_new = IsotonicRegression(out_of_bounds="clip").fit(p, tgt).predict(p)
                if np.max(np.abs(m_new - m_cur)) < 1e-4:
                    m_cur = m_new
                    break
                m_cur = m_new

            u = np.ones(len(y)) / len(y)
            arms = {"аффинная пара": fit_apply(p, lo, hi, fi, u),
                    "монотонный потолок": m_cur,
                    "оракул": p + oracle}
            row = {nm: float(strae(y, q, y_true_upper=hi, y_true_lower=lo))
                   for nm, q in arms.items()}
            row["R2"] = float(r2)
            row["тождество"] = ident
            row["сырое"] = float(strae(y, p, y_true_upper=hi, y_true_lower=lo))
            acc.setdefault(c, []).append(row)
        print(f"  сид {seed} готов", flush=True)

    print()
    print(f"{'фермент':8s} {'сырое':>8s} {'пара':>8s} {'монот. потолок':>15s} {'оракул':>8s}"
          f" {'R2 оракула':>11s} {'бюджет слоя':>12s} {'тождество':>10s}")
    tot = {k: [] for k in ("сырое", "аффинная пара", "монотонный потолок", "оракул",
                           "R2", "тождество")}
    for c in CYPS:
        d = pd.DataFrame(acc[c]).mean()
        for k in tot:
            tot[k].append(d[k])
        print(f"{c:8s} {d['сырое']:8.4f} {d['аффинная пара']:8.4f} "
              f"{d['монотонный потолок']:15.4f} {d['оракул']:8.4f} {d['R2']:11.3f}"
              f" {d['аффинная пара'] - d['монотонный потолок']:12.4f}"
              f" {d['тождество']:10.3f}")
    print(f"{'МАКРО':8s} {np.mean(tot['сырое']):8.4f} {np.mean(tot['аффинная пара']):8.4f} "
          f"{np.mean(tot['монотонный потолок']):15.4f} {np.mean(tot['оракул']):8.4f} "
          f"{np.mean(tot['R2']):11.3f} "
          f"{np.mean(tot['аффинная пара']) - np.mean(tot['монотонный потолок']):12.4f}")

    print("""
Как читать. Столбец «бюджет слоя» --- сколько ещё можно взять, не выходя из монотонного класса,
который аффинная пара приближает двумя параметрами. Он и говорит, стоит ли вообще улучшать
постобработку.

Столбец «тождество» --- корреляция оракульной поправки с остатком. Близко к единице
означает, что улучшать постобработку нечем в принципе: её идеальная мишень --- это ошибка самой
модели, а предсказать ошибку модели из того, что модели уже показали, нельзя по определению.
Это и есть общая форма потолка 77, и она сильнее монотонности: дело не в том, что поправка
обязана быть монотонной, а в том, что её мишень нельзя вычислить.

Столбец R2 читается вместе с этим, а не отдельно. Значение около нуля НЕ означает большого
запаса для посоединённого слоя --- оно означает, что поправка не является функцией предсказания
вовсе, а значит недостижима не только для монотонных преобразований, но и для любых функций
одного лишь p.

Это процедура, а не результат: предложения вида «посоединённая поправка» проверяются здесь
ДО постройки, и проверяются на величину, а не на принадлежность классу.""")


if __name__ == "__main__":
    main()
