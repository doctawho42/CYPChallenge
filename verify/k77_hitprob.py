"""P(попасть в полосу) --- объект неопределённости, специфичный ИМЕННО этой метрике.

Обычная количественная оценка неопределённости строит интервал вокруг предсказания. Здесь это
не тот объект. ST-RAE даёт РОВНО НОЛЬ за предсказание, попавшее внутрь опубликованной полосы
соединения, и линейную цену за промах. Значит родная величина --- не интервал вокруг y-крышка, а

    P(попадание_i) = P( lo_i <= y-крышка_i <= hi_i )

то есть поштучная вероятность НЕ ЗАПЛАТИТЬ. Она наблюдаема на обучении вне фолда, предсказуема
из того, что известно на тесте, и складывается в предсказание собственного счёта.

**Почему это вообще решаемо, и это не наша заслуга, а свойство их данных.** Пункт 114: ширина
полосы --- детерминированная функция метки с точностью до 3 процентов дисперсии (изотоника от y
на ширину, R^2 0.928--0.970), и корреляция метки с шириной -0.885..-0.928. У СЛАБЫХ соединений
полосы ШИРОКИЕ. Поэтому доля попаданий определяется прежде всего потенцией, а её мы предсказываем:
на CYP3A4 по квинтилям предсказанной потенции доля попаданий идёт 0.788 -> 0.161 при ширине
2.00 -> 0.21. Пятикратный разброс, весь видимый из выхода самой модели.

**Пара к мёртвой зоне.** Мёртвая зона учит предсказание ВНУТРЬ полосы; этот файл говорит, попали
ли. Две половины одного объекта, а не два приёма.

**Признаки берутся только те, что существуют на тесте.** Предсказание ансамбля; максимальное
танимото до обучающего набора (расстояние до многообразия); индикатор фермента. Разброс ансамбля
(пункт 222, лучший нормировщик из имеющихся) требует предсказаний отдельных членов, которых на
диске нет; как ВЕРХНЮЮ ОЦЕНКУ того, что он мог бы добавить, считается разногласие шести наборов
признаков из oof.json --- это не деплоируемый признак и он помечен как справочный.

**Один баг стоит отдельной строки, потому что его поймал контроль, а не глаз.** Первая версия
признака z проходила через `IsotonicRegression().fit(y, w)`, а он по умолчанию требует
НЕУБЫВАНИЯ. Связь потенции с шириной УБЫВАЮЩАЯ (корреляция -0.46 и -0.71), поэтому подгонка
схлопывалась в константу --- размах ровно 0.000 --- и z оказывался обратным преобразованием одного
лишь масштаба остатка. Поймано на том, что корреляция предсказанной ширины с попаданием вышла
+-0.01 на всех четырёх ферментах, чего не может быть у монотонной функции информативного признака.
Нужен `increasing=False`.

**Что меряется.** Различение (AUC), калибровка (Брайер, ECE, кривая надёжности) и агрегат:
сходится ли сумма предсказанных вероятностей с настоящим числом попаданий. Плюс вторая голова ---
E[L] --- дающая предсказанный ST-RAE, который сравнивается с настоящим вне фолда. Последнее и есть
проверка, что вся конструкция работает: предсказание собственного счёта, а не только рассказ о нём.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES

import argparse, json
import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import roc_auc_score, brier_score_loss

from cypsplit import butina_folds, fingerprints
from rdkit import DataStructs

CY = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=200, learning_rate=0.06, max_leaf_nodes=15,
          l2_regularization=1.0, random_state=0)


def maxsim(smiles, fold):
    """Максимальное танимото до обучающей части СВОЕГО фолда --- то, что видно на тесте."""
    fps = fingerprints(smiles)
    out = np.zeros(len(fps))
    for f in range(5):
        te = np.where(fold == f)[0]
        trn = np.where(fold != f)[0]
        ref = [fps[i] for i in trn]
        for i in te:
            out[i] = max(DataStructs.BulkTanimotoSimilarity(fps[i], ref))
    return out


def ece(p, y, bins=10):
    e, n = 0.0, len(y)
    for k in range(bins):
        m = (p >= k / bins) & (p < (k + 1) / bins if k < bins - 1 else p <= 1.0)
        if m.sum():
            e += m.sum() / n * abs(y[m].mean() - p[m].mean())
    return e


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=RES + "preds/hitprob.json")
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    P = json.load(open(RES + "preds/oof_dzens.json"))["preds"]
    V = json.load(open(RES + "preds/oof.json"))
    SETS = ["FP", "DESC", "MECH", "FP+DESC", "DESC+MECH", "FP+DESC+MECH"]
    fold, _ = butina_folds(list(rows.SMILES))
    print("считаю максимальное танимото вне фолда", flush=True)
    sim_all = maxsim(list(rows.SMILES), fold)

    recs, curves = [], {}
    for c in CY:
        y = tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float)
        m = np.isfinite(y)
        lo = tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float)[m]
        hi = tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float)[m]
        yv, fi, sim = y[m], fold[m], sim_all[m]
        yhat = np.asarray(P[f"0|мёртвая зона везде|{c}"], float)
        hit = ((yhat >= lo) & (yhat <= hi)).astype(int)
        L = np.maximum(0.0, np.maximum(lo - yhat, yhat - hi))
        disagree = np.std(np.stack([np.asarray(V[f"{s}|{c}"], float) for s in SETS]), 0)

        FEAT = {"ŷ": yhat[:, None],
                "ŷ + сходство": np.c_[yhat, sim],
                "ŷ + сходство + разногласие (справочно)": np.c_[yhat, sim, disagree]}
        print(f"\n=== {c}: n {len(yv)}, попаданий {hit.mean():.3f}, "
              f"ST-RAE вне фолда {L.sum()/np.maximum(0,np.maximum(lo-yv.mean(),yv.mean()-hi)).sum():.4f}",
              flush=True)
        for name, X in FEAT.items():
            ph = np.zeros(len(hit)); el = np.zeros(len(hit))
            for f in range(5):
                trn, te = fi != f, fi == f
                ph[te] = HistGradientBoostingClassifier(**KW).fit(X[trn], hit[trn]).predict_proba(X[te])[:, 1]
                el[te] = np.maximum(0.0, HistGradientBoostingRegressor(**KW).fit(X[trn], L[trn]).predict(X[te]))
            den = np.maximum(0, np.maximum(lo - yv.mean(), yv.mean() - hi)).sum()
            print(f"  {name:40s} AUC {roc_auc_score(hit,ph):.4f}  Брайер {brier_score_loss(hit,ph):.4f}  "
                  f"ECE {ece(ph,hit):.4f}  попаданий предск. {ph.sum():.0f} против {hit.sum()}  "
                  f"ST-RAE предск. {el.sum()/den:.4f} против {L.sum()/den:.4f}", flush=True)
            recs.append({"cyp": c, "features": name, "auc": float(roc_auc_score(hit, ph)),
                         "brier": float(brier_score_loss(hit, ph)), "ece": float(ece(ph, hit)),
                         "hits_pred": float(ph.sum()), "hits_true": int(hit.sum()),
                         "strae_pred": float(el.sum()/den), "strae_true": float(L.sum()/den),
                         "n": int(len(hit)), "base_rate": float(hit.mean())})
            if name == "ŷ + сходство":
                curves[c] = {"p": ph.tolist(), "hit": hit.tolist()}
        # контроль: константа на базовой доле
        const = np.full(len(hit), hit.mean())
        print(f"  {'КОНТРОЛЬ: константа = базовая доля':40s} AUC 0.5000  "
              f"Брайер {brier_score_loss(hit,const):.4f}  ECE {ece(const,hit):.4f}", flush=True)
        recs.append({"cyp": c, "features": "КОНТРОЛЬ константа", "auc": 0.5,
                     "brier": float(brier_score_loss(hit, const)), "ece": float(ece(const, hit)),
                     "hits_pred": float(const.sum()), "hits_true": int(hit.sum()),
                     "strae_pred": np.nan, "strae_true": np.nan,
                     "n": int(len(hit)), "base_rate": float(hit.mean())})

    df = pd.DataFrame(recs)
    print("\n\nСВОДКА")
    print(df.pivot_table(index="features", columns="cyp", values="auc").round(4).to_string())
    print("\nкалибровка (Брайер, меньше лучше)")
    print(df.pivot_table(index="features", columns="cyp", values="brier").round(4).to_string())
    d = df[df.features == "ŷ + сходство"]
    print(f"\nМАКРО: AUC {d.auc.mean():.4f}, ECE {d.ece.mean():.4f}, "
          f"предсказанный ST-RAE {d.strae_pred.mean():.4f} против настоящего {d.strae_true.mean():.4f}")
    json.dump({"table": recs, "curves": curves}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}", flush=True)


if __name__ == "__main__":
    main()
