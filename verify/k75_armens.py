"""Пятичленную машинку подачи --- на плечо преинкубации. Единственный незакрытый маршрут.

Пункт 238 оставил ровно один вывод. Ворота метки TDI --- это порог `pi_TDI > 4.301` на величине,
которую мы предсказываем напрямую; их оракул даёт MCC +0.5667 на CYP3A4, больше всего подаваемого
классификатора, а восстанавливаем мы 53.7 процента. Но предсказаны они были ГОЛЫМ HistGB, тогда
как регрессионный трек подачи --- четыре члена плюс мёртвая зона, +0.058 ранга над той же базой,
--- на плечо преинкубации не направляли ни разу.

Что переносится, проверено по данным, а не предположено:

  * полосы у плеча TDI ЕСТЬ (`{CYP}_pIC50_TDI_condition_conf_low/_conf_high`), значит мёртвая
    зона определена --- это было главное сомнение;
  * пулирование возможно, но по ДВУМ эндпоинтам вместо четырёх (2334 и 1493 строки), поэтому
    индикатор здесь двухпозиционный, а не четырёх --- числа с таблицами регрессионного трека
    построчно не сравнимы;
  * GP и гребневая переносятся как есть;
  * ствол НЕ переносится: его переподгонка живёт в src/trunk.py под torch и обучена на четырёх
    прямых мишенях. Пункт 120 оценил его вклад в +0.0054, так что это четырёхчленная версия
    пятичленной машинки, и так она и называется.

**Метрика здесь не RMSE и даже не ранг целиком.** Ворота --- это решение на пороге, поэтому
меряется AUC предсказанного плеча ПРОТИВ ИСТИННЫХ ВОРОТ (у голого HistGB 0.8709 на CYP3A4) и
итоговый MCC метки. Пункт 238 показал, что срез на усаженном предсказании нельзя ставить в
книжную точку, поэтому срез подбирается вне фолда всюду.

Members composed from `src/submit.py` itself --- `gbm_reg`, `_dz_design`, `_dz_fit`, `gp_prepare`,
`gp_predict` --- по образцу k58, чтобы члены не разъехались с подаваемыми.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES

import argparse
import json
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import matthews_corrcoef, roc_auc_score

from cypsplit import butina_folds
from submit import gbm_reg, _dz_design, _dz_fit, DZ_KW
from sklearn.ensemble import HistGradientBoostingRegressor

TDI = ["CYP3A4", "CYP2D6"]
L = np.log10(2.0)
GATE = 4.0 + L


def pooled2(X, e):
    """Индикатор на ДВА эндпоинта. Не pooled_design из подачи: там их четыре."""
    ind = np.zeros((len(X), 2), np.float32)
    ind[:, e] = 1.0
    return np.hstack([X, ind])


def _fit_one(kind, Xs, ys, folds, e, f, KW):
    """Одна подгонка одного члена на одном фолде. Общая для обычного прохода и для
    мёртвой зоны --- отличаются они только мишенью и пинами."""
    X, y, fi = Xs[e], ys[e], folds[e]
    trn, te = fi != f, fi == f
    if kind == "пул":
        # Строки ОБОИХ эндпоинтов, из каждого исключён тот же фолд по молекулам, так что
        # ни одна молекула тестового фолда не попадает в обучение ни с какой мишенью.
        Xtr = np.vstack([pooled2(Xs[k][folds[k] != f], k) for k in range(len(TDI))])
        ytr = np.concatenate([ys[k][folds[k] != f] for k in range(len(TDI))])
        return HistGradientBoostingRegressor(**KW).fit(Xtr, ytr).predict(pooled2(X[te], e))
    if kind in ("GP", "гребневая"):
        A, B = _dz_design(kind, X[trn], e, X[te])
        return _dz_fit(kind, A, y[trn], B)
    return HistGradientBoostingRegressor(**KW).fit(X[trn], y[trn]).predict(X[te])


PLAIN = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
             l2_regularization=1.0, random_state=0)


def member_oof(kind, Xs, ys, folds, KW=None):
    """Вне фолда для одного члена по обоим эндпоинтам.

    Пулированный член обучается ОДИН раз на фолд: его обучающая выборка --- строки обоих
    эндпоинтов без этого фолда --- от эндпоинта не зависит, зависит только предсказание.
    Подгонять её дважды было бы вдвое дороже при том же ответе."""
    KW = KW or PLAIN
    out = [np.zeros(len(ys[e])) for e in range(len(TDI))]
    if kind == "пул":
        for f in range(5):
            Xtr = np.vstack([pooled2(Xs[k][folds[k] != f], k) for k in range(len(TDI))])
            ytr = np.concatenate([ys[k][folds[k] != f] for k in range(len(TDI))])
            mdl = HistGradientBoostingRegressor(**KW).fit(Xtr, ytr)
            for e in range(len(TDI)):
                te = folds[e] == f
                if te.sum():
                    out[e][te] = mdl.predict(pooled2(Xs[e][te], e))
        return out
    for e in range(len(TDI)):
        for f in range(5):
            if (folds[e] == f).sum() == 0:
                continue
            out[e][folds[e] == f] = _fit_one(kind, Xs, ys, folds, e, f, KW)
    return out


def dz_refit(kind, Xs, targets, folds):
    """Переподгонка члена под спроецированную мишень. Пины DZ_KW --- absolute_error."""
    return member_oof(kind, Xs, targets, folds, KW=DZ_KW)


def fit_cut(pred, truth, fold):
    g = np.quantile(pred, np.linspace(0.02, 0.98, 97))
    out = np.zeros(len(pred), bool)
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0:
            continue
        b = max(g, key=lambda t: matthews_corrcoef(truth[trn], pred[trn] >= t))
        out[te] = pred[te] >= b
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/armens.json")
    a = ap.parse_args()

    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tdi = pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv").set_index("Molecule_Name").reindex(rows.Molecule_Name)
    inh = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv").set_index("Molecule_Name").reindex(rows.Molecule_Name)

    recs = []
    for seed in [int(s) for s in a.seeds.split(",")]:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        print(f"\n=== сид {seed} ===", flush=True)
        Xs, ys, folds, LOs, HIs, labs, shifts, gates = [], [], [], [], [], [], [], []
        for c in TDI:
            lab = tdi[f"{c}_is_TDI"]
            ta = tdi[f"{c}_pIC50_TDI_condition"].to_numpy(float)
            dr = inh[f"{c}_pIC50_direct_inhibition"].to_numpy(float)
            lo = tdi[f"{c}_pIC50_TDI_condition_conf_low"].to_numpy(float)
            hi = tdi[f"{c}_pIC50_TDI_condition_conf_high"].to_numpy(float)
            m = (lab.notna() & np.isfinite(ta) & np.isfinite(dr) &
                 np.isfinite(lo) & np.isfinite(hi)).to_numpy()
            Xs.append(X[m]); ys.append(ta[m]); folds.append(fold[m])
            LOs.append(lo[m]); HIs.append(hi[m])
            labs.append(lab[m].astype(bool).to_numpy())
            shifts.append((ta - dr)[m]); gates.append((ta > GATE)[m])
            print(f"  {c}: n {int(m.sum())}, ворота открыты {gates[-1].mean():.3f}", flush=True)

        # предсказание сдвига, чтобы собрать полное решение метки
        sh_hat = []
        for e in range(len(TDI)):
            p = np.zeros(len(shifts[e]))
            for f in range(5):
                trn, te = folds[e] != f, folds[e] == f
                p[te] = gbm_reg().fit(Xs[e][trn], shifts[e][trn]).predict(Xs[e][te])
            sh_hat.append(p)

        KINDS = ["поферментно", "пул", "GP", "гребневая"]
        members, dzmembers = {}, {}
        for k in KINDS:
            t0 = time.time()
            members[k] = member_oof(k, Xs, ys, folds)
            tg = [np.clip(members[k][e], LOs[e], HIs[e]) for e in range(len(TDI))]
            dzmembers[k] = dz_refit(k, Xs, tg, folds)
            print(f"    {k}: {time.time()-t0:.0f} с", flush=True)

        scores = {k: members[k] for k in KINDS}
        scores["АНСАМБЛЬ4"] = [np.mean([members[k][e] for k in KINDS], 0) for e in range(len(TDI))]
        scores["АНСАМБЛЬ4 + мёртвая зона"] = [
            np.mean([dzmembers[k][e] for k in KINDS], 0) for e in range(len(TDI))]

        for name, P in scores.items():
            for e, c in enumerate(TDI):
                p, g, y, fi = P[e], gates[e], labs[e], folds[e]
                rmse = float(np.sqrt(np.mean((p - ys[e]) ** 2)))
                rho = spearmanr(ys[e], p).statistic
                auc = roc_auc_score(g, p)
                gh = fit_cut(p, g, fi)
                mcc_gate = matthews_corrcoef(y, gh)
                mcc_full = matthews_corrcoef(y, gh & (sh_hat[e] > L))
                print(f"      {name:26s} {c}  RMSE {rmse:.4f} rho {rho:.4f} "
                      f"AUC(ворота) {auc:.4f}  MCC ворот {mcc_gate:+.4f}  "
                      f"MCC метки {mcc_full:+.4f}", flush=True)
                recs.append({"seed": seed, "cyp": c, "score": name, "rmse": rmse,
                             "rho": float(rho), "auc": float(auc),
                             "mcc_gate": float(mcc_gate), "mcc_full": float(mcc_full)})
        json.dump(recs, open(a.out, "w"))

    df = pd.DataFrame(recs)
    for met in ("auc", "mcc_gate", "mcc_full", "rmse"):
        print(f"\n{met}, среднее по сидам")
        print(df.pivot_table(index="score", columns="cyp", values=met).round(4).to_string())
    print(f"\nсохранено: {a.out}", flush=True)


if __name__ == "__main__":
    main()
