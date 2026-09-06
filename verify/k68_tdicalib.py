"""Порог классификационного трека: plug-in против откалиброванного и против подогнанного.

Что здесь проверяется и почему это не повтор f10_calib.py.

Подача ставит порог правилом plug-in: максимизирует MCC, ОЖИДАЕМЫЙ при собственных
вероятностях модели. Правило не требует меток --- и ровно поэтому верно лишь настолько,
насколько вероятности откалиброваны. Пункт 165 назвал это неоговорённым допущением,
пункт 229 показал, что допущение нарушено: средняя предсказанная доля 0.079 против
истинной 0.217 на CYP2D6.

Три вещи здесь сделаны иначе, чем в f10.

  1. Калибровка подгоняется на фолдах Бутины, а не на случайных группах
     (`rng.integers(0,5,n)` в f10). Случайные группы разводят близкие аналоги по разные
     стороны, так что измеренная там выгода калибровки завышена.
  2. Вложение полное: вероятности внешнего фолда получены моделью, не видевшей его, а
     калибратор подогнан на вневыборочных вероятностях ВНУТРЕННЕГО разбиения обучающей
     части. В f10 калибратор подгонялся на вневыборочных вероятностях всей выборки, в
     производстве которых участвовали метки оцениваемого фолда.
  3. Добавлено плечо, которого не было ни там, ни в подаче: порог, ПОДОГНАННЫЙ по меткам
     вне фолда.

Третье --- главное. Платт --- строго возрастающее отображение, изотоническая ---
неубывающее. Порядок они не меняют: AUC тот же, и MCC при ЭМПИРИЧЕСКИ оптимальном пороге
тот же с точностью до слипания узлов. Значит, калибровка может помочь ровно одним каналом:
подвинуть plug-in-порог туда, где и так стоял бы подогнанный. Подогнанный порог её
мажорирует по построению --- с точностью до собственного переобучения, которое здесь и
меряется. А метки у нас есть: 1495 и 2346 обучающих молекул.

Плечи:
    сырые + plug-in        --- то, что подаётся сегодня
    Платт + plug-in        --- предложение из пункта 229
    изотон + plug-in
    сырые + порог вне фолда --- подогнан по меткам внутреннего разбиения
    оракул порога          --- лучший эмпирический порог на самом внешнем фолде (потолок)

Цена: 5 внешних * (1 + 4 внутренних) обучений на фермент на сид. Около часа на сид.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial

import argparse
import json
import time

import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import matthews_corrcoef, roc_auc_score, brier_score_loss

from cypsplit import butina_folds
from submit import gbm_clf, plugin_threshold, TDI_CYPS

GRID = np.linspace(0.05, 0.95, 91)      # та же сетка, что у plugin_threshold


def logit(p):
    q = np.clip(p, 1e-6, 1 - 1e-6)
    return np.log(q / (1 - q)).reshape(-1, 1)


def platt(p_fit, y_fit, p_new):
    """Возвращает откалиброванные вероятности И наклон.

    Платт монотонно ВОЗРАСТАЕТ только при положительном наклоне. На классификаторе без
    порядка наклон может выйти отрицательным, и тогда "калибровка" молча переворачивает
    ранжирование. У CYP2D6 AUC 0.59 --- достаточно слабо, чтобы это случилось в отдельном
    фолде, поэтому наклон здесь не предполагается, а печатается."""
    lr = LogisticRegression(C=1e6).fit(logit(p_fit), y_fit)
    return lr.predict_proba(logit(p_new))[:, 1], float(lr.coef_[0, 0])


def isotonic(p_fit, y_fit, p_new):
    ir = IsotonicRegression(out_of_bounds="clip", y_min=0, y_max=1)
    return ir.fit(p_fit, y_fit).predict(p_new)


def fitted_threshold(p, y):
    """Порог по ЭМПИРИЧЕСКОМУ MCC на размеченных строках."""
    return float(GRID[int(np.argmax([matthews_corrcoef(y, (p >= t).astype(int))
                                     for t in GRID]))])


def run_enzyme(Xi, yi, fi, tag):
    """Полное вложение. Возвращает решения всех плеч по всем строкам."""
    n = len(yi)
    out = {k: np.zeros(n, bool) for k in
           ["сырые+plug", "Платт+plug", "изотон+plug", "сырые+подогн", "оракул"]}
    p_raw = np.zeros(n)
    p_pla = np.zeros(n)
    p_iso = np.zeros(n)
    thrs = {k: [] for k in out}
    slopes = []
    for f in range(5):
        trn, te = fi != f, fi == f
        if te.sum() == 0 or yi[trn].sum() == 0:
            continue
        t0 = time.time()
        clf = gbm_clf().fit(Xi[trn], yi[trn])
        p_raw[te] = clf.predict_proba(Xi[te])[:, 1]

        # Вневыборочные вероятности ВНУТРИ обучающей части: калибратор и подогнанный
        # порог не видят ни одной строки внешнего фолда, даже через модель.
        inner = fi[trn]
        p_in = np.zeros(int(trn.sum()))
        Xtr, ytr = Xi[trn], yi[trn]
        for g in np.unique(inner):
            a, b = inner != g, inner == g
            if b.sum() == 0 or ytr[a].sum() == 0:
                continue
            p_in[b] = gbm_clf().fit(Xtr[a], ytr[a]).predict_proba(Xtr[b])[:, 1]

        p_pla[te], slope = platt(p_in, ytr, p_raw[te])
        slopes.append(slope)
        p_iso[te] = isotonic(p_in, ytr, p_raw[te])

        for key, pp in [("сырые+plug", p_raw), ("Платт+plug", p_pla), ("изотон+plug", p_iso)]:
            t = plugin_threshold(pp[te])
            out[key][te] = pp[te] >= t
            thrs[key].append(t)
        t = fitted_threshold(p_in, ytr)
        out["сырые+подогн"][te] = p_raw[te] >= t
        thrs["сырые+подогн"].append(t)
        t = fitted_threshold(p_raw[te], yi[te])          # потолок, не для подачи
        out["оракул"][te] = p_raw[te] >= t
        thrs["оракул"].append(t)
        print(f"      {tag} фолд {f}: n_те {int(te.sum())}, {time.time()-t0:.0f} с, "
              f"пороги plug {thrs['сырые+plug'][-1]:.2f} / Платт "
              f"{thrs['Платт+plug'][-1]:.2f} / подогн {thrs['сырые+подогн'][-1]:.2f}"
              f", наклон Платта {slope:+.3f}", flush=True)
    return out, p_raw, p_pla, p_iso, thrs, slopes


def boot(y, a, b, n=3000, seed=0):
    """Парный бутстрап разности MCC по молекулам. Решения фиксированы."""
    rng = np.random.default_rng(seed)
    d = np.empty(n)
    for i in range(n):
        k = rng.integers(0, len(y), len(y))
        d[i] = matthews_corrcoef(y[k], b[k]) - matthews_corrcoef(y[k], a[k])
    return float(np.mean(d)), float(np.quantile(d, 0.025)), float(np.quantile(d, 0.975))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--out", default=RES + "preds/tdicalib.json")
    a = ap.parse_args()

    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tdi = pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv").set_index("Molecule_Name")
    keep = rows.Molecule_Name.isin(tdi.index).to_numpy()
    T = tdi.loc[rows.Molecule_Name[keep]].reset_index()
    Xt = X[keep]
    smi = list(rows.SMILES[keep])

    recs = []
    for seed in [int(s) for s in a.seeds.split(",")]:
        fold, ncl = butina_folds(smi, seed=seed)
        print(f"\n=== сид {seed}: {ncl} кластеров на {len(smi)} строках TDI ===", flush=True)
        for c in TDI_CYPS:
            lab = T[f"{c}_is_TDI"]
            m = lab.notna().to_numpy()
            yi = lab[m].astype(int).to_numpy()
            Xi, fi = Xt[m], fold[m]
            print(f"  {c}: n {len(yi)}, положительных {int(yi.sum())} "
                  f"({yi.mean():.3f})", flush=True)
            out, p_raw, p_pla, p_iso, thrs, slopes = run_enzyme(Xi, yi, fi, c)
            if min(slopes) <= 0:
                print(f"    ОТРИЦАТЕЛЬНЫЙ наклон Платта в {sum(s_ <= 0 for s_ in slopes)} "
                      f"фолде(ах): калибровка там переворачивает порядок", flush=True)

            # Платт монотонен ВНУТРИ фолда, но карта у каждого фолда своя, поэтому
            # объединённый AUC после калибровки --- уже другая величина: значения из
            # разных фолдов стали несравнимы. Порядок проверяется там, где он и обещан.
            within = [(roc_auc_score(yi[fi == f], p_raw[fi == f]),
                       roc_auc_score(yi[fi == f], p_pla[fi == f]))
                      for f in range(5) if (fi == f).sum() and 0 < yi[fi == f].mean() < 1]
            dmax = max(abs(x - y) for x, y in within)
            print(f"    AUC сырых {roc_auc_score(yi, p_raw):.4f} (объединённый); "
                  f"внутри фолдов Платт совпадает с сырым до {dmax:.2e}", flush=True)
            if dmax > 1e-9:
                print(f"    порядок внутри фолда сдвинулся на {dmax:.2e} --- смотри "
                      f"наклоны: {' '.join(f'{s_:+.2f}' for s_ in slopes)}", flush=True)
            for tag, pp in [("сырые", p_raw), ("Платт", p_pla), ("изотон", p_iso)]:
                print(f"    {tag:8s} E[p] {pp.mean():.3f} против истинной {yi.mean():.3f}"
                      f"   Брайер {brier_score_loss(yi, pp):.4f}", flush=True)
            base = out["сырые+plug"]
            for key in out:
                mcc = matthews_corrcoef(yi, out[key])
                pos = int(out[key].sum())
                line = (f"    {key:14s} MCC {mcc:.4f}  положительных {pos} "
                        f"({pos/len(yi):.3f})")
                if key != "сырые+plug":
                    mu, lo, hi = boot(yi, base, out[key], seed=seed)
                    line += f"   против подачи {mu:+.4f} [{lo:+.4f}, {hi:+.4f}]"
                print(line, flush=True)
                recs.append({"seed": seed, "cyp": c, "arm": key, "mcc": float(mcc),
                             "pos_rate": pos / len(yi), "n": int(len(yi)),
                             "thr": thrs[key], "platt_slopes": slopes})
        json.dump(recs, open(a.out, "w"))

    df = pd.DataFrame(recs)
    print("\n\nMCC по плечам (строки --- плечи, столбцы --- ферменты, среднее по сидам)")
    print(df.pivot_table(index="arm", columns="cyp", values="mcc").round(4))
    print("\nдоля положительных")
    print(df.pivot_table(index="arm", columns="cyp", values="pos_rate").round(3))
    print("\nзнак против подачи по клеткам (сид x фермент):")
    piv = df.pivot_table(index=["seed", "cyp"], columns="arm", values="mcc")
    for arm in piv.columns:
        if arm == "сырые+plug":
            continue
        d = piv[arm] - piv["сырые+plug"]
        print(f"  {arm:14s} среднее {d.mean():+.4f}, положительных "
              f"{int((d > 0).sum())} из {len(d)}")
    print(f"\nсохранено: {a.out}", flush=True)


if __name__ == "__main__":
    main()
