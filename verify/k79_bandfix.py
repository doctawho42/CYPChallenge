"""Полоса фальсификации, пересчитанная на ПОДАВАЕМОЙ руке. Пункт 246 отозван, это замена.

Пункт 246 посчитал полосу на `results/preds/oof_dzens.json`, где рука называется «мёртвая зона
везде» и является ЧЕТЫРЁХЧЛЕННЫМ ансамблем пункта 149 при сиде 0, макро 0.6599. Имя совпало с
подаваемой конфигурацией, состав --- нет.

Проверяя это, нашлась третья конфигурация: `verify/k58_dzsubmit.py` даёт 0.6459 для пяти членов с
проходом во всех, но **не применяет `SOLO`** --- поферментный отбор пункта 218, по которому на
CYP3A4 подаётся ОДИН гауссовский процесс, а не среднее пяти. Значит ни 0.6599, ни 0.6459 не
описывают то, что отправляется, и счёта подаваемой руки не было нигде.

Здесь он считается через `src/submit.py:oof_predictions`, которая применяет `_keep` и потому
собирает ровно ту композицию, что уходит в файл: пять членов с мёртвой зоной на CYP1A2, CYP2C9 и
CYP2D6, один GP на CYP3A4. Аффинная пара --- `src/shrinkchoice.py:fit_apply`, поферментно и по
фолдам, то есть та же, которой меряются все остальные числа проекта.

Полос две, потому что организаторы публикуют два внешних замера: живой лидерборд на ПОЛОВИНЕ
теста, нарезанной по химическим сериям, и промежуточное раскрытие 25 сентября на ПОЛНЫХ 750.

Это ВЫБОРОЧНАЯ полоса, а не предсказательный интервал: она покрывает розыгрыш, но не сдвиг
распределения, который пункты 123, 129 и 147 тремя независимыми маршрутами объявили непроверяемым
изнутри. Счёт вне полосы опровергает названное допущение, а не модель; счёт внутри не подтверждает
ничего.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES

import argparse, json, time
import numpy as np
import pandas as pd

from cypsplit import butina_folds
import submit as SB
from shrinkchoice import fit_apply

CYPS = SB.CYPS


def strae(p, y, lo, hi):
    num = np.maximum(0.0, np.maximum(lo - p, p - hi)).sum()
    mu = y.mean()
    den = np.maximum(0.0, np.maximum(lo - mu, mu - hi)).sum()
    return num / den if den > 0 else np.nan


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--draws", type=int, default=3000)
    ap.add_argument("--out", default=RES + "preds/band.json")
    a = ap.parse_args()

    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)
    fold, _ = butina_folds(list(rows.SMILES))

    cache = RES + "preds/oof_submitted.json"
    if _pl.Path(cache).exists():
        C = json.load(open(cache))
        P = [np.asarray(v, float) for v in C["P"]]
        print(f"вневыборочные предсказания подаваемой руки взяты из {cache}", flush=True)
    else:
        print("вневыборочный проход подаваемой композиции (пять членов + мёртвая зона, SOLO на 3A4)",
              flush=True)
        t0 = time.time()
        P = SB.oof_predictions(X, y, mask, fold, "ансамбль5", None, (LO, HI))
        print(f"  готово за {time.time()-t0:.0f} с", flush=True)
        # СОХРАНЯЕМ. Отсутствие этих предсказаний --- причина, по которой пункт 246 посчитал
        # полосу на чужой руке: счёта подаваемой конфигурации не существовало нигде.
        json.dump({"P": [v.tolist() for v in P],
                   "что": "oof_predictions, ансамбль5, мёртвая зона, SOLO на CYP3A4, сид 0"},
                  open(cache, "w"))
        print(f"  сохранено: {cache}", flush=True)

    # Подаваемое преобразование --- аффинная пара, подогнанная под НАКЛОНЁННУЮ цель при
    # разворачиваемых дельтах, а не обычная пара. Полоса должна относиться к нему.
    SB.LO, SB.HI = LO, HI
    DELTA = [float(x) for x in SB.DELTA_DEFAULT.split(",")]   # ЧИТАЕТСЯ ИЗ ПОДАЧИ
    lams = SB.fit_shrinkage(P, y, mask, DELTA)
    print(f"подаваемые дельты {DELTA}; (lambda, mu, sh) по ферментам: "
          + ", ".join(f"{c[3:]} ({l:.3f}, {m:.2f}, {sh:+.3f})"
                      for c, (l, m, sh) in zip(CYPS, lams)), flush=True)

    per, packs, per_plain = {}, {}, {}
    for e, c in enumerate(CYPS):
        m = mask[:, e]
        w = np.ones(int(m.sum()))
        plain = fit_apply(P[e], LO[m, e], HI[m, e], fold[m], w)
        L_, mu_, sh_ = lams[e]
        depl = L_ * P[e] + (1.0 - L_) * mu_ + sh_
        ok = np.isfinite(plain) & np.isfinite(depl)
        per_plain[c] = float(strae(plain[ok], y[m, e][ok], LO[m, e][ok], HI[m, e][ok]))
        per[c] = float(strae(depl[ok], y[m, e][ok], LO[m, e][ok], HI[m, e][ok]))
        packs[c] = dict(idx=np.where(m)[0][ok], q=depl[ok], y=y[m, e][ok],
                        lo=LO[m, e][ok], hi=HI[m, e][ok])
        print(f"  {c}: подаваемое преобразование {per[c]:.4f}   обычная пара {per_plain[c]:.4f}"
              f"   цена ставки на сдвиг {per[c]-per_plain[c]:+.4f}  (n {int(ok.sum())})", flush=True)
    macro = float(np.mean([per[c] for c in CYPS]))
    macro_plain = float(np.mean([per_plain[c] for c in CYPS]))
    print(f"\nМАКРО: подаваемое {macro:.4f}, обычная пара {macro_plain:.4f}, "
          f"цена ставки {macro-macro_plain:+.4f}", flush=True)
    print(f"\nМАКРО ПОДАВАЕМОЙ РУКИ: {macro:.4f}", flush=True)
    print("  для сравнения: 0.6599 --- четырёхчленная (пункт 149), 0.6459 --- пять членов без "
          "SOLO (пункт 213)", flush=True)

    pos = {c: {v: k for k, v in enumerate(packs[c]["idx"])} for c in CYPS}
    rng = np.random.default_rng(0)
    out = {"macro_oof": macro, "macro_plain": macro_plain, "per": per,
           "per_plain": per_plain, "delta": DELTA, "bands": {}}
    for n in (750, 375):
        draws = []
        for _ in range(a.draws):
            pick = rng.integers(0, len(rows), n)
            vals = []
            for c in CYPS:
                s = np.array([pos[c][j] for j in pick if j in pos[c]])
                d = packs[c]
                vals.append(strae(d["q"][s], d["y"][s], d["lo"][s], d["hi"][s])
                            if len(s) >= 20 else np.nan)
            draws.append(np.nanmean(vals))
        dr = np.array(draws)
        q = np.quantile(dr, [0.025, 0.975])
        out["bands"][n] = dict(mean=float(dr.mean()), sd=float(dr.std()),
                               lo=float(q[0]), hi=float(q[1]))
        print(f"n={n}: среднее {dr.mean():.4f}  sd {dr.std():.4f}  "
              f"95 % [{q[0]:.4f}, {q[1]:.4f}]  полуширина {(q[1]-q[0])/2:.4f}", flush=True)
    json.dump(out, open(a.out, "w"))
    print(f"\nсохранено: {a.out}", flush=True)


if __name__ == "__main__":
    main()
