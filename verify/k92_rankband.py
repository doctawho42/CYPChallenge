"""Полоса выборки для РАНГА на подаваемом плече. Пункт 293 сделал её обязательной.

Зачем именно теперь. До 12 сентября табло утверждало, что «ранг на лидерборде не показывают», и
поэтому предрегистрировать ранговое предсказание было незачем. Пункт 293 это опроверг: pinned
README организаторов говорит, что вторичные метрики (MAE, R^2, Spearman rho, Kendall tau)
публикуются С БУТСТРЭП-ИНТЕРВАЛАМИ. Значит собственная валюта проекта --- ранг, в которой измерены
почти все его результаты, --- впервые становится ВНЕШНЕ проверяемой, и 25 сентября её можно
falsify. Но предсказание без полосы не фальсифицируемо, а полосы на ранге не считал никто: все
четыре предрегистрированные полосы (пункты 246, 253, 255, 256) --- по ST-RAE, а пункт 291/k91 ---
по MCC.

Вход --- `results/preds/oof_submitted.json`: вневыборочные предсказания ПОДАВАЕМОГО плеча
(ансамбль5 + мёртвая зона + SOLO 1A2=пофе+GP+ствол, 2C9/3A4=GP+ствол, сид 0), тот самый файл, по
которому считалась полоса ST-RAE в `k79_bandfix`. Никакого обучения: секунды.

ПОЧЕМУ РАНГ НЕ НАДО ПРОГОНЯТЬ ЧЕРЕЗ АФФИННУЮ ПАРУ. Пара строго возрастающая (lambda > 0), поэтому
Спирмен до и после неё совпадает ТОЧНО --- пункт 277 проверил это прямо, получив 1.0000 на всех
четырёх ферментах. Значит ранг подаваемых предсказаний равен рангу вневыборочных, и полоса,
посчитанная здесь, относится к тому, что уйдёт на лидерборд.

ДВЕ РАЗНЫЕ ВЕЛИЧИНЫ, как в k91, и по той же причине (пункт 276 --- смешение стоило табло):

  1. CI ОРГАНИЗАТОРОВ --- бутстрэп С ВОЗВРАЩЕНИЕМ на своём n, 1000 повторов, сид 0
     (`evaluation/config.py:BOOTSTRAP_SAMPLES`, `evaluation/utils.py:BOOTSTRAP_SEED`). Это то, что
     они напечатают рядом с нашим rho.
  2. ПОЛОСА РАСКРЫТИЯ --- подвыборка БЕЗ возвращения до тестового размера: насколько rho гуляет от
     того, КАКИЕ соединения раскроют. Только она делает предсказание проверяемым.

Размеры. Тест --- 750 соединений, прогнанных против всех четырёх изоформ, поэтому для прямого
ингибирования n=750 на фермент правдоподобен (в отличие от трека TDI, где метка разрежена и
k91 берёт n=375 как более честный). n=375 оставлен как размер живого лидерборда (табло, пункт 256).

ЧЕГО ЭТА ПОЛОСА НЕ ПОКРЫВАЕТ, и это надо сказать до раскрытия. Она --- полоса ВЫБОРКИ: разброс от
того, какие молекулы попали в набор. Она НЕ покрывает сдвиг распределения: тест лежит там, где у
кросс-валидации почти нет массы (пункт 123, chi^2 = 2.838, потолок ESS 26.1%), медиана
ближайшего соседа на тесте 0.587 против <=0.450 на любом переразбиении обучения (пункт 129). Счёт
ВНЕ полосы фальсифицирует конкретное утверждение «тест --- такая же выборка, как наша»; счёт
ВНУТРИ полосы не подтверждает ничего.

Пишет results/logs/k92_rankband.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES

import argparse, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
N_BOOT = 1000          # evaluation/config.py:BOOTSTRAP_SAMPLES
BOOT_SEED = 0          # evaluation/utils.py:BOOTSTRAP_SEED
SIZES = [750, 375]


def rho(y, p):
    return float(spearmanr(y, p).statistic)


def organiser_ci(y, p, rng):
    n = len(y)
    v = np.empty(N_BOOT)
    for b in range(N_BOOT):
        i = rng.integers(0, n, n)
        v[b] = rho(y[i], p[i])
    return v


def reveal_band(y, p, n, rng):
    N = len(y)
    if n >= N:
        return None
    v = np.empty(N_BOOT)
    for b in range(N_BOOT):
        i = rng.choice(N, n, replace=False)
        v[b] = rho(y[i], p[i])
    return v


def summarise(v):
    return {"mean": float(v.mean()), "sd": float(v.std(ddof=1)),
            "lo95": float(np.percentile(v, 2.5)), "hi95": float(np.percentile(v, 97.5))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=RES + "logs/k92_rankband.json")
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(Y)

    C = json.load(open(RES + "preds/oof_submitted.json"))
    print(f"плечо: {C['что']}\n", flush=True)
    P = [np.asarray(v, float) for v in C["P"]]

    out, per_enzyme_draws = {"плечо": C["что"]}, {}
    for e, c in enumerate(CYPS):
        y, p = Y[mask[:, e], e], P[e]
        if len(y) != len(p):
            raise SystemExit(f"{c}: метки {len(y)} против предсказаний {len(p)} --- маска разъехалась")
        point = rho(y, p)
        ci = organiser_ci(y, p, np.random.default_rng(BOOT_SEED))
        rec = {"n": int(len(y)), "rho": point, "ci_organisers": summarise(ci), "reveal": {}}
        print(f"{c}: n={len(y)}  rho={point:+.4f}", flush=True)
        print(f"    CI организаторов (с возвр., n={len(y)}): "
              f"[{rec['ci_organisers']['lo95']:+.4f}, {rec['ci_organisers']['hi95']:+.4f}]  "
              f"sd {rec['ci_organisers']['sd']:.4f}", flush=True)
        draws = {"ci": ci}
        for n in SIZES:
            v = reveal_band(y, p, n, np.random.default_rng(BOOT_SEED + n))
            if v is None:
                print(f"    полоса раскрытия n={n}: пропущена (n >= {len(y)})", flush=True)
                continue
            s = summarise(v); rec["reveal"][str(n)] = s
            draws[str(n)] = v
            print(f"    полоса раскрытия n={n} (без возвр.): "
                  f"[{s['lo95']:+.4f}, {s['hi95']:+.4f}]  sd {s['sd']:.4f}  "
                  f"полуширина {(s['hi95'] - s['lo95']) / 2:.4f}", flush=True)
        out[c] = rec
        per_enzyme_draws[c] = draws

    # Макро: усредняем ПО ФЕРМЕНТАМ внутри каждого ресэмпла. Наборы соединений у ферментов
    # разные, поэтому каждый ресэмплится своим потоком; это независимые розыгрыши, и их
    # поферментное среднее --- корректное распределение макро-ранга.
    print(flush=True)
    for key, label in [("ci", "макро CI организаторов")] + [(str(n), f"макро полоса n={n}") for n in SIZES]:
        vs = [per_enzyme_draws[c][key] for c in CYPS if key in per_enzyme_draws[c]]
        if len(vs) != len(CYPS):
            continue
        m = np.mean(vs, axis=0)
        out[label] = summarise(m)
        s = out[label]
        print(f"{label}: [{s['lo95']:+.4f}, {s['hi95']:+.4f}]  sd {s['sd']:.4f}  "
              f"среднее {s['mean']:+.4f}  полуширина {(s['hi95'] - s['lo95']) / 2:.4f}", flush=True)

    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump(out, open(a.out, "w"), ensure_ascii=False, indent=1)
    print(f"\nсохранено: {a.out}", flush=True)
    print("""
Как читать. Это полоса ВЫБОРКИ по рангу на подаваемом плече, и она сравнима с тем, что организаторы
напечатают, потому что аффинная пара ранг не меняет (пункт 277: Спирмен 1.0000 до и после).

Полосу раскрытия надо сравнивать с РАЗНОСТЯМИ, которые проект заявляет по рангу: траектория
+0.0579 от базовой модели, состав +0.0059, мёртвая зона +0.0197, макро-пол 0.0036. Если полоса
одиночного счёта окажется шире этих разностей --- это не противоречие: лидерборд считает все подачи
на ОДНИХ И ТЕХ ЖЕ молекулах, поэтому поштучный шум в разности сокращается, и сравнивать наши
приросты надо с парным полом, а не с этой полосой. Полоса нужна для другого: чтобы 25 сентября
можно было сказать, лежит ли наш АБСОЛЮТНЫЙ ранг там, где мы его ждали.""")


if __name__ == "__main__":
    main()
