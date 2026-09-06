"""Четвёртый маршрут к сдвигу delta: прогнать саму процедуру отбора. Модели в контуре нет.

Три существующих маршрута оценивают сдвиг метки на тесте косвенно и все три опираются на модель
или на распределение её выходов: перевзвешивание маргинала (`src/reweight.py`, +0.4 при delta=0 и
+0.9 при delta=0.5), перцентили якорей (+1.05, ВЕРХНЯЯ граница, потому что соседи якорей выбраны
по сходству и регрессируют к среднему) и сдвиг предсказаний (`verify/k5_shift.py`, +0.09, НИЖНЯЯ
граница, потому что интерполирующая модель проносит лишь часть входного сдвига). Докстринг
`src/submit.py` объявляет живым диапазон +0.1..+0.6 и называет `--shrink` единственным ещё
открытым решением.

Четвёртый маршрут не оценивает сдвиг --- он его ВОСПРОИЗВОДИТ. Процедура отбора опубликована:
двадцать пять лучших по CYP1A2, двадцать пять по CYP2C9, двадцать пять по CYP3A4, к каждому якорю
ближайшие соседи. Её можно прогнать на обучающем наборе, где метки известны, и прямо измерить, как
сдвинулось распределение получившегося псевдотеста. Ни одной подгонки, ни одного предсказания.

**Контроль, без которого число не читается.** Расширение по соседям само по себе двигает маргинал:
у любой молекулы соседи похожи на неё, а плотная область набора отличается от разрежённой. Поэтому
рядом с настоящей процедурой прогоняется та же самая, но с ВЫБОРОМ ЯКОРЕЙ СЛУЧАЙНО. Разность двух
арм и есть вклад отбора по потентности; сырой сдвиг псевдотеста его завышает.

**Что здесь заведомо смещено, и обе стороны вниз.** Псевдотест выходит меньше 750 уникальных из-за
пересечения окрестностей, так что якоря в нём перевешены; и пул соседей --- наши 4905 обучающих, а
не вся библиотека организаторов, поэтому соседи ближе и обогащение сильнее. Число соседей на якорь
поэтому варьируется до размера, дающего 750 уникальных.

**И одна клетка этим маршрутом НЕ решается, что важнее прочего.** Дефолт на CYP2D6 отрицательный
(-0.5) и обоснован химией: тест несёт втрое меньше оснований при pH 7.4, чем маска меток 2D6, а
CYP2D6 --- единственный фермент, связывающий через солевой мостик. Здесь проверяется, воспроизводит
ли обеднение основаниями сама процедура отбора; если да, то это её следствие, а не независимый факт
о тесте. Но если симуляция недобирает обеднение, знак этой клетки она не решает, и файл это скажет.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES

import argparse, json
import numpy as np
import pandas as pd
from rdkit import DataStructs
from cypsplit import fingerprints

CY = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
SEL = ["CYP1A2", "CYP2C9", "CYP3A4"]     # три фермента, по которым отбирали якоря


def simulate(Y, fps, sel_idx, n_anchor, n_nb, rng, jitter, random_anchors=False):
    """Один розыгрыш процедуры. Возвращает индексы псевдотеста."""
    chosen = []
    for c in SEL:
        y = Y[c]
        obs = np.where(np.isfinite(y))[0]
        if random_anchors:
            chosen += list(rng.choice(obs, n_anchor, replace=False))
            continue
        # верхушка с джиттером: берём из top-(n_anchor*jitter) случайные n_anchor
        top = obs[np.argsort(-y[obs])[:int(n_anchor * jitter)]]
        chosen += list(rng.choice(top, n_anchor, replace=False))
    chosen = list(dict.fromkeys(chosen))
    out = set(chosen)
    for a in chosen:
        sims = np.asarray(DataStructs.BulkTanimotoSimilarity(fps[a], fps))
        sims[a] = -1
        out.update(np.argsort(-sims)[:n_nb].tolist())
    return np.array(sorted(out)), np.array(chosen)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--reps", type=int, default=30)
    ap.add_argument("--nb", default="10,16,24")
    ap.add_argument("--out", default=RES + "preds/designsim.json")
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = {c: tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CY}
    mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    M = np.load(D + "feats.npz")["MECH"]
    Mt = np.load(D + "test_feats.npz")["MECH"]
    ib, ifp = mn.index("is_base_74"), mn.index("frac_prot_74")
    print("отпечатки", flush=True)
    fps = fingerprints(list(rows.SMILES))

    base = {c: (np.nanmedian(Y[c]), np.nanmean(Y[c]), np.nanstd(Y[c])) for c in CY}
    d6 = np.isfinite(Y["CYP2D6"])
    print(f"опора: is_base_74 обучающий {M[:, ib].mean():.3f}, маска 2D6 {M[d6, ib].mean():.3f}, "
          f"ТЕСТ {Mt[:, ib].mean():.3f}\n", flush=True)

    recs = []
    for nb in [int(x) for x in a.nb.split(",")]:
        for tag, rand in [("отбор по потентности", False), ("КОНТРОЛЬ: случайные якоря", True)]:
            rng = np.random.default_rng(0)
            sizes, shifts, sds, bases, fprot = [], {c: [] for c in CY}, {c: [] for c in CY}, [], []
            for _ in range(a.reps):
                idx, _anch = simulate(Y, fps, None, 25, nb, rng, 3.0, rand)
                sizes.append(len(idx)); bases.append(M[idx, ib].mean()); fprot.append(M[idx, ifp].mean())
                for c in CY:
                    v = Y[c][idx]; v = v[np.isfinite(v)]
                    if len(v) < 30: continue
                    shifts[c].append(np.median(v) - base[c][0])
                    sds[c].append(v.std() / base[c][2])
            print(f"--- соседей {nb}, {tag}: псевдотест {np.mean(sizes):.0f} уникальных", flush=True)
            print(f"    is_base_74 {np.mean(bases):.3f} (тест 0.104, маска 2D6 {M[d6, ib].mean():.3f}), "
                  f"frac_prot_74 {np.mean(fprot):.3f} (тест {Mt[:, ifp].mean():.3f})", flush=True)
            for c in CY:
                if not shifts[c]: continue
                print(f"      {c}: сдвиг медианы {np.mean(shifts[c]):+.3f} "
                      f"[{np.quantile(shifts[c],.1):+.3f}, {np.quantile(shifts[c],.9):+.3f}]   "
                      f"отношение sd {np.mean(sds[c]):.3f}", flush=True)
                recs.append({"nb": nb, "arm": tag, "cyp": c, "shift": float(np.mean(shifts[c])),
                             "lo": float(np.quantile(shifts[c], .1)), "hi": float(np.quantile(shifts[c], .9)),
                             "sdratio": float(np.mean(sds[c])), "size": float(np.mean(sizes)),
                             "is_base": float(np.mean(bases))})
    df = pd.DataFrame(recs)
    print("\n\nСДВИГ МЕДИАНЫ: отбор минус контроль (вклад именно отбора по потентности)")
    p = df.pivot_table(index=["nb", "cyp"], columns="arm", values="shift")
    p["чистый вклад отбора"] = p["отбор по потентности"] - p["КОНТРОЛЬ: случайные якоря"]
    print(p.round(3).to_string())
    print("\nдефолты submit.py: 1A2 0.0, 2C9 +0.3, 2D6 -0.5, 3A4 +0.7")
    json.dump(recs, open(a.out, "w"))
    print(f"\nсохранено: {a.out}", flush=True)


if __name__ == "__main__":
    main()
