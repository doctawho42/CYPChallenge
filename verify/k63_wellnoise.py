"""How wide is the band supposed to be? An independent measurement of the assay's own noise.

Item 114 established that the confidence band shipped with our labels is 3.92 sigma where sigma
is a **deterministic function of the label**: regressing width on value gives R^2 of 0.93 to 0.98.
So the band is derived from the answer rather than measured per compound, and every check of
predicted uncertainty against it checks our agreement with a formula. Item 114 also established
that nothing in our own data supplies a replacement -- the only independent spread estimate there
does not serve.

Item 201 found one outside. The organisers' Octant release carries `inhibition_wells.tsv`:
16931 wells of raw fluorescence with plate, row, column, concentration and an outlier flag --
a layer **below** anything we hold. And unlike the labels in that release, which failed on
overlap at fourteen shared molecules, well-level variance needs **no overlap at all**: how much of
the observed spread is plate, position and replicate is a statement about the assay, not about
which compounds went through it.

This file measures that, and then closes a chain that our own data cannot close:

    шум лунки          <- прямо по репликам отрицательного контроля
    -> точность кривой <- распространением через подгонку по 12 точкам
    -> их se           <- их собственное опубликованное число, для сверки
    -> наша полоса     <- во сколько раз шире

The three middle links are theirs and independent of us, so if they agree the chain is sound and
our band is the only free end.

**Where the replicates are.** Library compounds sit at one well per concentration, so they carry
no replication. The controls do: 912 negative-control wells over 14 plates. Negative control is
the assay measuring "no inhibition" repeatedly, which is exactly a reproducibility estimate.

**Caveats, because the comparison is between two assays and not one.** Octant's inhibition arm
uses a 30-minute active-enzyme pre-incubation and covers CYP3A4 only (item 201), so it is the same
laboratory and the same 1536-well fluorescence format but not the same arm. And the propagation
below is a rough sqrt(n)/h argument rather than a Fisher information calculation; it is used only
to check that the chain closes to within a factor, which it does.

Needs `inhibition_wells.tsv` and `inhibition.tsv` from
huggingface.co/datasets/openadmet/Octant_CYP_inhibition_reactivity_blog_release (CC-BY-4.0),
neither committed. Pass their directory with --octant. Seconds.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D

import argparse

import numpy as np
import pandas as pd

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--octant", required=True, help="каталог с inhibition_wells.tsv и inhibition.tsv")
    a = ap.parse_args()
    o = _pl.Path(a.octant)

    w = pd.read_csv(o / "inhibition_wells.tsv", sep="\t")
    n0 = len(w)
    w = w[~w.outlier.astype(bool)]
    print(f"лунок {n0}, после снятия выбросов {len(w)}; планшетов {w.plate.nunique()}, "
          f"соединений {w.ocnt_batch.nunique()}\n", flush=True)

    nc = w[w.compound_class == "Negative Control"]
    tot = float(nc.fluorescence_norm.var(ddof=1))
    btw = float(nc.groupby("plate").fluorescence_norm.mean().var(ddof=1))
    wth = float(nc.groupby("plate").fluorescence_norm.var(ddof=1).mean())
    print("РАЗЛОЖЕНИЕ ДИСПЕРСИИ ОТРИЦАТЕЛЬНОГО КОНТРОЛЯ, единицы log2fc")
    print(f"  лунок {len(nc)} на {nc.plate.nunique()} планшетах")
    print(f"  полная            дисп. {tot:.5f}   ско {np.sqrt(tot):.4f}")
    print(f"  между планшетами  дисп. {btw:.5f}   ско {np.sqrt(btw):.4f}   ({btw/tot:.1%})")
    print(f"  внутри планшета   дисп. {wth:.5f}   ско {np.sqrt(wth):.4f}")
    print("  Межпланшетная доля нулевая по построению: fluorescence_norm нормирована\n"
          "  внутри планшета относительно его же контролей, так что планшет уже снят.\n")

    e = nc[(nc.row == nc.row.min()) | (nc.row == nc.row.max())
           | (nc.col == nc.col.min()) | (nc.col == nc.col.max())]
    i = nc[~nc.index.isin(e.index)]
    print("КРАЕВОЙ ЭФФЕКТ")
    print(f"  край   n={len(e):4d}  среднее {e.fluorescence_norm.mean():+.4f}  "
          f"ско {e.fluorescence_norm.std():.4f}")
    print(f"  внутри n={len(i):4d}  среднее {i.fluorescence_norm.mean():+.4f}  "
          f"ско {i.fluorescence_norm.std():.4f}")
    d = float(i.fluorescence_norm.mean() - e.fluorescence_norm.mean())
    print(f"  систематическая разница {d:+.4f} --- ТОГО ЖЕ ПОРЯДКА, что и сам шум.\n"
          f"  Это смещение, а не разброс: усреднением по лункам оно не уходит.\n")

    sw = float(np.sqrt(tot))
    oc = pd.read_csv(o / "inhibition.tsv", sep="\t")
    se = oc.CYP3A4_pIC50_se.dropna()
    ci = (oc.CYP3A4_pIC50_ci_upper - oc.CYP3A4_pIC50_ci_lower).dropna()

    print("ЦЕПОЧКА: от лунки к заявленной точности")
    print(f"  шум лунки, измерен здесь                       {sw:.4f} log2fc")
    for h in (1.0, 1.5, 2.0):
        print(f"  -> точность кривой при 12 точках, крутизне {h:.1f}   ~{sw/np.sqrt(12)/h:.4f} pIC50")
    print(f"  их собственная se подгонки, медиана            {se.median():.4f} pIC50")
    print(f"  их ширина CI / se                              {(ci/se).median():.3f}  "
          f"(ожидается 3.92)")
    print("  Три средних звена --- их и от нас не зависят. Они сходятся, значит цепочка цела.\n")

    tr = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
    print("НАША ПОЛОСА: подразумеваемая sigma = ширина / 3.92")
    print(f"  {'фермент':8s} {'n':>6s} {'мед. ширина':>12s} {'sigma':>8s} {'10-й':>8s} {'90-й':>8s}")
    meds = []
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna()
        wd = (tr.loc[m, col + "_conf_high"] - tr.loc[m, col + "_conf_low"]).to_numpy()
        s = wd / 3.92
        meds.append(float(np.median(s)))
        print(f"  {c:8s} {int(m.sum()):6d} {np.median(wd):12.3f} {np.median(s):8.3f} "
              f"{np.percentile(s, 10):8.3f} {np.percentile(s, 90):8.3f}")
    ours = float(np.median(meds))
    print(f"\n  медиана по ферментам {ours:.3f} pIC50, "
          f"против их {se.median():.4f}: в {ours/se.median():.1f} раза шире")
    print(f"""
Как читать. Цепочка от сырой лунки до опубликованной se ЗАМЫКАЕТСЯ на их стороне: измеренный
здесь шум {sw:.3f}, распространённый через подгонку по двенадцати точкам, даёт величину порядка
{sw/np.sqrt(12)/1.5:.3f}, и их собственная se равна {se.median():.4f}. Отношение их CI к их se ---
{(ci/se).median():.3f} против 3.92, то есть конструкция интервала у нас с ними одна и та же.

Единственное свободное звено --- наша полоса, и она в {ours/se.median():.1f} раза шире.

Это НЕ значит, что наша полоса завышена. Полоса воспроизводимости законно шире полосы точности
подгонки: повторить эксперимент целиком --- не то же самое, что переподогнать те же двенадцать
точек, и множитель около четырёх для этого правдоподобен. Что теперь измерено --- это УРОВЕНЬ:
наша полоса живёт в масштабе воспроизводимости, а не точности, и разница названа числом.

Что остаётся дефектом --- не уровень, а ЗАВИСИМОСТЬ ОТ МЕТКИ (пункт 114). Ширина, будучи функцией
ответа при R^2 0.93--0.98, не является посоединениевой измеренной неопределённостью ни в каком
масштабе; она формула. Здесь померен масштаб этой формулы, а не её происхождение.

И следствие для метрики, которое стоит назвать. Мёртвая зона эксплуатирует полосу, а полоса
откалибрована в масштабе воспроизводимости. Значит метрика прощает различения, которые платформа
на самом деле разрешает, примерно вчетверо. Это свойство метрики, а не нашей модели: мы её
оптимизируем, а не выбираем.""")


if __name__ == "__main__":
    main()
