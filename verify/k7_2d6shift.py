"""Two measurements of the test-set shift disagree on CYP2D6. Which is wrong, and why.

The disagreement, as it stood. src/reweight.py reasons from how the test was built: the
organisers picked anchors from the upper percentiles of our activity distribution and added
structural neighbours, which puts the shift of the label mean somewhere between +0.1 and +0.6.
verify/k5_shift.py measures instead: it runs one model over both sets and compares the
distributions of its own output. On three enzymes the two roughly agree in direction. On
CYP2D6 they do not - the construction argument gives +0.19 and the measurement gives -0.19.

The first thing to check is whether they are measuring the same quantity, and they are not.
reweight.py's anchor percentiles are hard-coded at 93 / 98 / 62 / 90, and its own text says
why 2D6 is 62: that enzyme took no part in anchor selection. Its number is the internal
CONTROL - an enzyme not used to pick the anchors should show no enrichment, and it does not.
It was never an estimate of 2D6's shift, and the +0.1 to +0.6 bracket was derived from the
three enzymes that were used. Applying that bracket to 2D6 is a category error, and this
script exists because the document made it.

That leaves the measurement, which needs its own scrutiny before it can be believed. It is
tight: the shift of -0.190 has a bootstrap interval of [-0.234, -0.145] over the 750 test
compounds, nowhere near zero. But a model can shift its own output without the labels
shifting, simply by meeting chemistry it has not seen, so the number needs a mechanism.

The mechanism is in the biochemistry section. CYP2D6 recognises its substrates through a salt
bridge from an active-site aspartate and glutamate to a protonated basic nitrogen; the other
three bind by lipophilicity, planarity or an anion. So a set depleted in basic amines should
have genuinely lower CYP2D6 activity - and should, if anything, look slightly MORE active on
the other three, because their actives are the compounds a basic centre displaces.

Both halves of that are checkable without any label from the test set, and both are checked
below: whether the test really is depleted, and whether the depletion moves activity in the
direction and roughly the amount observed.

Reads data/feats.npz, data/test_feats.npz and data/test_pred.npz. The last two are written by
verify/k5_shift.py, so run that first. Writes nothing.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()

import json

import numpy as np
import pandas as pd

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
# WITHDRAWN. b = 0.789 / 0.898 / 0.743 / 0.999 is the slope of y on yhat - the share of the
# OUTPUT deviation that is real, which is regression dilution. The share of an INPUT shift
# that reaches the output is the other slope, Cov/Var(y) = R^2/b = 0.29 / 0.41 / 0.22 / 0.59.
# Dividing by b gives neither, so the numbers it produced (+0.018 / +0.164 / -0.256 / +0.437)
# had no reading under which they were right, and the claim that they were a lower bound on
# |delta| was unfounded. The raw shift is printed below without any division; converting it
# to a label shift needs a propagation factor that is not one number, because it depends on
# which axis the composition moved along.
# The 2D6-relevant chemistry, in the order a chemist would ask about it.
KEYS = ["pharm_2d6", "is_base_74", "n_basicN_ali", "ph74_n_cation",
        "n_amine_tert", "n_piperidine", "pka_max_basic"]


def main():
    for f in ("test_feats.npz", "test_pred.npz"):
        if not _pl.Path(D + f).exists():
            raise SystemExit(f"нет data/{f} - сначала verify/k5_shift.py")

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    M = np.load(D + "feats.npz")["MECH"].astype(float)
    MT = np.load(D + "test_feats.npz")["MECH"].astype(float)
    PT = np.load(D + "test_pred.npz")
    oof = json.load(open(RES + "preds/oof.json"))

    print("=" * 92)
    print("1. Насколько прочно измерение сдвига предсказаний")
    print("=" * 92)
    rng = np.random.default_rng(0)
    shift = {}
    print(f"{'фермент':8s} {'сред OOF':>9s} {'сред тест':>10s} {'сдвиг':>8s} "
          f"{'95% интервал (по молекулам)':>28s}")
    for c in CYPS:
        p = np.asarray(oof[f"FP+DESC+MECH|{c}"], float)
        pt = PT[c]
        d = float(pt.mean() - p.mean())
        shift[c] = d
        bs = np.array([rng.choice(pt, len(pt), True).mean()
                       - rng.choice(p, len(p), True).mean() for _ in range(4000)])
        lo, hi = np.percentile(bs, [2.5, 97.5])
        note = "" if lo * hi > 0 else "   ноль внутри"
        print(f"{c:8s} {p.mean():9.3f} {pt.mean():10.3f} {d:+8.3f} "
              f"{'[' + f'{lo:+.3f}' + ',' + f'{hi:+.3f}' + ']':>28s}{note}")

    print("\n  Интервал выше --- по молекулам, и он занижен: тест состоит из якорей с\n"
          "  аналогами, а не из независимых соединений. При кластеризации теста порогом 0.50\n"
          "  выходит 172 группы на 750 соединений, и кластерный бутстрап даёт design effect\n"
          "  по дисперсии от 4.4 до 6.7, то есть ошибку вдвое-втрое шире:\n"
          "    1A2 [-0.111,+0.124]   2C9 [+0.025,+0.257]   2D6 [-0.281,-0.093]   "
          "3A4 [+0.295,+0.576]\n"
          "  На 2D6 ноль не накрывается и там, и это единственный фермент, где так.")

    print()
    print("=" * 92)
    print("2. Обеднён ли тест тем, что узнаёт CYP2D6")
    print("=" * 92)
    print(f"{'признак':18s} {'обучение':>10s} {'тест':>9s} {'отношение':>10s} {'z':>7s}")
    for k in KEYS:
        i = mn.index(k)
        a, b = M[:, i], MT[:, i]
        se = np.sqrt(a.var() / len(a) + b.var() / len(b))
        z = (b.mean() - a.mean()) / max(se, 1e-12)
        print(f"{k:18s} {a.mean():10.3f} {b.mean():9.3f} "
              f"{(b.mean() / a.mean() if a.mean() else np.nan):10.2f} {z:+7.1f}")

    print()
    print("=" * 92)
    print("3. Двигает ли это активность туда, куда наблюдается")
    print("=" * 92)
    i = mn.index("is_base_74")
    b, bt = M[:, i] > 0.5, MT[:, i] > 0.5
    fte = float(bt.mean())
    # Доля оснований берётся ВНУТРИ МАСКИ ФЕРМЕНТА, а не по всем 4905 строкам. Это не
    # придирка: маски сильно различаются по составу, потому что кривые доза-эффект снимали
    # по результатам скрининга, и в маске CYP2D6 оснований 0.355 против 0.11-0.12 у трёх
    # остальных. Считать разность долей от общей доли 0.175 значит сравнивать тест не с той
    # выборкой, на которой измерен контраст активности. Первая версия этого скрипта так и
    # делала и занизила объяснённую долю сдвига 2D6 вчетверо.
    print(f"доля оснований при pH 7.4 на тесте: {fte:.3f}\n")
    print(f"{'фермент':8s} {'в маске':>8s} {'основания':>10s} {'прочие':>9s} {'разрыв':>8s} "
          f"{'ожидаемо':>10s} {'измерено':>10s} {'доля':>6s}")
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy()
        bb = b[m]
        ftr = float(bb.mean())
        gap = float(y[bb].mean() - y[~bb].mean())
        exp = (fte - ftr) * gap
        frac = exp / shift[c] * 100 if shift[c] else float("nan")
        print(f"{c:8s} {ftr:8.3f} {y[bb].mean():10.3f} {y[~bb].mean():9.3f} {gap:+8.3f} "
              f"{exp:+10.3f} {shift[c]:+10.3f} {frac:5.0f}%")

    print("""
Что из этого следует.

CYP2D6 --- единственный фермент, у которого основания АКТИВНЕЕ прочих соединений. У
остальных трёх наоборот. Это ровно то, чего требует биохимия раздела о ферментах: 2D6 держит
лиганда солевым мостиком к протонированному азоту, а три других узнают липофильность,
плоскостность или анион. Тест обеднён основаниями примерно вдвое, и поэтому один и тот же
сдвиг состава опускает активность на 2D6 и слегка поднимает её на остальных.

Знак совпадает на всех четырёх ферментах, а на 2D6 сходится и величина: один бинарный
признак объясняет около трёх четвертей сдвига. На трёх остальных он объясняет единицы
процентов, и там основную часть даёт обогащение активными, ради которого набор и собирали.

Разница между ферментами тут не случайна и стоит того, чтобы её назвать. В маске CYP2D6
оснований 0.355, а в масках остальных трёх --- от 0.114 до 0.122. Кривые доза-эффект снимали
по результатам скрининга, поэтому у 2D6 в размеченной части сидят преимущественно те
соединения, которые он и узнаёт. Тест же размечен целиком, и доля оснований в нём 0.104.
Отсюда и разрыв в четверть: тест разбавлен именно тем, чем маска 2D6 обогащена.

Практический вывод один и он существенный: сдвиг маргинали НЕ ОДИН на все ферменты. Вилка
от +0.1 до +0.6 выведена из трёх ферментов, по которым отбирали якоря, и к CYP2D6 не
относится --- у него сдвиг отрицательный. Подгонять постобработку под общий сдвиг значит
подгонять её для 2D6 в заведомо неверную сторону.""")


if __name__ == "__main__":
    main()
