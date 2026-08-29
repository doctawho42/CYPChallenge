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
# Attenuation of the affine recalibration, from src/reweight.py: a model that interpolates
# passes only part of an input shift into its output, so dividing by it turns a shift in
# predictions into a lower bound on the shift in labels.
ATT = {"CYP1A2": 0.789, "CYP2C9": 0.898, "CYP2D6": 0.743, "CYP3A4": 0.999}
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
          f"{'95% интервал':>19s} {'/ затухание':>12s}")
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
              f"[{lo:+6.3f},{hi:+6.3f}] {d / ATT[c]:+12.3f}{note}")

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
    ftr, fte = float(b.mean()), float(bt.mean())
    print(f"доля оснований при pH 7.4: обучение {ftr:.3f}, тест {fte:.3f}, "
          f"разница {fte - ftr:+.3f}\n")
    print(f"{'фермент':8s} {'основания':>10s} {'прочие':>9s} {'разрыв':>8s} "
          f"{'ожидаемо':>10s} {'измерено':>10s} {'знак':>6s}")
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy()
        bb = b[m]
        gap = float(y[bb].mean() - y[~bb].mean())
        exp = (fte - ftr) * gap
        ok = "сходится" if np.sign(exp) == np.sign(shift[c]) else "РАЗНЫЙ"
        print(f"{c:8s} {y[bb].mean():10.3f} {y[~bb].mean():9.3f} {gap:+8.3f} "
              f"{exp:+10.3f} {shift[c]:+10.3f} {ok:>6s}")

    print("""
Что из этого следует.

CYP2D6 --- единственный фермент, у которого основания АКТИВНЕЕ прочих соединений. У
остальных трёх наоборот. Это ровно то, чего требует биохимия раздела о ферментах: 2D6 держит
лиганда солевым мостиком к протонированному азоту, а три других узнают липофильность,
плоскостность или анион. Тест обеднён основаниями примерно вдвое, и поэтому один и тот же
сдвиг состава опускает активность на 2D6 и слегка поднимает её на остальных.

Знак совпадает на всех четырёх ферментах. Величина --- нет: один бинарный признак объясняет
около пятой части сдвига на 2D6 и малую долю на 2C9 и 3A4, где основную часть даёт
обогащение активными, ради которого набор и собирали. Это оценка направления, а не модель
сдвига, и выдавать её за вторую нельзя.

Практический вывод один и он существенный: сдвиг маргинали НЕ ОДИН на все ферменты. Вилка
от +0.1 до +0.6 выведена из трёх ферментов, по которым отбирали якоря, и к CYP2D6 не
относится --- у него сдвиг отрицательный. Подгонять постобработку под общий сдвиг значит
подгонять её для 2D6 в заведомо неверную сторону.""")


if __name__ == "__main__":
    main()
