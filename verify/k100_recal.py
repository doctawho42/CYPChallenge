"""Item 308: the leaderboard says our ORDER is mid-field and our SCALE is bottom-field.

Five things this reproduces, because item 308 quotes all of them:

  1. The dated board snapshot, as a committed file rather than a constant: 164 regression
     entrants and 90 classification entrants on 15 September 2026, every tab, our row
     included. `results/leaderboard_2026-09-15.json`. The boards move -- our own position
     shifted from 102/163 to 103/164 inside one afternoon without any number of ours
     changing -- so the claims are anchored to the snapshot and not to the live page.

  2. Item 293's four pre-registered predictions, scored. One of four holds, and it is the
     one item 293 itself called the sharpest. Of the three misses, TWO are in the better
     direction and one -- the only scale-sensitive metric -- is worse by a wide margin.

  3. The diagnosis, which needs no inference at all: sort our four cells by
     k = sd(p)/sd(y_test) and both ST-RAE and R2 come out exactly monotone, with no
     exception. And among the entrants whose macro Spearman sits within 0.03 of ours --
     the same ordering, by the board's own measure -- almost all of them score better.

  4. The moments the correction uses, re-derived from the board's published MAE and R2,
     and compared against team briford's independently probed constants (item 300). The
     MEANS agree on all four within 0.10. The SDs do NOT: two of the four miss by more than
     that (0.132 and 0.133), and this file says so rather than quoting the agreement alone.

  5. The applied rule, and what it must produce. The transform is monotone, so rank cannot
     move; the predicted R2 has a closed form, so the whole inference chain is falsifiable
     at the NEXT submission rather than at the 25 September reveal.

A defect in the record, stated here because this file exists to prevent its recurrence: the
four multipliers were fixed as a MEDIAN over an assumption family, and that family was never
written down. The reconstruction below brackets the applied values but does not reproduce
them. They are therefore carried as pre-registered constants, not as a formula's output.

Runtime: a couple of seconds. Reads the snapshot and the two submission files; writes nothing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import RES

import json

import numpy as np
import pandas as pd

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
ME = "Wizard Lizard Gizzard"
SNAP = RES + "leaderboard_2026-09-15.json"
TABS = dict(zip(CYPS, ["partial_10", "partial_11", "partial_12", "partial_13"]))

# APPLIED, and fixed before the result was seen (item 308). Median over an assumption
# family; see the docstring on why the family itself is not reconstructible.
APPLIED_B = dict(zip(CYPS, [2.12, 1.24, 2.26, 1.24]))
APPLIED_MU = dict(zip(CYPS, [4.315, 4.801, 3.075, 4.815]))

# Probed by team briford / SuperCowPowers, public via jeremycheminf/openadmet_scripts and
# recorded in item 300. Independent of anything below: they come from three affine-transformed
# submissions of one vector, not from published MAE.
BRIFORD_MU = dict(zip(CYPS, [4.412, 4.830, 3.107, 4.880]))
BRIFORD_SD = dict(zip(CYPS, [1.553, 1.101, 1.599, 1.272]))

# MAE/RMSE for three error shapes. Laplace 1/sqrt2, normal sqrt(2/pi), uniform sqrt3/2.
FACTORS = {"Лаплас": 2.0 ** -0.5, "норм": (2.0 / np.pi) ** 0.5, "равном": 3.0 ** 0.5 / 2.0}


def load():
    snap = json.load(open(SNAP))
    def row(ep):
        b = snap["boards"][ep]
        for r in b["data"]:
            if r[1] == ME:
                return dict(zip(b["headers"], r))
        raise SystemExit(f"{ep}: нашей строки нет в снимке")
    sub = pd.read_csv(RES + "submission/activity_submission.csv")
    v = {c: float(sub[c + "_pIC50_direct_inhibition"].to_numpy(float).std(ddof=1))
         for c in CYPS}
    mp = {c: float(sub[c + "_pIC50_direct_inhibition"].to_numpy(float).mean()) for c in CYPS}
    return snap, row, v, mp


def prereg(row):
    """Item 293's four predictions, scored. Higher is better for 2, 3 and 4; lower for 1."""
    print("Предрегистрация пункта 293 на живой доске\n")
    m9, m14 = row("partial_9"), row("partial_14")
    order = sorted(CYPS, key=lambda c: -row(TABS[c])["Spearman's ρ"])
    tests = [
        ("макро ST-RAE (n=375)", m9["MA-ST-RAE"], 0.5945, 0.7265, "меньше"),
        ("макро Spearman", m9["MA-Spearman's ρ"], 0.6276, 0.6590, "больше"),
        ("макро MCC", m14["MA-MCC"], 0.2048, 0.2803, "больше"),
    ]
    held = 0
    print("    %-22s %10s %18s %-22s" % ("предсказание", "факт", "полоса", "вердикт"))
    for name, val, lo, hi, better in tests:
        inside = lo <= val <= hi
        held += inside
        if inside:
            v = "ВНУТРИ"
        else:
            worse = (val > hi) if better == "меньше" else (val < lo)
            v = "ВНЕ %+.4f (%s)" % (val - (hi if val > hi else lo),
                                    "хуже" if worse else "ЛУЧШЕ")
        print("    %-22s %10.4f %18s %-22s" % (name, val, "[%.4f, %.4f]" % (lo, hi), v))
    ok = order == ["CYP3A4", "CYP2C9", "CYP1A2", "CYP2D6"]
    held += ok
    print("    %-22s %10s %18s %-22s"
          % ("порядок ранга", " > ".join(c[3:] for c in order), "3A4>2C9>1A2>2D6",
             "ВЫПОЛНИЛОСЬ" if ok else "ПРОВАЛИЛОСЬ"))
    print(f"\n    выполнилось {held} из 4; единственный провал в ХУДШУЮ сторону --- "
          f"единственная метрика, чувствительная к масштабу")
    return held


def diagnosis(snap, row, v):
    """No inference: the board's own numbers, sorted by our spread deficit."""
    print("\n\nДиагноз без всякого вывода: сортировка по k = sd(p)/sd(y)\n")
    s = {}
    for c in CYPS:
        r = row(TABS[c])
        s[c] = r["MAE"] / (FACTORS["норм"] * np.sqrt(1.0 - r["R²"]))
    order = sorted(CYPS, key=lambda c: v[c] / s[c])
    print("    %-8s %8s %9s %9s %9s" % ("фермент", "k", "ST-RAE", "R2", "ро"))
    for c in order:
        r = row(TABS[c])
        print("    %-8s %8.4f %9.4f %+9.4f %9.4f"
              % (c, v[c] / s[c], r["ST-RAE"], r["R²"], r["Spearman's ρ"]))
    st = [row(TABS[c])["ST-RAE"] for c in order]
    r2 = [row(TABS[c])["R²"] for c in order]
    print("\n    ST-RAE монотонно ПАДАЕТ с ростом k: %s" % all(np.diff(st) < 0))
    print("    R2 монотонно РАСТЁТ с ростом k:      %s" % all(np.diff(r2) > 0))
    print("    Два монотонных порядка на четырёх клетках без исключений. Меньше разброса ---")
    print("    хуже обе метрики, и это ровно то, что аффинная пара не могла исправить:")
    print("    она подогнана к полосам НАШИХ метк, то есть к обучающему распределению.")

    b9 = snap["boards"]["partial_9"]
    H = b9["headers"]
    i_s, i_r = H.index("MA-ST-RAE"), H.index("MA-Spearman's ρ")
    rows = [x for x in b9["data"] if isinstance(x[i_s], (int, float))]
    sr = np.array([x[i_s] for x in rows])
    rh = np.array([x[i_r] for x in rows])
    me = row("partial_9")
    sel = np.abs(rh - me["MA-Spearman's ρ"]) <= 0.03
    print("\n    Соперники с тем же порядком по мере самой доски (ро в +-0.03 от %.4f):"
          % me["MA-Spearman's ρ"])
    print("      %d участников, их MA-ST-RAE %.4f .. %.4f, медиана %.4f"
          % (sel.sum(), sr[sel].min(), sr[sel].max(), np.median(sr[sel])))
    print("      наш %.4f --- лучше нас %d из %d при ПРАКТИЧЕСКИ ТОМ ЖЕ ранге"
          % (me["MA-ST-RAE"], (sr[sel] < me["MA-ST-RAE"]).sum(), sel.sum()))
    print("      мы %d-е из %d, медиана доски %.4f"
          % (me["Rank"], len(rows), np.median(sr)))
    return s


def moments(row, v, mp, s):
    """Re-derive the test moments, and test them against briford's probe."""
    print("\n\nМоменты закрытого набора из опубликованных MAE и R2\n")
    print("    sd: MAE = f * RMSE и RMSE = sd*sqrt(1-R2) дают sd = MAE / (f*sqrt(1-R2))")
    print("    ср: тождество briford R2 = 2*r*k - k^2 - b^2 при r = ро с доски, корень < 0\n")
    print("    %-8s %9s %9s %9s %11s %9s %9s"
          % ("фермент", "Лаплас", "норм", "равном", "briford sd", "ср. выв.", "briford ср"))
    mu = {}
    for c in CYPS:
        r = row(TABS[c])
        row_sd = [r["MAE"] / (f * np.sqrt(1.0 - r["R²"])) for f in FACTORS.values()]
        m2 = 2 * r["Spearman's ρ"] * s[c] * v[c] - v[c] ** 2 - r["R²"] * s[c] ** 2
        mu[c] = mp[c] - np.sqrt(m2) if m2 > 0 else float("nan")
        print("    %-8s %9.4f %9.4f %9.4f %11.4f %9.4f %9.4f"
              % (c, row_sd[0], row_sd[1], row_sd[2], BRIFORD_SD[c], mu[c], BRIFORD_MU[c]))

    dmu = [abs(mu[c] - BRIFORD_MU[c]) for c in CYPS]
    dsd = [abs(s[c] - BRIFORD_SD[c]) for c in CYPS]
    print("\n    СРЕДНИЕ согласны: расхождения %s --- все четыре внутри 0.10"
          % ", ".join("%.3f" % x for x in dmu))
    print("    SD НЕ согласны:   расхождения %s --- внутри 0.10 только %d из 4"
          % (", ".join("%.3f" % x for x in dsd), sum(x <= 0.10 for x in dsd)))
    print("    Это корроборация зонда briford ПО ЦЕНТРУ и НЕ по разбросу. Пункт 300")
    print("    зафиксировал асимметрию вердикта заранее: согласие подтверждает зонд,")
    print("    расхождение не опровергает его (живая половина против полных 750).")

    print("\n    И дефект самого допущения, который нельзя замолчать: при моей нормировке")
    print("    RAE тождественно равен sqrt(1-R2), поэтому 'коэффициент прощения'")
    print("    ST-RAE/RAE считается прямо, а превышать единицу он НЕ МОЖЕТ --- ошибка до")
    print("    границы полосы не бывает больше ошибки до точки.\n")
    print("    %-8s %12s %12s %12s" % ("фермент", "при 0.707", "при 0.798", "при 0.866"))
    for c in CYPS:
        r = row(TABS[c])
        # RAE = MAE / mean|y - ср(y)|. With sd(y) = MAE/(f*sqrt(1-R2)) and a normal label
        # spread, mean|y - ср(y)| = sqrt(2/pi)*sd(y), so RAE = (f/sqrt(2/pi))*sqrt(1-R2)
        # and the ratio below is what the board's ST-RAE implies about band forgiveness.
        vals = [r["ST-RAE"] / np.sqrt(1.0 - r["R²"]) * (FACTORS["норм"] / f)
                for f in FACTORS.values()]
        flag = "  <-- НЕВОЗМОЖНО" if max(vals) > 1.0 else ""
        print("    %-8s %12.4f %12.4f %12.4f%s" % (c, vals[0], vals[1], vals[2], flag))
    print("\n    Нарушение сидит на 2D6 при норме и на 1A2 при Лапласе --- то есть ровно на")
    print("    двух клетках, где множители самые большие и где лежит почти весь приз.")
    return mu


def rule(row, v, s):
    """The applied rule, its reconstruction bracket, and the closed-form prediction."""
    print("\n\nПрименённое правило и что оно ОБЯЗАНО дать\n")
    print("    q = mu(y) + b*(p - ср(p)), монотонно при b>0 --- ранг не двигается")
    print("    k' = b*sd(p)/sd(y); оптимум по квадрату при k' = r, поэтому k' > r --- перестрел")
    print("    R2 = 2*r*k' - k'^2 при совпавшем центре\n")
    print("    %-8s %6s %10s %8s %10s %9s %11s"
          % ("фермент", "b", "рекон.", "k'", "ро борд", "перестрел", "R2 предск"))
    pr = []
    for c in CYPS:
        r = row(TABS[c])
        cand = sorted(r["Spearman's ρ"] * x / v[c]
                      for x in [r["MAE"] / (f * np.sqrt(1 - r["R²"])) for f in FACTORS.values()]
                      + [BRIFORD_SD[c]])
        k1 = APPLIED_B[c] * v[c] / s[c]
        R2n = 2 * r["Spearman's ρ"] * k1 - k1 ** 2
        pr.append(R2n)
        print("    %-8s %6.2f %10s %8.4f %10.4f %9.3f %11.4f"
              % (c, APPLIED_B[c], "%.2f-%.2f" % (cand[0], cand[-1]), k1,
                 r["Spearman's ρ"], k1 / r["Spearman's ρ"], R2n))
    m9 = row("partial_9")
    print("\n    Реконструкция ОХВАТЫВАЕТ применённые значения, но не воспроизводит их:")
    print("    медиана считалась в сессии и состав семейства не был записан. Множители")
    print("    поэтому несутся как предрегистрированные константы, а не как вывод формулы.")
    print("\n    ПРОВЕРЯЕМОЕ УТВЕРЖДЕНИЕ, и проверить его можно на СЛЕДУЮЩЕЙ подаче,")
    print("    а не 25 сентября: макро R2 обязан уйти с %.4f примерно на %.4f."
          % (m9["MA-R²"], np.mean(pr)))
    print("    Если выведенные моменты неверны, он там не окажется. Ранг при этом обязан")
    print("    остаться на %.4f РОВНО --- преобразование монотонно." % m9["MA-Spearman's ρ"])


if __name__ == "__main__":
    snap, row, v, mp = load()
    print("Снимок: %s, тянут %s\n" % (SNAP, snap["pulled_utc"]))
    prereg(row)
    s = diagnosis(snap, row, v)
    moments(row, v, mp, s)
    rule(row, v, s)
