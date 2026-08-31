"""Rank counted only over the pairs the metric can pay for.

Where this came from. The four-seed own-boosting sweep put a pairwise objective against a squared
one on the same trees, and the pairwise arm won overall (0.5692 against 0.5615 of Spearman) while
**losing CYP3A4**: 0.738 against 0.754. That is the wrong shape for a loss function that only ever
discards pairs the bands cannot order, and it points at the criterion rather than the arm.

The argument. ST-RAE, after the affine pair, pays for getting an order right only when the two
bands are disjoint -- if they overlap, no prediction can be penalised for putting the pair either
way round. Spearman against the point label pays for **every** pair. So on an enzyme with many
overlapping bands, rho is partly earned on pairs the competition will never score, and an
objective that deliberately drops those pairs must look worse under rho while being no worse under
the metric. CYP3A4 is exactly that enzyme: item 138 measured 51.4 per cent of its labels below the
instrument's resolution floor, the highest of the four.

Item 80 made rank the criterion because raw ST-RAE was untrustworthy, and that was right. But the
rank it chose counts invisible pairs, and this is the first arm sharp enough to expose it.

The statistic. A Kendall tau restricted to distinguishable pairs,

    tau_vis = (C - D) / V   over ordered pairs with lo_i > hi_j,

where V is how many such pairs exist. It is the same criterion as rho -- order against truth, not
scale -- with the pairs the metric cannot see removed. Nothing else changes.

Pre-registered readings, and both are informative.

  If CYP3A4's apparent loss shrinks or reverses under tau_vis, the loss was an artefact of the
  criterion. Then item 80's rule needs the qualifier that rho over-counts on wide-banded enzymes,
  and every per-enzyme comparison in this file involving CYP3A4 carries that bias.

  If the loss survives under tau_vis, the pairwise objective really does give up CYP3A4 order that
  the metric would have paid for, the trade is real, and the arm is a per-enzyme choice rather
  than an improvement.

Reads whatever prediction files exist under results/preds/. Writes nothing. ~2 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from cypsplit import butina_folds

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]

# Файл -> какие ключи из него брать. Ключи различаются по скриптам, поэтому берётся
# всё, что кончается на имя фермента, а рука вынимается из середины.
FILES = [("oof.json", "FP+DESC+MECH"), ("oof_ownboost.json", None),
         ("oof_pairloss.json", None), ("oof_dead.json", None), ("oof_dead123.json", None),
         ("oof_pool.json", None), ("oof_aux.json", None), ("oof_deadpair.json", None),
         ("oof_ncgc.json", None)]


def tau_visible(y_lo, y_hi, p, chunk=512):
    """Kendall tau по парам с непересекающимися полосами. Возвращает (tau, доля пар)."""
    n = len(p)
    C = Dd = V = 0
    for a in range(0, n, chunk):
        b = min(a + chunk, n)
        vis = y_lo[a:b, None] > y_hi[None, :]      # истина i строго выше j
        if not vis.any():
            continue
        d = p[a:b, None] - p[None, :]
        V += int(vis.sum())
        C += int((vis & (d > 0)).sum())
        Dd += int((vis & (d < 0)).sum())
    if V == 0:
        return np.nan, 0.0
    return (C - Dd) / V, V / (n * (n - 1) / 2)


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())

    print(f"{'фермент':9s} {'строк':>6s} {'видимых пар':>12s} {'доля меток < pC0':>17s}")
    truth = {}
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy()
        lo = tr.loc[m, col + "_conf_low"].to_numpy()
        hi = tr.loc[m, col + "_conf_high"].to_numpy()
        truth[c] = (y, lo, hi)
        _, share = tau_visible(lo, hi, y)
        print(f"{c:9s} {len(y):6d} {100*share:11.1f} % {100*(y <= 4.305).mean():16.1f} %")
    print()

    table = []
    for fn, only in FILES:
        path = RES + "preds/" + fn
        if not os.path.exists(path):
            print(f"(нет {fn}, пропускаю)")
            continue
        blob = json.load(open(path))
        preds = blob.get("preds", blob)
        # Ключи бывают "рука|фермент" и "сид|рука|фермент".
        arms = {}
        for k in preds:
            parts = k.split("|")
            if parts[-1] not in CYPS:
                continue
            arm = "|".join(parts[:-1])
            if only is not None and parts[0] != only:
                continue
            arms.setdefault(arm, {})[parts[-1]] = np.asarray(preds[k], float)
        for arm, d in arms.items():
            if len(d) < len(CYPS):
                continue
            r = {"файл": fn, "рука": arm}
            for c in CYPS:
                y, lo, hi = truth[c]
                p = d[c]
                if len(p) != len(y):
                    r = None
                    break
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
                r[f"{c} tau"] = round(float(tau_visible(lo, hi, p)[0]), 4)
            if r is None:
                continue
            for tag in ("rho", "tau"):
                r[f"MACRO {tag}"] = round(float(np.mean([r[f"{c} {tag}"] for c in CYPS])), 4)
            table.append(r)

    df = pd.DataFrame(table)
    if df.empty:
        print("нечего сравнивать")
        return
    # Усреднение по сидам: рука без ведущего номера сида.
    df["группа"] = df["рука"].str.replace(r"^\d+\|", "", regex=True)
    g = df.groupby(["файл", "группа"], sort=False)[
        ["MACRO rho", "MACRO tau"] + [f"{c} {t}" for c in CYPS for t in ("rho", "tau")]].mean()
    print(g.round(4).to_string())

    print("""
Как читать. Две колонки на фермент: rho считает ВСЕ пары, tau --- только те, у которых полосы
не пересекаются, то есть только те, за которые метрика вообще платит.

Смотреть надо не на уровни, а на РАЗНОСТЬ между руками в каждой из двух колонок. Если рука
проигрывает по rho и не проигрывает по tau, её проигрыш был на парах, которых в зачёте нет ---
и тогда критерий пункта 80 требует оговорки, а не рука отклонения.

CYP3A4 --- диагностический фермент: у него 51.4 % меток ниже порога прибора, значит доля
невидимых пар наибольшая, и расхождение двух колонок должно быть заметнее всего именно там.""")


if __name__ == "__main__":
    main()
