"""Dry run for item 203: what does a manufactured label cost where manufacture is best conditioned?

Item 203 found 1238 CYP3A4 measurements in `data/cyp-challenge-TRAIN_TDI.csv` that no script in
this repository reads -- molecules carrying a pre-incubation arm and no direct label, dropped
everywhere because every consumer realigns against `data/rows.csv`. Using them means MANUFACTURING
a direct label from the TDI arm and the fitted shift, and item 203 pre-registered the test that
decides whether that is worth doing, before any feature is built for them.

**The test.** Take the 2334 CYP3A4 molecules that carry BOTH arms -- where the manufacture is best
conditioned that it will ever be, since the shift is estimated on exactly this population --
**throw their real labels away**, replace them with manufactured ones, and measure what that costs.

    если подмена ЗДЕСЬ стоит больше пола CYP3A4 (0.0033), то 1238 таких же меток на
    молекулах, где сдвиг экстраполируется на 0.79 логединицы вверх и доля TDI-положительных
    невычислима, не дадут ничего, и источник закрыт одним числом

This is deliberately harsher than the real use, and the asymmetry is the point: replacing measures
label quality in isolation, while the real use ADDS rows and would confound quality with sample
size. A third arm replaces only a random 35 per cent, the ratio 1238 would make against 2335, as
the nearest honest simulation of the real thing.

**Why the shift is regressed rather than subtracted.** Item 203 measured the shift at +0.3388 with
sd 0.3638 and se 0.0075 -- the best-conditioned offset in this file, 111 times NCGC's molecules at
28 times its precision (item 144). But it is not constant: regressed on the TDI arm the slope is
+0.0852 at correlation 0.265, and the 1238 sit 0.79 log units more potent than the paired ones,
where a constant would carry a systematic +0.067. Both forms are run, because if the constant does
as well as the regression then the slope is noise and one fewer thing needs estimating.

**Everything is fitted on training folds only.** The shift regression sees the held-out fold's
labels through neither side. Scoring is always against the REAL labels and the REAL bands -- only
the training target is manufactured, which is the whole construction.

**The permutation arm is what separates two different failures.** If manufacturing from a SHUFFLED
TDI arm scores like manufacturing from the real one, then the pre-incubation reading carries
nothing molecule-specific and the source is dead for a reason that has nothing to do with the
shift. If the real arm beats the shuffled one and still loses to the true labels, the reading is
informative and merely too noisy -- which is what item 210 found for the barrier, and would be the
second instance of the same shape.

CYP3A4 only: it is the only enzyme with a TDI arm worth speaking of, and item 203's rows are all
its. Same folds, learner and metric as src/ablate.py. Writes results/preds/oof_dry.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

C = "CYP3A4"
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--frac", type=float, default=1238 / (1238 + 2335),
                    help="доля подменяемых в смешанной руке; умолчание --- доля 1238 среди 3573")
    ap.add_argument("--out", default=RES + "preds/oof_dry.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    tdi = pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv").set_index("Molecule_Name")
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)

    col = f"{C}_pIC50_direct_inhibition"
    m = tr[col].notna().to_numpy()
    y = tr[col].to_numpy(float)
    lo = tr[col + "_conf_low"].to_numpy(float)
    hi = tr[col + "_conf_high"].to_numpy(float)
    arm = tdi.reindex(rows.Molecule_Name)[f"{C}_pIC50_TDI_condition"].to_numpy(float)
    both = m & ~np.isnan(arm)
    print(f"{C}: прямых {int(m.sum())}, плечо TDI {int((~np.isnan(arm)).sum())}, "
          f"ОБА {int(both.sum())}", flush=True)
    print(f"подменяются только строки с ОБОИМИ плечами; строк без плеча "
          f"{int((m & np.isnan(arm)).sum())} и они всегда настоящие\n", flush=True)

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        rng = np.random.default_rng(7000 + seed)
        half = rng.random(len(rows)) < a.frac
        shuf = arm.copy()
        idx = np.where(~np.isnan(arm))[0]
        shuf[idx] = arm[idx][rng.permutation(len(idx))]

        for name in ("настоящие", "изготовленные, регрессия", "изготовленные, константа",
                     f"изготовленные у {a.frac:.0%}", "изготовленные, плечо перемешано"):
            t0 = time.time()
            src = shuf if "перемешано" in name else arm
            sel = both if "%" not in name else (both & half)
            p = np.zeros(int(m.sum()))
            for f in range(5):
                trn_mol, te_mol = fold != f, fold == f
                te = (m & te_mol)[m]
                if te.sum() == 0:
                    continue
                # Мишень: настоящая метка везде, кроме подменяемых строк.
                t = y.copy()
                if name != "настоящие":
                    fit = both & trn_mol                      # оценка сдвига --- ТОЛЬКО обучение
                    d = src[fit] - y[fit]
                    if "константа" in name:
                        est = np.full(len(src), float(d.mean()))
                    else:
                        b, c0 = np.polyfit(src[fit], d, 1)
                        est = b * src + c0
                    rep = sel & trn_mol
                    t[rep] = src[rep] - est[rep]
                trn = m & trn_mol
                mdl = HistGradientBoostingRegressor(**KW).fit(X[trn], t[trn])
                p[te] = mdl.predict(X[(m & te_mol)])
            yy, ll, hh, fi = y[m], lo[m], hi[m], fold[m]
            q = fit_apply(p, ll, hh, fi, np.ones(len(yy)) / len(yy))
            r = {"seed": seed, "рука": name,
                 "подменено": int((sel & (fold >= 0)).sum()),
                 "пара": round(float(strae(yy, q, y_true_upper=hh, y_true_lower=ll)), 4),
                 "rho": round(float(spearmanr(yy, p).statistic), 4)}
            out[f"{seed}|{name}"] = p.tolist()
            table.append(r)
            print(f"  сид {seed} {name:34s} пара {r['пара']:.4f} rho {r['rho']:.4f} "
                  f"(подменено {r['подменено']}, {time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    g = df.groupby("рука", sort=False)[["пара", "rho"]].mean()
    base = g.loc["настоящие"]
    g["ранг минус настоящие"] = g["rho"] - base["rho"]
    g["пара минус настоящие"] = g["пара"] - base["пара"]
    print(g.round(4).to_string())
    print()
    for nm in g.index:
        if nm == "настоящие":
            continue
        d = np.array([df[(df.seed == s) & (df["рука"] == nm)].rho.iloc[0]
                      - df[(df.seed == s) & (df["рука"] == "настоящие")].rho.iloc[0]
                      for s in seeds])
        print(f"{nm:34s} ранг {d.mean():+.4f} знак {int((d < 0).sum())}/4 против  "
              + " ".join(f"{v:+.4f}" for v in d))
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print(f"""
Как читать. Пол CYP3A4 --- 0.0033, самый низкий из четырёх.

Предрегистрация пункта 203: если подмена ЗДЕСЬ, на популяции, где сдвиг оценивается лучше
всего, что когда-либо будет, стоит больше пола, то 1238 таких же меток на молекулах, где
сдвиг экстраполируется на 0.79 логединицы вверх и доля TDI-положительных невычислима, не
дадут ничего. Источник закрывается одним числом.

Рука с перемешанным плечом разделяет два разных провала. Если она идёт наравне с настоящей,
преинкубационное показание не несёт ничего посоединениевого, и дело не в сдвиге вовсе. Если
настоящая её обыгрывает и всё равно проигрывает настоящим меткам, показание информативно и
лишь слишком шумно --- та же форма, что пункт 210 нашёл для барьера.""")


if __name__ == "__main__":
    main()
