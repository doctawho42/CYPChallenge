"""The TDI arms of item 337 scored on the SHIPPED path, which is the quantity the board pays for.

WHY THIS FILE EXISTS. `verify/k113_tdiemb.py` measures each arm's AUC and its best MCC over a
threshold grid. The AUC is honest and decides the question it was built for -- a frozen encoder
orders CYP2D6 TDI better than our own features, +0.0403 of AUC with sign 4/4 against a
discriminating control that fails. The MCC it reports is an ORACLE: `best_mcc` sweeps the grid and
takes the maximum against the true labels. That is the right summary for "did the ordering improve"
and the wrong one for choosing which arm to submit, because the board scores MCC at OUR threshold,
and our threshold comes from the shipped rule: Platt calibration out of fold, then `plugin_threshold`
on the calibrated probabilities (items 233, 235). Comparing an oracle gain against item 235's floors,
which were measured on the shipped path, substitutes one quantity for another.

So this file runs the shipped rule on each arm's cached out-of-fold probabilities and reports the MCC
that results.

THE ONE PLACE THIS IS NOT A COPY OF `src/submit.py`. `tdi_calibrate` is hard-wired to arm A: it
recomputes the out-of-fold probabilities itself with `gbm_clf()` on `Xtr`, so it cannot be handed
another arm's vector. The Platt step is reimplemented here for an arbitrary probability vector, and
in one respect more conservatively than the shipped version: shipped fits the Platt map on ALL
training rows and applies it to the blind set, while here the map is fitted OUT OF FOLD (on the other
four folds) and applied to the held-out fold, because there is no blind set to apply it to and
fitting on the rows one then scores would flatter every arm. `plugin_threshold` is used unchanged ---
it is labels-free, reading only the probabilities, so computing it on all rows is what ships and
leaks nothing.

CONTROLS, each of which must succeed:

  * arm A's AUC recomputed from the cached probabilities must reproduce the committed
    `results/logs/k113_tdiemb.json` exactly. If the cache does not hold what the log summarises,
    nothing below describes the same arms.
  * every Platt slope must be positive. Item 235 measured all 40 slopes of the four-seed sweep
    positive with a minimum of +0.084, and `submit.tdi_calibrate` raises on a non-positive slope
    because at CYP2D6's AUC of 0.59 an order reversal is a live failure, not a theoretical one.
  * the calibration must MOVE the predicted positive rate toward the truth. Item 235's reason for
    existing is that the uncalibrated mean predicted rate was 0.099 against a true 0.217 on CYP2D6;
    if the calibrated rate does not close that gap, the step is not doing what it is for.

    uv run python verify/k116_tdiship.py

Reads results/preds/k113_probs_seed{0..3}.json; writes results/logs/k116_tdiship.json.
Exit: 0 ran, 2 a control failed, 3 the probability cache is missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES

import json

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import matthews_corrcoef, roc_auc_score

from cypsplit import butina_folds, fold_digest
from submit import plugin_threshold

TDI = ["CYP3A4", "CYP2D6"]
SEEDS = [0, 1, 2, 3]
ARMS = ["A наши признаки", "B замороженный энкодер", "C среднее A и B", "D контроль DESC+MECH"]
# Item 235, measured on the shipped path over four seeds. These are the floors an arm must clear.
FLOOR = {"CYP3A4": 0.0281, "CYP2D6": 0.0419, "макро": 0.0076}
LOG = RES + "logs/k116_tdiship.json"


def _logit(q):
    q = np.clip(q, 1e-6, 1 - 1e-6)
    return np.log(q / (1 - q)).reshape(-1, 1)


def platt_oof(p, y, fold):
    """Platt map fitted out of fold and applied to the held-out fold.

    `submit.tdi_calibrate` fits on all training rows and applies to the blind set. There is no
    blind set here, so the map is nested instead; fitting it on the rows it then scores would
    flatter every arm by the same unknown amount, which is worse than being conservative.
    """
    out = np.zeros(len(y))
    slopes = []
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0 or y[trn].sum() == 0 or y[trn].sum() == trn.sum():
            out[te] = p[te]
            continue
        lr = LogisticRegression(C=1e6).fit(_logit(p[trn]), y[trn])
        slopes.append(float(lr.coef_[0, 0]))
        out[te] = lr.predict_proba(_logit(p[te]))[:, 1]
    return out, slopes


def main():
    missing = [s for s in SEEDS if not _pl.Path(RES + f"preds/k113_probs_seed{s}.json").exists()]
    if missing:
        print("нет кэша вероятностей для сидов", missing,
              "--- сначала uv run python verify/k113_tdiemb.py --seeds 0,1,2,3")
        return 3

    rows = pd.read_csv(D + "rows.csv")
    tdi = (pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv")
             .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    ref = json.load(open(RES + "logs/k113_tdiemb.json"))

    out = {"meta": {"floors": FLOOR, "seeds": SEEDS}, "cells": {}}
    all_slopes = []
    res = {c: {a: {"oracle": [], "ship": [], "auc": []} for a in ARMS} for c in TDI}

    for s in SEEDS:
        store = json.load(open(RES + f"preds/k113_probs_seed{s}.json"))
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        if store.get("digest") != dg:
            print(f"ОТКАЗ: кэш сида {s} собран на разбиении {store.get('digest')}, "
                  f"а не {dg}")
            return 2
        print(f"\n  сид {s}: дайджест {dg}", flush=True)

        for c in TDI:
            yv = tdi[f"{c}_is_TDI"].to_numpy(float)
            m = ~np.isnan(yv)
            y = yv[m].astype(int)
            fm = fold[m]
            cell = store[c]
            pa, pb, pdd = (np.asarray(cell["A"], float), np.asarray(cell["B"], float),
                           np.asarray(cell["D"], float))
            P = {"A наши признаки": pa, "B замороженный энкодер": pb,
                 "C среднее A и B": (pa + pb) / 2, "D контроль DESC+MECH": pdd}

            # CONTROL: the cache must hold what the committed log summarises.
            i = ref[c]["seeds"].index(s)
            want = ref[c]["арм"]["A наши признаки"]["пары"][i][0]
            got = float(roc_auc_score(y, pa))
            if abs(got - want) > 1e-9:
                print(f"ОТКАЗ: {c} сид {s}: AUC арма A из кэша {got:.6f} против {want:.6f} "
                      f"в логе k113")
                return 2

            print(f"    {c}: истинная доля положительных {y.mean():.3f}", flush=True)
            for a in ARMS:
                p = P[a]
                cal, slopes = platt_oof(p, y, fm)
                all_slopes += slopes
                thr = plugin_threshold(cal)
                mcc_ship = float(matthews_corrcoef(y, cal >= thr))
                grid = np.arange(0.05, 0.96, 0.01)
                mcc_oracle = max(float(matthews_corrcoef(y, p >= t)) for t in grid)
                auc = float(roc_auc_score(y, p))
                res[c][a]["ship"].append(mcc_ship)
                res[c][a]["oracle"].append(mcc_oracle)
                res[c][a]["auc"].append(auc)
                out["cells"].setdefault(f"{c}|{a}|{s}", {}).update(
                    auc=auc, mcc_ship=mcc_ship, mcc_oracle=mcc_oracle, thr=thr,
                    rate_raw=float((p >= thr).mean()), rate_cal=float((cal >= thr).mean()),
                    mean_p_raw=float(p.mean()), mean_p_cal=float(cal.mean()),
                    true_rate=float(y.mean()), slope_min=min(slopes) if slopes else None)
                print(f"      %-26s AUC %.4f  MCC отгруж. %.4f  оракул %.4f  порог %.2f  "
                      "сред.p %.3f -> %.3f" % (a, auc, mcc_ship, mcc_oracle, thr,
                                               p.mean(), cal.mean()), flush=True)

    # ---- CONTROL: every Platt slope positive ---------------------------------------------
    print(f"\n  КОНТРОЛЬ наклонов Платта: {len(all_slopes)} штук, минимум "
          f"{min(all_slopes):+.4f} (пункт 235: все 40 положительны, минимум +0.084)", flush=True)
    if min(all_slopes) <= 0:
        print("ОТКАЗ: неположительный наклон --- карта переворачивает порядок, а не калибрует")
        return 2

    # ---- CONTROL: calibration moves the rate toward the truth ----------------------------
    print("\n  КОНТРОЛЬ калибровки: средняя предсказанная вероятность против истинной доли",
          flush=True)
    worse = 0
    for c in TDI:
        for a in ARMS:
            k = [v for key, v in out["cells"].items() if key.startswith(f"{c}|{a}|")]
            t = np.mean([v["true_rate"] for v in k])
            r0 = np.mean([v["mean_p_raw"] for v in k])
            r1 = np.mean([v["mean_p_cal"] for v in k])
            ok = abs(r1 - t) <= abs(r0 - t)
            worse += 0 if ok else 1
            print("      %-26s %-8s истина %.3f   сырое %.3f -> калибр. %.3f   %s"
                  % (a, c, t, r0, r1, "ближе" if ok else "ДАЛЬШЕ"))
    if worse:
        print(f"ОТКАЗ: калибровка в {worse} клетках уводит долю ОТ истины --- "
              "шаг делает не то, для чего он есть")
        return 2

    # ---- the table -----------------------------------------------------------------------
    print("\n" + "=" * 96)
    print("MCC НА ОТГРУЖЕННОМ ПУТИ (Платт вне фолда + plug-in), четыре сида")
    print("=" * 96)
    base = "A наши признаки"
    print("  %-26s %8s %8s %9s %8s %9s %8s"
          % ("арм", "AUC", "MCC отгр", "Δ к A", "знак", "оракул", "Δ оракул"))
    for c in TDI:
        print(f"  --- {c} (пол пункта 235: {FLOOR[c]:.4f}) ---")
        b = np.array(res[c][base]["ship"])
        for a in ARMS:
            v = np.array(res[c][a]["ship"])
            o = np.array(res[c][a]["oracle"])
            ob = np.array(res[c][base]["oracle"])
            print("  %-26s %8.4f %8.4f %+9.4f %6d/%d %9.4f %+8.4f"
                  % (a, np.mean(res[c][a]["auc"]), v.mean(), v.mean() - b.mean(),
                     int((v > b).sum()), len(v), o.mean(), o.mean() - ob.mean()))
    print("  --- макро (пол %.4f) ---" % FLOOR["макро"])
    for a in ARMS:
        mv = np.mean([np.mean(res[c][a]["ship"]) for c in TDI])
        mb = np.mean([np.mean(res[c][base]["ship"]) for c in TDI])
        mo = np.mean([np.mean(res[c][a]["oracle"]) for c in TDI])
        mob = np.mean([np.mean(res[c][base]["oracle"]) for c in TDI])
        per_seed = [np.mean([res[c][a]["ship"][i] for c in TDI]) for i in range(len(SEEDS))]
        per_base = [np.mean([res[c][base]["ship"][i] for c in TDI]) for i in range(len(SEEDS))]
        sign = sum(1 for x, z in zip(per_seed, per_base) if x > z)
        print("  %-26s %8s %8.4f %+9.4f %6d/%d %9.4f %+8.4f"
              % (a, "", mv, mv - mb, sign, len(SEEDS), mo, mo - mob))
    out["summary"] = {c: {a: {k: list(map(float, v)) for k, v in res[c][a].items()}
                          for a in ARMS} for c in TDI}
    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
    print(f"\n  записано: {LOG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
