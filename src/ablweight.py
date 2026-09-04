"""Down-weighting the rows whose curve fit fell apart -- and the control that decides whether
that is a new thing or the dead zone wearing a hat.

Section 3 of the document observes that label reliability is 0.90 to 0.95 but the error is
distributed very unevenly, and says outright that weighting compounds by their individual error
makes more sense than it looks. Measured here, the top five per cent by sigma carry 46, 42, 57 and
22 per cent of the total squared error -- real, and far less uniform than the document's "almost
sixty per cent" suggests.

**The objection that has to be answered first, and it is fatal to the naive version.** Item 114
established that sigma is a deterministic function of the label: isotonic regression of the band
half-width on the label value gives R^2 of 0.928 to 0.970, reproduced here. So weighting by
1/sigma^2 is weighting by potency, and the dead zone already gives exactly those rows a wide free
zone. The naive arm cannot be new; it can only be a re-expression, and item 181 is what that looks
like when it is measured properly.

**What survives the objection.** After removing the isotonic fit, 3.0 to 7.2 per cent of sigma's
variance remains, and its correlation with the label is 0.03, -0.04, 0.11 and -0.02 -- the label
dependence is gone. Inside that residual the tail persists: the top five per cent of residuals
carry 50 to 60 per cent of residual variance. **That residual is the honest signal**: the part of a
compound's stated error that its potency does not predict, which is what "the fit fell apart" means
once potency is accounted for.

Arms, all measured WITH the dead zone, because that is the submitted configuration and because
item 181 is the standing lesson that an intervention overlapping it must be measured on top of it
rather than against a bare squared loss:

    контроль                  без взвешивания
    1/sigma^2                 наивное; ожидается перевыражение мёртвой зоны
    1/остаток^2               честное; веса из ОСТАТКА sigma после изотоники
    отбросить верхние 5 %     жёсткая версия честного
    1/остаток^2, перемешан    контроль: тот же маргинал весов, связь со строкой разорвана

Every weight is fitted on training folds only -- the isotonic that defines the residual never sees
the fold it is applied to.

Same folds, learner and metric as src/ablate.py. Writes results/preds/oof_weight.json.
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
from sklearn.isotonic import IsotonicRegression
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
DZ = dict(KW, loss="absolute_error")


def weights(arm, y, sig, trn, rng):
    """Веса строк. Всё, что оценивается, оценивается ТОЛЬКО на обучающих фолдах."""
    w = np.ones(len(y))
    if arm == "контроль":
        return w
    if arm == "1/sigma^2":
        w = 1.0 / np.maximum(sig, 1e-3) ** 2
    else:
        iso = IsotonicRegression(increasing="auto", out_of_bounds="clip")
        iso.fit(y[trn], sig[trn])
        r = np.abs(sig - iso.predict(y))
        if "перемешан" in arm:
            r = r[rng.permutation(len(r))]
        if "отбросить" in arm:
            thr = np.quantile(r[trn], 0.95)
            w = (r <= thr).astype(float)
        else:
            w = 1.0 / np.maximum(r, np.quantile(np.maximum(r[trn], 1e-6), 0.10)) ** 2
    w = w / w[trn].mean()          # нормировка, чтобы сила регуляризации не поехала
    return np.clip(w, 0.0, 50.0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--arms", default="контроль|1/sigma^2|1/остаток^2|"
                                      "отбросить верхние 5%|1/остаток^2, перемешан")
    ap.add_argument("--out", default=RES + "preds/oof_weight.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    arms = a.arms.split("|")

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)

    out, table, TGT = {}, [], {}
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for arm in arms:
            t0 = time.time()
            rng = np.random.default_rng(3000 + seed)
            r = {"seed": seed, "рука": arm}
            for e, c in enumerate(CYPS):
                te0 = time.time()
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                sig = (hi - lo) / 3.92
                Xi, fi = X[m], fold[m]

                # Мишень мёртвой зоны считается ОДИН раз, БЕЗ весов, и переиспользуется
                # всеми руками. Иначе руки различались бы не только весами, но и мишенью,
                # и сравнение мерило быдве вещи сразу. Кэш по (сид, фермент).
                ck = (seed, c)
                if ck not in TGT:
                    p0 = np.zeros(len(y))
                    for f in range(5):
                        a_, b_ = fi != f, fi == f
                        if b_.sum() == 0:
                            continue
                        p0[b_] = HistGradientBoostingRegressor(**KW).fit(
                            Xi[a_], y[a_]).predict(Xi[b_])
                    TGT[ck] = np.clip(p0, lo, hi)
                    print(f"      мишень {c} посчитана ({time.time()-te0:.0f} с)", flush=True)
                tgt = TGT[ck]
                p = np.zeros(len(y))
                for f in range(5):
                    a_, b_ = fi != f, fi == f
                    if b_.sum() == 0:
                        continue
                    w = weights(arm, y, sig, a_, rng)
                    p[b_] = HistGradientBoostingRegressor(**DZ).fit(
                        Xi[a_], tgt[a_], sample_weight=w[a_]).predict(Xi[b_])

                print(f"      {arm:24s} {c} ({time.time()-te0:.0f} с)", flush=True)
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{arm}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for t in ("пара", "rho"):
                r[f"MACRO {t}"] = round(float(np.mean([r[f"{c} {t}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {arm:26s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                         + [f"{c} rho" for c in CYPS]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Всё сравнивается с рукой «контроль», и все руки идут С мёртвой зоной.

«1/sigma^2» против контроля --- перевыражение или нет. Ожидание: около нуля, потому что sigma
есть функция метки при R^2 0.93--0.97, а мёртвая зона ту же информацию уже забрала жёстким
допуском. Заметный выигрыш здесь означал бы, что градуированный вес добавляет к жёсткому
допуску, чего пункт 200 не нашёл в другой форме.

«1/остаток^2» --- честная версия: вес по той части заявленной погрешности, которую потенция НЕ
предсказывает. Это единственная рука, у которой есть право быть новой.

«перемешан» --- решает. Тот же маргинал весов, разорвана связь веса со строкой. Если честная
рука идёт наравне с перемешанной, работает распределение весов, а не то, какой строке какой вес
достался, и направление закрыто.""")


if __name__ == "__main__":
    main()
