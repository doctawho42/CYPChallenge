"""Свип max_features у HistGB на двух бустинговых членах. Предрегистрация --- пункт 270.

Пункт 140 померил субсэмплинг столбцов как +0.0123 ранга, монотонно по ручке, знак 4/4 ---
и закрыл вопрос строкой «у пинованного HistGradientBoostingRegressor нет max_features вовсе,
проверено». Проверено было против УСТАНОВЛЕННОЙ 1.3.2, тогда как пин --- >=1.3,<1.9, а в 1.8.0
параметр есть (пункт 266). За 269 пунктов ни один гиперпараметр HistGB не свипался ни разу.

Запускать ТОЛЬКО под оверлеем:  uv run --with 'scikit-learn==1.8.0' python verify/k82_maxfeat.py

Два контроля, без которых числа не читаются:
  * max_features=1.0 бит в бит совпадает с отсутствием параметра (проверено, max|d| = 0.00e+00),
    поэтому умолчание --- no-op и ничто опубликованное не двигается;
  * поферментный член при 1.0 ОБЯЗАН воспроизвести results/preds/oof.json: макро-4 0.565046.
    Скрипт это утверждает и падает, если нет.

Условия приёмки записаны в пункте 270 до прогона и здесь не повторяются. Пишет
results/logs/k82_maxfeat.json после каждого сида.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse, json, time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
# Пины gbm_reg() из src/submit.py, скопированные, а не импортированные: max_features
# добавляется здесь, и подмешивать его в общую фабрику до решения нельзя.
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
GRID = [1.0, 0.3, 0.1, 0.03]
ANCHOR = 0.565046          # макро-4 поферментного члена при mf=1.0, из oof.json


def mdl(mf):
    return HistGradientBoostingRegressor(max_features=mf, **KW)


def pooled_design(Xr, e):
    ind = np.zeros((len(Xr), len(CYPS)), np.float32)
    ind[:, e] = 1.0
    return np.hstack([Xr, ind])


def perenz(X, Y, mask, fold, mf):
    P = [np.zeros(int(mask[:, e].sum())) for e in range(len(CYPS))]
    for e in range(len(CYPS)):
        m = mask[:, e]
        Xi, yy, fi = X[m], Y[m, e], fold[m]
        for f in range(5):
            a, b = fi != f, fi == f
            if b.sum() == 0:
                continue
            P[e][b] = mdl(mf).fit(Xi[a], yy[a]).predict(Xi[b])
    return P


def pooled(X, Y, mask, fold, mf):
    P = [np.zeros(int(mask[:, e].sum())) for e in range(len(CYPS))]
    for f in range(5):
        Xs = [pooled_design(X[mask[:, e] & (fold != f)], e) for e in range(len(CYPS))]
        ys = [Y[mask[:, e] & (fold != f), e] for e in range(len(CYPS))]
        m_ = mdl(mf).fit(np.vstack(Xs), np.concatenate(ys))
        for e in range(len(CYPS)):
            mk = mask[:, e]
            b = fold[mk] == f
            if b.sum():
                P[e][b] = m_.predict(pooled_design(X[mk][b], e))
    return P


def score(P, Y, LO, HI, mask, fold):
    rk, pr = [], []
    for e in range(len(CYPS)):
        m = mask[:, e]
        yy, lo, hi = Y[m, e], LO[m, e], HI[m, e]
        u = np.ones(len(yy)) / len(yy)
        rk.append(float(spearmanr(yy, P[e]).statistic))
        pr.append(float(strae(yy, fit_apply(P[e], lo, hi, fold[m], u),
                              y_true_upper=hi, y_true_lower=lo)))
    return float(np.mean(rk)), float(np.mean(pr)), rk, pr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--members", default="поферментно,пул")
    a = ap.parse_args()
    import sklearn
    if tuple(int(x) for x in sklearn.__version__.split(".")[:2]) < (1, 4):
        raise SystemExit(f"scikit-learn {sklearn.__version__} без max_features. Запускайте под "
                         f"оверлеем: uv run --with 'scikit-learn==1.8.0' python {__file__}")
    print(f"scikit-learn {sklearn.__version__}\n", flush=True)

    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(Y)
    members = {"поферментно": perenz, "пул": pooled}
    want = [m.strip() for m in a.members.split(",") if m.strip()]
    dst = RES + "logs/k82_maxfeat.json"
    out = json.load(open(dst)) if _pl.Path(dst).exists() else {}

    print(f"  {'сид':>3s} {'член':>12s} {'mf':>5s} {'макро-ранг':>11s} {'макро-пара':>11s} {'c':>6s}",
          flush=True)
    for s in [int(x) for x in a.seeds.split(",") if x.strip()]:
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        out.setdefault(str(s), {})
        for nm in want:
            out[str(s)].setdefault(nm, {})
            for mf in GRID:
                t0 = time.time()
                sc = score(members[nm](X, Y, mask, fold, mf), Y, LO, HI, mask, fold)
                out[str(s)][nm][str(mf)] = sc
                print(f"  {s:3d} {nm:>12s} {mf:5.2f} {sc[0]:11.4f} {sc[1]:11.4f} "
                      f"{time.time() - t0:6.0f}", flush=True)
                if s == 0 and nm == "поферментно" and mf == 1.0 and abs(sc[0] - ANCHOR) > 1e-5:
                    raise SystemExit(f"ЯКОРЬ НЕ СОШЁЛСЯ: {sc[0]:.6f} против {ANCHOR}. "
                                     f"max_features=1.0 обязан быть no-op.")
            base = out[str(s)][nm]["1.0"][0]
            print("      " + "   ".join(f"mf{mf}: {out[str(s)][nm][str(mf)][0]-base:+.4f}"
                                        for mf in GRID[1:]), flush=True)
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(dst, "w"), ensure_ascii=False, indent=1)

    print(f"\n{'член':>12s} {'mf':>5s} {'Δранг':>9s} {'sd':>8s} {'знак':>6s} {'Δпара':>9s}  пол 0.0052",
          flush=True)
    for nm in want:
        ks = [k for k in out if nm in out[k] and all(str(mf) in out[k][nm] for mf in GRID)]
        for mf in GRID[1:]:
            d = np.array([out[k][nm][str(mf)][0] - out[k][nm]["1.0"][0] for k in ks])
            q = np.array([out[k][nm][str(mf)][1] - out[k][nm]["1.0"][1] for k in ks])
            print(f"{nm:>12s} {mf:5.2f} {d.mean():+9.4f} {d.std(ddof=1) if len(d)>1 else 0:8.4f} "
                  f"{int((d>0).sum()):3d}/{len(d):<2d} {q.mean():+9.4f}", flush=True)


if __name__ == "__main__":
    main()
