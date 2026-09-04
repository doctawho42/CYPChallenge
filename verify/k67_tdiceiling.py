"""Where is the TDI ceiling? An oracle decomposition of a label defined by a difference.

MCC of 0.32 on CYP3A4 and 0.12 on CYP2D6 is bad, and classification is a third of the leaderboard.
Before trying anything, this file asks where the ceiling is, because the label's definition makes
one answer arithmetically likely and it should be checked rather than assumed.

**The rule, verified at 100.0 per cent against the measured arms (item 197):**

    is_TDI  =  (pi_tdi - pi_dir) > log10(2) = 0.301,   если pi_dir > 4
               pi_tdi > 4.301,                          иначе

The first branch is a **difference of two quantities against a threshold of 0.301**, and our
root-mean-square error on each arm is near 0.7. If the two errors were independent the difference
would carry an error near 1.0 against a signal of 0.301 -- three times more noise than signal, and
an MCC of 0.3 would be roughly what the arithmetic permits rather than a failure of effort.

The premise is checkable and is not obviously true: the two arms are measured on the same plate from
the same compound, so their errors may be strongly correlated and cancel in the difference. Section 7
already measured a consequence of this -- predicting the shift with its own head beats differencing
two regressions, 0.323 against 0.255 of MCC on CYP3A4 -- without measuring the correlation itself.

**Five oracles, each substituting truth for one part**, which is item 128's procedure carried to the
classification track for the first time:

    правило на ИСТИННЫХ плечах        1.000 по построению --- проверка стенда
    истинное прямое + предсказанная Delta   сколько стоит ошибка ТОЛЬКО в Delta
    предсказанное прямое + истинная Delta   сколько стоит ошибка ТОЛЬКО в прямом плече
    оба предсказаны                          что даёт маршрут через правило
    структурный классификатор                что подаётся сейчас

Reading them together says which half to work on, and whether either is worth working on at all.

**Delta is predicted by its own model**, not by differencing, for the reason section 7 gives. Its
error is reported beside the arms' so that the correlation that makes differencing worse is visible
as a number rather than inferred.

Reads data/feats.npz, data/rows.csv and both label files. CYP2D6 and CYP3A4 only. Same folds and
learner as src/ablate.py.
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
from sklearn.ensemble import HistGradientBoostingClassifier, HistGradientBoostingRegressor
from sklearn.metrics import matthews_corrcoef

from cypsplit import butina_folds

TDI_CYPS = ["CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)


def rule(direct, tdiarm):
    return np.where(direct > 4.0, (tdiarm - direct) > np.log10(2.0), tdiarm > 4.301)


def oof(X, y, fold, clf=False):
    p = np.full(len(y), np.nan)
    obs = np.isfinite(y) if not clf else np.isfinite(y.astype(float))
    for f in range(5):
        trn, te = obs & (fold != f), fold == f
        if te.sum() == 0 or trn.sum() < 50:
            continue
        if clf:
            m = HistGradientBoostingClassifier(**KW).fit(X[trn], y[trn].astype(int))
            p[te] = m.predict_proba(X[te])[:, 1]
        else:
            p[te] = HistGradientBoostingRegressor(**KW).fit(X[trn], y[trn]).predict(X[te])
    return p


def best_mcc(score, truth):
    """Наилучший достижимый MCC по порогу --- потолок маршрута, а не его подача."""
    qs = np.quantile(score[np.isfinite(score)], np.linspace(0.02, 0.98, 97))
    return max(matthews_corrcoef(truth, score >= t) for t in qs)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--out", default=RES + "preds/oof_tdiceil.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    inh = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
             .set_index("Molecule_Name").reindex(rows.Molecule_Name).reset_index())
    tdi = (pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv")
             .set_index("Molecule_Name").reindex(rows.Molecule_Name).reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)

    table = []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for c in TDI_CYPS:
            t0 = time.time()
            d = inh[f"{c}_pIC50_direct_inhibition"].to_numpy(float)
            t = tdi[f"{c}_pIC50_TDI_condition"].to_numpy(float)
            lab = tdi[f"{c}_is_TDI"]
            m = lab.notna().to_numpy() & np.isfinite(d) & np.isfinite(t)
            yb = lab[m].astype(bool).to_numpy()
            dlt = t - d

            pd_ = oof(X, np.where(np.isfinite(d), d, np.nan), fold)
            pt_ = oof(X, np.where(np.isfinite(t), t, np.nan), fold)
            pdl = oof(X, np.where(np.isfinite(dlt), dlt, np.nan), fold)
            pcl = oof(X, lab.to_numpy(), fold, clf=True)

            ed, et, el = d[m] - pd_[m], t[m] - pt_[m], dlt[m] - pdl[m]
            r = {"seed": seed, "фермент": c, "n": int(m.sum()),
                 "положительных": int(yb.sum()),
                 "ско прям": round(float(np.sqrt((ed ** 2).mean())), 3),
                 "ско tdi": round(float(np.sqrt((et ** 2).mean())), 3),
                 "ско Delta прямо": round(float(np.sqrt((el ** 2).mean())), 3),
                 "ско разности плеч": round(float(np.sqrt(((et - ed) ** 2).mean())), 3),
                 "corr ошибок плеч": round(float(np.corrcoef(ed, et)[0, 1]), 3)}

            r["оракул: истинные плечи"] = round(
                float(matthews_corrcoef(yb, rule(d[m], t[m]))), 3)
            r["истин.прям + предск.Delta"] = round(
                float(matthews_corrcoef(yb, rule(d[m], d[m] + pdl[m]))), 3)
            r["предск.прям + истин.Delta"] = round(
                float(matthews_corrcoef(yb, rule(pd_[m], pd_[m] + dlt[m]))), 3)
            r["оба предсказаны"] = round(
                float(matthews_corrcoef(yb, rule(pd_[m], pd_[m] + pdl[m]))), 3)
            r["классификатор, потолок"] = round(float(best_mcc(pcl[m], yb)), 3)
            r["Delta как счёт, потолок"] = round(float(best_mcc(pdl[m], yb)), 3)
            table.append(r)
            print(f"  сид {seed} {c} посчитан ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print("ОШИБКИ ПРЕДСКАЗАНИЯ, ско в единицах pIC50")
    print(df.set_index(["seed", "фермент"])[
        ["n", "положительных", "ско прям", "ско tdi", "ско Delta прямо",
         "ско разности плеч", "corr ошибок плеч"]].to_string())
    print()
    print("ОРАКУЛЫ, MCC")
    print(df.set_index(["seed", "фермент"])[
        ["оракул: истинные плечи", "истин.прям + предск.Delta",
         "предск.прям + истин.Delta", "оба предсказаны",
         "классификатор, потолок", "Delta как счёт, потолок"]].to_string())
    json.dump({"table": table}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Первая строка оракулов обязана быть единицей: правило воспроизводит метку из
истинных плеч по построению (пункт 197). Если нет --- стенд сломан и дальше читать нечего.

Дальше две строки отвечают на вопрос, ради которого файл написан. «истинное прямое +
предсказанная Delta» показывает, сколько стоит ошибка ТОЛЬКО в сдвиге; «предсказанное прямое +
истинная Delta» --- сколько стоит ошибка только в уровне. Та, что ниже, и есть узкое место, и
работать надо там.

«ско разности плеч» против «ско Delta прямо» --- цена дифференцирования. Если разность плеч
заметно хуже прямого предсказания сдвига, ошибки плеч коррелируют слабее, чем нужно, и это
объясняет, почему маршрут через два плеча проигрывает.

«классификатор, потолок» --- наилучший MCC по любому порогу, а не по подстановочному, так что
он сравним с оракулами и отделяет вопрос порога от вопроса упорядочения.""")


if __name__ == "__main__":
    main()
