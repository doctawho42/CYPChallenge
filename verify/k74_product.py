"""Умножать два фактора конъюнкции против обучения на метку напрямую. Другой РАНГ, не порог.

Всё, что делалось с классификационным треком до сих пор, упиралось в одно упорядочение --- то,
которое даёт классификатор, обученный на `is_TDI`. Калибровка (пункт 235) монотонна и ранг не
меняет; оценщики порога (k72) тем более. Оракул порога +0.0472 --- это потолок ТОГО ЖЕ ранга.

Произведение --- другая упорядочивающая функция. Метка есть конъюнкция (пункт 234):

    is_TDI  <=>  (Delta > log10 2)  AND  (pi_TDI > 4.301)

Обучаем классификатор на каждый конъюнкт отдельно и перемножаем вероятности. При корреляции
ошибок плеч 0.84 (пункт 229) произведение --- не совместная вероятность, поэтому это скор для
порога, а не калиброванная величина. Но оно использует структуру, которую прямой классификатор
выбрасывает, и первые два сида в k71 дали ему 0.3629 против 0.3394 у двух жёстких порогов и
~0.315 у подаваемого.

**AUC здесь первоклассная величина, а не справка.** Если AUC произведения выше AUC прямого
классификатора, то заявление не про порог: структура метки даёт лучший ранг, чем обучение на
метку. Если AUC равен, а MCC выше --- значит выигрыш всё-таки пороговый и его съест калибровка.

**Два контроля, без которых число ничего не значит.** На CYP3A4 одни ворота дают ~0.30, а
произведение 0.36, то есть фактор сдвига добавляет ~0.06 --- вчетверо больше, чем те +0.027,
которые сдвиг давал при жёстких порогах (пункт 238). Утверждение сильное, поэтому каждый фактор
перемешивается по строкам по очереди: если произведение с перемешанным сдвигом не падает, значит
оно просто переупакованные ворота.

Плюс `min(gp, dp)` --- другой оператор конъюнкции, совпадающий с совместной вероятностью при
полной зависимости, тогда как произведение отвечает независимости. Истина между ними, и какой
ближе --- вопрос измерения, а не выбора.

Пороги ставятся тремя способами (plug-in, argmax, центроид) на каждом скоре, потому что k72
показал: argmax --- плохой оценщик, и сравнивать скоры надо не через него одного.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES

import argparse
import json
import time

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import matthews_corrcoef, roc_auc_score

from cypsplit import butina_folds
from submit import gbm_clf, plugin_threshold

L = np.log10(2.0)
GATE = 4.0 + L
GRID = np.linspace(0.001, 0.999, 200)


def oof_clf(X, y, fold, tag):
    p = np.zeros(len(y))
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0 or y[trn].sum() == 0:
            continue
        t0 = time.time()
        p[te] = gbm_clf().fit(X[trn], y[trn].astype(int)).predict_proba(X[te])[:, 1]
        print(f"      {tag} фолд {f}: {int(trn.sum())} обуч, {time.time()-t0:.0f} с", flush=True)
    return p


def cut_argmax(p, y, g):
    return float(g[int(np.argmax([matthews_corrcoef(y, p >= t) for t in g]))])


def cut_centroid(p, y, g, n=200, seed=0):
    c = np.array([matthews_corrcoef(y, p >= t) for t in g])
    rng = np.random.default_rng(seed)
    peaks = [np.max([matthews_corrcoef(y[k], p[k] >= t) for t in g])
             for k in (rng.integers(0, len(y), len(y)) for _ in range(n))]
    ok = c >= c.max() - float(np.std(peaks))
    return float(g[ok].mean())


def apply_oof(score, y, fold, how):
    """Порог подбирается на обучающих фолдах, применяется к тестовому."""
    out = np.zeros(len(y), bool)
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0:
            continue
        g = np.quantile(score[trn], np.linspace(0.02, 0.98, 97))
        if how == "plug-in":
            t = plugin_threshold(score[trn] / max(score[trn].max(), 1e-9), GRID) * \
                max(score[trn].max(), 1e-9)
        elif how == "argmax":
            t = cut_argmax(score[trn], y[trn], g)
        else:
            t = cut_centroid(score[trn], y[trn], g)
        out[te] = score[te] >= t
    return out


def platt_oof(p, y, fold):
    def lg(q):
        q = np.clip(q, 1e-6, 1 - 1e-6)
        return np.log(q / (1 - q)).reshape(-1, 1)
    out = np.zeros(len(p))
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0:
            continue
        lr = LogisticRegression(C=1e6).fit(lg(p[trn]), y[trn])
        out[te] = lr.predict_proba(lg(p[te]))[:, 1]
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/product.json")
    a = ap.parse_args()

    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tdi = pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv").set_index("Molecule_Name").reindex(rows.Molecule_Name)
    inh = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv").set_index("Molecule_Name").reindex(rows.Molecule_Name)

    recs = []
    for seed in [int(s) for s in a.seeds.split(",")]:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        print(f"\n=== сид {seed} ===", flush=True)
        for c in ["CYP3A4", "CYP2D6"]:
            lab = tdi[f"{c}_is_TDI"]
            ta = tdi[f"{c}_pIC50_TDI_condition"].to_numpy(float)
            dr = inh[f"{c}_pIC50_direct_inhibition"].to_numpy(float)
            m = (lab.notna() & np.isfinite(ta) & np.isfinite(dr)).to_numpy()
            y = lab[m].astype(bool).to_numpy()
            Xi, fi = X[m], fold[m]
            gate_t, shift_t = ta[m] > GATE, (ta[m] - dr[m]) > L
            print(f"  {c}: n {len(y)}, положительных {y.mean():.3f}", flush=True)

            lp = oof_clf(Xi, y, fi, f"{c} метка")
            gp = oof_clf(Xi, gate_t, fi, f"{c} ворота")
            dp = oof_clf(Xi, shift_t, fi, f"{c} сдвиг")
            lpc = platt_oof(lp, y, fi)

            rng = np.random.default_rng(seed)
            scores = {
                "метка напрямую":            lp,
                "метка + Платт":             lpc,
                "ворота одни":               gp,
                "сдвиг один":                dp,
                "ПРОИЗВЕДЕНИЕ":              gp * dp,
                "минимум":                   np.minimum(gp, dp),
                "произв., сдвиг перемешан":  gp * dp[rng.permutation(len(dp))],
                "произв., ворота перемешаны": gp[rng.permutation(len(gp))] * dp,
            }
            au = {k: roc_auc_score(y, v) for k, v in scores.items()}
            print(f"    AUC: " + "  ".join(f"{k} {v:.4f}" for k, v in au.items()
                                           if k in ("метка напрямую", "ПРОИЗВЕДЕНИЕ",
                                                    "минимум", "ворота одни")), flush=True)
            for k, v in scores.items():
                for how in ("plug-in", "argmax", "центроид"):
                    dec = apply_oof(v, y, fi, how)
                    mcc = matthews_corrcoef(y, dec)
                    recs.append({"seed": seed, "cyp": c, "score": k, "thr": how,
                                 "mcc": float(mcc), "auc": float(au[k]),
                                 "rate": float(dec.mean())})
                best = max(r["mcc"] for r in recs[-3:])
                print(f"      {k:26s} AUC {au[k]:.4f}  MCC "
                      + " ".join(f"{r['thr']} {r['mcc']:+.4f}" for r in recs[-3:])
                      + f"   лучший {best:+.4f}", flush=True)
        json.dump(recs, open(a.out, "w"))

    df = pd.DataFrame(recs)
    print("\n\nMCC, среднее по сидам (строки --- скор, столбцы --- порог)")
    for c in ["CYP3A4", "CYP2D6"]:
        print(f"\n{c}:")
        print(df[df.cyp == c].pivot_table(index="score", columns="thr", values="mcc").round(4).to_string())
        print("AUC:", df[df.cyp == c].groupby("score").auc.mean().round(4).to_dict())
    print(f"\nсохранено: {a.out}", flush=True)


if __name__ == "__main__":
    main()
