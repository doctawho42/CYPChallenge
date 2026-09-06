"""Связка двух подпольных выигрышей, померенная как ОДИН объект. Предрегистрация --- пункт 244.

Вопрос, ради которого файл существует: пол --- свойство ИЗМЕРЕНИЯ, а не компонента, поэтому
несколько малых вмешательств можно вводить вместе и переходить пол связкой, даже если ни одно
из них поодиночке пол не переходит. Мёртвая зона --- существующее доказательство этой формы.

Но связывать можно только эффекты малые и НАСТОЯЩИЕ. У отсутствующего эффекта складывать нечего:
связка платит его дисперсию и не получает средней. Различает их устойчивость знака, а не величина
(пункт 213). По этому правилу из закрытого за двое суток выживают ровно два кандидата --- ансамбль
на плече преинкубации (+0.0104, знак 7/8) и конъюнктивный скор (+0.0162..+0.0187, 6-7/8), --- а
квантовый блок исключён: его ПЕРЕМЕШАННЫЙ контроль набирает больше обоих настоящих плеч.

Связка:

    P(ворота)  <- четырёхчленный ансамбль с мёртвой зоной на pi_TDI, затем одномерный
                  калибратор pi_hat - 4.301 -> вероятность, подогнанный вне фолда
    P(сдвиг)   <- классификатор на (Delta > log10 2)
    скор        = P(ворота) * P(сдвиг),  калибруется, порог plug-in

Меряется против ПОДАВАЕМОГО (структурный классификатор + Платт), четыре сида, оба эндпоинта.
Условия приёмки записаны в пункте 244 до прогона и здесь не повторяются, чтобы их нельзя было
подправить по дороге.

Контрольные плечи, без которых число ничего не значит: каждый фактор поодиночке, произведение с
ПЕРЕМЕШАННЫМ вторым фактором, и голое плечо вместо ансамблевого --- чтобы видеть, что именно из
связки работает.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES

import argparse, json, time
import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import matthews_corrcoef, roc_auc_score

from cypsplit import butina_folds
from submit import gbm_clf, gbm_reg, plugin_threshold
from k75_armens import member_oof, dz_refit, TDI, GATE, L

GRID = np.linspace(0.001, 0.999, 200)


def lg(q):
    q = np.clip(q, 1e-6, 1 - 1e-6)
    return np.log(q / (1 - q)).reshape(-1, 1)


def platt_oof(score, y, fold, raw=False):
    """Одномерный калибратор вне фолда. raw=True: вход не вероятность, а действительное число."""
    out = np.zeros(len(score))
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0 or y[trn].sum() in (0, trn.sum()):
            continue
        A = score.reshape(-1, 1) if raw else lg(score)
        lr = LogisticRegression(C=1e6).fit(A[trn], y[trn])
        out[te] = lr.predict_proba(A[te])[:, 1]
    return out


def oof_clf(X, y, fold, tag):
    p = np.zeros(len(y))
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0 or y[trn].sum() == 0:
            continue
        t0 = time.time()
        p[te] = gbm_clf().fit(X[trn], y[trn].astype(int)).predict_proba(X[te])[:, 1]
        print(f"      {tag} фолд {f}: {time.time()-t0:.0f} с", flush=True)
    return p


def decide(score, y, fold):
    """Калибровка произведения, затем plug-in. Порог берётся на обучающих фолдах."""
    cal = platt_oof(score, y, fold)
    out = np.zeros(len(y), bool)
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0:
            continue
        out[te] = cal[te] >= plugin_threshold(cal[trn], GRID)
    return out, cal


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/bundle.json")
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
        Xs, ys, folds, LOs, HIs, labs, shifts, gates = [], [], [], [], [], [], [], []
        for c in TDI:
            lab = tdi[f"{c}_is_TDI"]
            ta = tdi[f"{c}_pIC50_TDI_condition"].to_numpy(float)
            dr = inh[f"{c}_pIC50_direct_inhibition"].to_numpy(float)
            lo = tdi[f"{c}_pIC50_TDI_condition_conf_low"].to_numpy(float)
            hi = tdi[f"{c}_pIC50_TDI_condition_conf_high"].to_numpy(float)
            m = (lab.notna() & np.isfinite(ta) & np.isfinite(dr)
                 & np.isfinite(lo) & np.isfinite(hi)).to_numpy()
            Xs.append(X[m]); ys.append(ta[m]); folds.append(fold[m])
            LOs.append(lo[m]); HIs.append(hi[m])
            labs.append(lab[m].astype(bool).to_numpy())
            shifts.append((ta - dr)[m]); gates.append((ta > GATE)[m])
            print(f"  {c}: n {int(m.sum())}", flush=True)

        KINDS = ["поферментно", "пул", "GP", "гребневая"]
        dz = {}
        for k in KINDS:
            t0 = time.time()
            P = member_oof(k, Xs, ys, folds)
            tg = [np.clip(P[e], LOs[e], HIs[e]) for e in range(len(TDI))]
            dz[k] = dz_refit(k, Xs, tg, folds)
            print(f"    плечо, член {k}: {time.time()-t0:.0f} с", flush=True)
        arm_ens = [np.mean([dz[k][e] for k in KINDS], 0) for e in range(len(TDI))]
        arm_bare = []
        for e in range(len(TDI)):
            p = np.zeros(len(ys[e]))
            for f in range(5):
                trn, te = folds[e] != f, folds[e] == f
                p[te] = gbm_reg().fit(Xs[e][trn], ys[e][trn]).predict(Xs[e][te])
            arm_bare.append(p)

        for e, c in enumerate(TDI):
            y, fi = labs[e], folds[e]
            print(f"  --- {c}", flush=True)
            lp = oof_clf(Xs[e], y, fi, f"{c} метка")
            dp = oof_clf(Xs[e], shifts[e] > L, fi, f"{c} сдвиг")
            gp_ens = platt_oof(arm_ens[e] - GATE, gates[e], fi, raw=True)
            gp_bare = platt_oof(arm_bare[e] - GATE, gates[e], fi, raw=True)
            rng = np.random.default_rng(seed)
            arms = {
                "подаётся: метка + Платт": lp,
                "СВЯЗКА: ворота(ансамбль) * сдвиг": gp_ens * dp,
                "контроль: ворота(голое плечо) * сдвиг": gp_bare * dp,
                "контроль: ворота(ансамбль) одни": gp_ens,
                "контроль: сдвиг один": dp,
                "контроль: ворота(ансамбль) * ПЕРЕМЕШАННЫЙ сдвиг": gp_ens * dp[rng.permutation(len(dp))],
            }
            for k, s in arms.items():
                dec, cal = decide(s, y, fi)
                mcc = matthews_corrcoef(y, dec)
                auc = roc_auc_score(y, s)
                print(f"      {k:46s} MCC {mcc:+.4f}  AUC {auc:.4f}  доля {dec.mean():.3f}", flush=True)
                recs.append({"seed": seed, "cyp": c, "arm": k, "mcc": float(mcc),
                             "auc": float(auc), "rate": float(dec.mean())})
        json.dump(recs, open(a.out, "w"))

    df = pd.DataFrame(recs)
    print("\n\nMCC, среднее по сидам")
    print(df.pivot_table(index="arm", columns="cyp", values="mcc").round(4).to_string())
    print("\nAUC, среднее по сидам")
    print(df.pivot_table(index="arm", columns="cyp", values="auc").round(4).to_string())
    piv = df.pivot_table(index=["seed", "cyp"], columns="arm", values="mcc")
    base = "подаётся: метка + Платт"
    print("\n=== ПРОВЕРКА ПРЕДРЕГИСТРАЦИИ (пункт 244), макро-пол 0.0076")
    for k in piv.columns:
        if k == base:
            continue
        d = piv[k] - piv[base]
        mac = d.groupby(level="seed").mean()
        ok = (mac.mean() > 0, int((d > 0).sum()) >= 6, mac.mean() > 0.0076)
        print(f"  {k:46s} макро {mac.mean():+.4f}  знак {int((d>0).sum())}/{len(d)}  "
              f"условия {'+'.join('да' if x else 'НЕТ' for x in ok)}")
    print(f"\nсохранено: {a.out}", flush=True)


if __name__ == "__main__":
    main()
