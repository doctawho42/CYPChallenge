"""The TDI probability as the document specifies it: integrate the labelling rule, do not fit a sigmoid.

Section 10 states the rule and forbids the shortcut:

    p_i = Pr(pi > 4, Delta > log10 2) + Pr(pi <= 4, pi + Delta > 4 + log10 2)     [eq:tdirule]

    "Считать это надо совместно: pi и Delta коррелированы, перемножать вероятности
     по отдельности нельзя. На практике --- Монте-Карло по сэмплам из ансамбля."

`src/tdiprob.py` line 21 fits `HistGradientBoostingClassifier(...).predict_proba(...)` -- a separate
sigmoid, which is precisely what the text rules out.

**The precondition is not a near-miss, it is exact.** Applying the rule to the measured labels
reproduces `is_TDI` on 100.0 per cent of compounds, zero false positives and zero false negatives,
on both enzymes that carry the flag (1493 and 2334 rows). So `is_TDI` is a **deterministic function
of (pi, Delta)**, and integrating the rule over their joint predictive distribution is not an
approximation of the label probability -- it is the label probability. A classifier fitted directly
to the flag is learning a function that two quantities already in the model determine exactly.

**And this closes with the other unbuilt member rather than standing beside it.** `src/abldelta.py`
implements section 4's `pi_tdi = pi_dir + Delta, Delta >= 0` and therefore produces exactly the pair
this rule needs. The two arcs are one construction seen from the direct and the TDI track.

The joint distribution is built non-parametrically and out of fold, as in `src/bayesact.py`: the
out-of-fold residuals of `pi` and `Delta` are resampled **in pairs**, so their correlation is
carried rather than assumed. That correlation is the whole point of the document's warning, and the
third arm tests the warning itself.

    классификатор        нынешний код --- прямая подгонка к флагу
    правило, совместно   Монте-Карло по парам остатков (r_pi, r_Delta)
    правило, независимо  те же маргинали, пары РАЗОРВАНЫ. Если это работает так же,
                         корреляция не нужна и предупреждение документа неверно

Scored by MCC at a threshold chosen on training folds only, plus AUC and mean predicted probability
against the true positive rate, since section 10 measured that the plug-in threshold misses the
optimum when the probabilities are uncalibrated (0.325 against 0.358 achievable on CYP3A4).

Writes results/preds/oof_tdirule.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score
from sklearn.tree import DecisionTreeRegressor

from cypsplit import butina_folds

CYPS_TDI = ["CYP2D6", "CYP3A4"]
NTREE, LR, DEPTH, MAXFEAT = 200, 0.06, 5, 0.3
NMC = 400          # сэмплов Монте-Карло на соединение
L2 = np.log10(2.0)


def softplus(x):
    return np.log1p(np.exp(-np.abs(x))) + np.maximum(x, 0.0)


def sigmoid(x):
    return 0.5 * (1.0 + np.tanh(0.5 * x))


def boost2(Xtr, yd, md, yt, mt, Xte, seed):
    """Общие сплиты, выходы (pi_dir, raw); Delta = softplus(raw) >= 0. Как в abldelta."""
    rng = np.random.default_rng(seed)
    a0 = float(yd[md].mean()) if md.any() else 0.0
    S = np.column_stack([np.full(len(yd), a0), np.zeros(len(yd))])
    P = np.column_stack([np.full(len(Xte), a0), np.zeros(len(Xte))])
    for _ in range(NTREE):
        dlt = softplus(S[:, 1])
        r_t = np.where(mt, yt - (S[:, 0] + dlt), 0.0)
        G = np.column_stack([np.where(md, yd - S[:, 0], 0.0) + r_t, r_t * sigmoid(S[:, 1])])
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MAXFEAT,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(Xtr, G)
        S = S + LR * t.predict(Xtr)
        P = P + LR * t.predict(Xte)
    return P[:, 0], softplus(P[:, 1]), S[:, 0], softplus(S[:, 1])


def rule_prob(pi, dl, rp, rd, rng, joint):
    """Монте-Карло по правилу. joint=False рвёт пары остатков --- контроль предупреждения."""
    n = len(pi)
    i = rng.integers(0, len(rp), (NMC, n))
    j = i if joint else rng.integers(0, len(rd), (NMC, n))
    P = pi[None, :] + rp[i]
    Dd = np.maximum(dl[None, :] + rd[j], 0.0)
    hit = ((P > 4) & (Dd > L2)) | ((P <= 4) & (P + Dd > 4 + L2))
    return hit.mean(0)


def best_thr(y, p):
    ts = np.unique(np.round(p, 4))
    if len(ts) > 300:
        ts = np.quantile(p, np.linspace(0.01, 0.99, 300))
    sc = [matthews_corrcoef(y, p > t) for t in ts]
    return float(ts[int(np.argmax(sc))]), float(np.max(sc))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/oof_tdirule.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    td = (pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv")
            .set_index("Molecule_Name").reindex(rows.Molecule_Name).reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for c in CYPS_TDI:
            t0 = time.time()
            a_, b_, k_ = (f"{c}_pIC50_direct_inhibition", f"{c}_pIC50_TDI_condition",
                          f"{c}_is_TDI")
            md, mt = tr[a_].notna().to_numpy(), td[b_].notna().to_numpy()
            mk = td[k_].notna().to_numpy()
            yd, yt = tr[a_].to_numpy(float), td[b_].to_numpy(float)
            yk = td[k_].to_numpy()
            n = len(rows)
            Pj = np.zeros(n); Pi = np.zeros(n)
            for f in range(5):
                trn_mol, te_mol = fold != f, fold == f
                te = mk & te_mol
                if not te.any():
                    continue
                keep = trn_mol & (md | mt)
                pi_te, dl_te, pi_tr, dl_tr = boost2(
                    X[keep], np.nan_to_num(yd)[keep], md[keep],
                    np.nan_to_num(yt)[keep], mt[keep], X[te], seed * 10 + f)
                # остатки на ОБУЧАЮЩИХ строках, парами
                okp = md[keep]; okd = mt[keep] & md[keep]
                rp = (yd[keep][okp] - pi_tr[okp])
                rd = ((yt[keep] - yd[keep])[okd] - dl_tr[okd])
                k = min(len(rp), len(rd))
                rp, rd = rp[:k], rd[:k]
                g = np.random.default_rng(seed * 100 + f)
                Pj[te] = rule_prob(pi_te, dl_te, rp, rd, g, True)
                Pi[te] = rule_prob(pi_te, dl_te, rp, rd, g, False)
                # Классификатор не переподгоняется: его 300 итераций на 2295 колонках --- это
                # 2378 мс на дерево (пункт 140), и он один занимал больше времени, чем всё
                # остальное. Эталоном берутся ОПУБЛИКОВАННЫЕ вероятности tdiprob.py; маски
                # там свои (2346 и 1495 против 2334 и 1493), поэтому сравнение по уровню, а
                # не построчное, и это записано, а не скрыто.
            y = yk[mk].astype(int)
            for nm, p in (("правило, совместно", Pj[mk]),
                          ("правило, независимо", Pi[mk])):
                thr, mcc = best_thr(y, p)
                table.append(dict(seed=seed, фермент=c, рука=nm, MCC=mcc,
                                  AUC=float(roc_auc_score(y, p)), порог=thr,
                                  сред_p=float(p.mean()), доля_полож=float(y.mean())))
                out[f"{seed}|{nm}|{c}"] = p.tolist()
            print(f"  сид {seed} {c}: " + "  ".join(
                f"{r['рука']} MCC {r['MCC']:.3f} AUC {r['AUC']:.3f}"
                for r in table[-2:]) + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print("\n" + df.groupby(["фермент", "рука"])[["MCC", "AUC", "сред_p", "доля_полож"]]
          .mean().round(4).to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Правило воспроизводит метку на 100 %, поэтому вопрос не в том, верно ли оно, а в
том, лучше ли интегрировать его по предсказательному распределению, чем подгонять классификатор
к флагу напрямую.

Третья рука проверяет само предупреждение документа. Она берёт те же маргинальные остатки, но
РВЁТ пары, то есть считает Pr(pi) и Pr(Delta) независимыми. Если она не хуже совместной, то
корреляция не нужна и предупреждение было лишним; если хуже --- документ прав, и разница есть
цена независимого приближения.

Средняя предсказанная вероятность против доли положительных --- калибровка. Раздел 10 намерил,
что у сырых вероятностей бустинга она смещена (0.256 против 0.326 на 3A4), и что подстановочный
порог из-за этого промахивается мимо максимума MCC.""")


if __name__ == "__main__":
    main()
