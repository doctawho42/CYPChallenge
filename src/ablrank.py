"""Two-view low-rank completion: the enzyme subspace from the dense view, the scores from the sparse one.

The structure the data actually has, and it has been in plain sight. The potency matrix is 4905 by
4 and it is observed **twice, differently**:

    кривые     разреженно и точно      6525 ячеек из 19620
    скрининг   плотно и шумно         17504 ячейки, все четыре фермента у 4375 молекул

Density and precision are complementary in exactly the way that makes low-rank completion work. The
dense view identifies the **column space** -- how the four enzymes relate -- because it sees every
enzyme for every screened molecule. The sparse view identifies **scale and the per-molecule scores**
because its entries are fitted curves rather than single points.

Why a rank restriction is not the proteochemometrics item 168 closed by arithmetic. That argument
was about the **input**: with four enzymes all seen in training, a one-hot indicator is a sufficient
statistic for enzyme identity, so any protein-descriptor vector is a change of basis on it and
`src/ablcoord.py` measured exactly that at +0.0007. A rank constraint on the **output** is not a
change of basis, it is a *restriction*, and restrictions help under sparse labels even when the
unrestricted model is perfectly identifiable. That is the ordinary reduced-rank regression argument
and this table is its textbook case: four correlated outputs, two thirds of the cells missing.

And it is the formal version of a mechanism this file already measured. Item 132 established that
pooling works through **contrast** -- the enzyme indicator lets the trees condition on which enzyme
a row is -- and item 152 that the curve table carries all four enzymes for **41** molecules while
the screen carries all four for **4375**. A rank-r factorisation is precisely the object whose
factors are "level" and "contrast", and the screen is where contrast is measured a hundredfold more
densely.

The construction. With `V` the 4-by-r decoder and `f_1..f_r` the latent models,

    y_ie ~= sum_k V_ek f_k(x_i),   и градиент по f_k равен sum_{e наблюдено} V_ek (y_ie - yhat_ie)

so missingness needs no imputation: a molecule contributes through whichever enzymes it has. `V`
comes from the top-r right singular vectors of the **calibrated, centred screening matrix**, fitted
on training-fold molecules only.

Arms, and the third is what makes the second readable:

    независимо           четыре модели, эталон
    ранг 4               V = единичная; ограничения нет, обязана воспроизвести эталон
    ранг 2, V скрининг   предложение
    ранг 2, V случайная  ортонормальная V того же ранга. Отделяет "помогает ОГРАНИЧЕНИЕ ранга"
                         от "помогает ИМЕННО подпространство, измеренное скринингом"
    ранг 3, V скрининг

Pre-registered. Rank 4 must reproduce the per-enzyme reference, or the harness is wrong rather than
the idea. Rank 2 with the screen's decoder must beat rank 2 with a random decoder, or the subspace
is not what is doing the work. And the gain should be largest where labels are fewest -- CYP2C9 at
1285 -- since a restriction buys most where there is least to fit on.

Writes results/preds/oof_rank.json.
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
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression
from sklearn.tree import DecisionTreeRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
NTREE, LR, DEPTH, MAXFEAT = 200, 0.06, 5, 0.3      # пины из src/ablpairloss.py
FLOOR = {"CYP1A2": 0.0061, "CYP2C9": 0.0071, "CYP2D6": 0.0049, "CYP3A4": 0.0033}


def decoder_from_screen(S, ok, r, seed, random=False):
    """Верхние r правых сингулярных векторов калиброванной центрированной матрицы скрининга."""
    if random:
        g = np.random.default_rng(seed)
        Q, _ = np.linalg.qr(g.normal(size=(len(CYPS), len(CYPS))))
        return Q[:, :r]
    # НЕ центрировать. Центрирование убирает направление (1,1,1,1), то есть общий уровень
    # потенции --- доминирующую компоненту, без которой ранг 2 не может представить даже
    # среднюю активность. Дымовой прогон с центрированием дал 0.3137 против эталонных
    # 0.5320. На нецентрированной матрице первый сингулярный вектор и есть уровень, второй
    # --- контраст, что ровно разложение пункта 132.
    A = S[ok]
    _, sv, Vt = np.linalg.svd(A, full_matrices=False)
    return Vt[:r].T, sv


def boost_latent(Xtr, Y, M, V, Xte, seed):
    """r латентных бустингов с общим линейным декодером V (4 x r).

    Пропуски не заполняются: ячейка без метки просто не входит в сумму градиента, поэтому
    молекула вносит вклад через те ферменты, которые у неё есть.
    """
    r = V.shape[1]
    rng = np.random.default_rng(seed)
    # Вес молекулы для латентной модели k --- диагональ гауссо-ньютоновского гессиана,
    # sum_e M_ie * V_ek^2: сколько информации молекула вообще несёт об этой латенте.
    # Без него молекула без метки по ферменту e входит в критерий сплита НУЛЁМ и стягивает
    # дерево к неделению. Дымовой прогон поймал это как расхождение руки "ранг 4" с
    # эталоном на 0.042 --- при V = I вес сводится к маске, и ранг 4 становится тождествен
    # поферментной модели, чем он и обязан быть.
    W = M.astype(float) @ (V ** 2)                      # (n, r)
    # Одинаковый ПОЛНЫЙ бюджет деревьев у всех рук: ранг r получает 4/r раундов, иначе
    # сравнение идёт про число деревьев, а не про ограничение.
    rounds = max(1, int(round(NTREE * len(CYPS) / r)))
    Fi = np.zeros((len(Xtr), r))
    Fe = np.zeros((len(Xte), r))
    # старт: наименьшие квадраты на наблюдённых средних по ферментам
    mu = np.array([Y[M[:, e], e].mean() if M[:, e].any() else 0.0 for e in range(len(CYPS))])
    f0 = np.linalg.lstsq(V, mu, rcond=None)[0]
    Fi += f0; Fe += f0
    for _ in range(rounds):
        Yh = Fi @ V.T
        R = np.where(M, Y - Yh, 0.0)                 # остаток только по наблюдённым
        G = R @ V                                    # (n, r): градиент по латентным
        for k in range(r):
            t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MAXFEAT,
                                      random_state=int(rng.integers(1 << 30)))
            wk = W[:, k]
            if not np.any(wk > 0):
                continue
            t.fit(Xtr, np.where(wk > 0, G[:, k] / np.maximum(wk, 1e-9), 0.0),
                  sample_weight=wk)
            Fi[:, k] += LR * t.predict(Xtr)
            Fe[:, k] += LR * t.predict(Xte)
    return Fe @ V.T


def boost1(Xtr, ytr, Xte, seed):
    rng = np.random.default_rng(seed)
    s = np.full(len(ytr), float(ytr.mean())); pred = np.full(len(Xte), s[0])
    for _ in range(NTREE):
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MAXFEAT,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(Xtr, ytr - s)
        s = s + LR * t.predict(Xtr); pred = pred + LR * t.predict(Xte)
    return pred


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--arms", default="независимо|ранг 4|ранг 2, V скрининг|"
                                      "ранг 2, V случайная|ранг 3, V скрининг")
    ap.add_argument("--out", default=RES + "preds/oof_rank.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    arms = a.arms.split("|")

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)
    n = len(rows)

    Y = np.zeros((n, 4)); M = np.zeros((n, 4), bool); L2 = np.full((n, 4), np.nan)
    LO, HI = {}, {}
    for e, c in enumerate(CYPS):
        col = f"{c}_pIC50_direct_inhibition"
        M[:, e] = tr[col].notna().to_numpy()
        Y[M[:, e], e] = tr.loc[M[:, e], col].to_numpy()
        LO[c] = tr[col + "_conf_low"].to_numpy(float)
        HI[c] = tr[col + "_conf_high"].to_numpy(float)
        sub = sc[sc.enzyme == c].set_index("Molecule_Name")["log2fc_estimate"]
        L2[:, e] = sub.reindex(rows.Molecule_Name).to_numpy(float)
    dense = np.isfinite(L2).all(1)
    print(f"молекул {n}; помечено по ферментам {M.sum(0).tolist()}; "
          f"скрининг по ВСЕМ четырём у {int(dense.sum())}; кривые по всем четырём у "
          f"{int(M.all(1).sum())}\n")

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for arm in arms:
            t0 = time.time()
            r = {"seed": seed, "рука": arm}
            P = {c: np.zeros(n) for c in CYPS}
            for f in range(5):
                trn_mol, te_mol = fold != f, fold == f
                if arm == "независимо":
                    for e, c in enumerate(CYPS):
                        a_, b_ = M[:, e] & trn_mol, M[:, e] & te_mol
                        if b_.any():
                            P[c][b_] = boost1(X[a_], Y[a_, e], X[b_], seed * 10 + f)
                    continue
                rk = int(arm.split("ранг")[1].split(",")[0].strip()) if "ранг" in arm else 4
                if arm == "ранг 4":
                    V = np.eye(4)
                elif "случайная" in arm:
                    V = decoder_from_screen(None, None, rk, seed * 10 + f, random=True)
                else:
                    # Калибровка log2fc -> шкала pIC50, поферментно, ТОЛЬКО на обучающих.
                    S = np.zeros((n, 4))
                    for e, c in enumerate(CYPS):
                        fitm = M[:, e] & np.isfinite(L2[:, e]) & trn_mol
                        iso = IsotonicRegression(increasing=False, out_of_bounds="clip")
                        iso.fit(L2[fitm, e], Y[fitm, e])
                        S[:, e] = iso.predict(np.nan_to_num(L2[:, e], nan=0.0))
                    V, sv = decoder_from_screen(S, dense & trn_mol, rk, seed)
                    if f == 0:
                        print(f"    {arm}: сингулярные числа скрининга "
                              f"{np.round(sv / sv.sum(), 3).tolist()}", flush=True)
                keep = trn_mol & M.any(1)
                Pe = boost_latent(X[keep], Y[keep], M[keep], V, X[te_mol], seed * 10 + f)
                for e, c in enumerate(CYPS):
                    P[c][te_mol] = Pe[:, e]

            for e, c in enumerate(CYPS):
                m = M[:, e]
                y, lo, hi, fi = Y[m, e], LO[c][m], HI[c][m], fold[m]
                p = P[c][m]
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{arm}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for t in ("пара", "rho"):
                r[f"MACRO {t}"] = round(float(np.mean([r[f"{c} {t}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {arm:22s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    g = df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                       + [f"{c} rho" for c in CYPS]].mean()
    print("\n" + g.round(4).to_string())
    if "независимо" in g.index:
        print(f"\n{'рука':22s} " + "".join(f"{c[3:]+' Δ':>12s}" for c in CYPS))
        for arm in g.index:
            if arm == "независимо":
                continue
            line = f"{arm:22s} "
            for c in CYPS:
                d = g.loc[arm, f"{c} rho"] - g.loc["независимо", f"{c} rho"]
                line += f"{d:+11.4f}{'*' if abs(d) > FLOOR[c] else ' '}"
            print(line)
        print("  * --- больше СОБСТВЕННОГО пола фермента (пункт 165)")
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Порядок обязателен. «ранг 4» неограничен и обязан лечь около «независимо» ---
если нет, сломан стенд, а не идея. Дальше «ранг 2, V скрининг» против «ранг 2, V случайная»:
если они равны, работает само ограничение ранга, а измеренное скринингом подпространство ни
при чём, и тогда вся конструкция про экономию параметров, а не про ферменты. И только если
скрининговая V выигрывает у случайной, можно говорить, что плотный вид дал столбцовое
пространство.

Сингулярные числа печатаются на первом фолде: они показывают, сколько структуры ферментов
вообще есть в ранге 2 --- если первое число близко к единице, все четыре фермента почти
пропорциональны, и рангу 2 нечего добавлять к рангу 1.""")


if __name__ == "__main__":
    main()
