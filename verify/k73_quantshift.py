"""Квантовый блок против СДВИГА TDI --- цели, против которой его ни разу не мерили.

Пункт 186 прогнал этот блок против ПОТЕНЦИИ: четыре сида, перестановочный контроль, ни один
фермент не проходит своего пола, настоящий блок неотличим от перемешанного. Закрыт. Но целью
там был pIC50 прямого ингибирования.

Пункт 238 установил, что весь остаток классификационного трека --- это СДВИГ
Delta = pi_TDI - pi_dir, и что физика у него другая: не электростатика обратимого связывания, а
реакционная способность, то есть превращает ли фермент лиганд в реакционноспособную частицу.
Ровно эти величины в блоке и лежат: `fukui_minus`, HOMO, LUMO, щель. Признаки посчитаны и лежат
на диске для всех 4905 обучающих и 750 тестовых молекул, ноль NaN.

**Прекон пункта 166 отработал и дал жёлтый свет.** Девять колонок из десяти восстанавливаются из
DESC+MECH вне фолда при R^2 > 0.5 (медиана 0.68; `n_arom_N` 0.97, `q_basicN` 0.76, `fukui_minus`
0.55), и только `dipole` не восстанавливается (0.14). Это не уровень третичного амина из пункта
236, где R^2 был 0.998 и колонка буквально уже лежала в матрице, но и не чистая новая величина.

Поэтому вопрос разделён на два плеча, и это главное в файле:

    +квант            сырой блок: помогает ли он вообще
    +квант ОСТАТОК    блок минус его вневыборочная реконструкция из DESC+MECH

Если помогает сырой, а остаток нет --- выигрыш есть переупаковка уже имеющегося, и строить нечего.
Если помогает остаток --- это та часть, которой в матрице не было, и тогда есть что обсуждать.

Плюс перестановочный контроль пункта 186 (блок, перемешанный по строкам) и контроль ширины: 10
колонок гауссова шума, чтобы отделить эффект от простого расширения матрицы на десять столбцов.

Мишени: Delta на обоих ферментах. Метрики --- R^2 и Спирмен по Delta, плюс MCC метки через
правило, потому что курс обмена R^2 в MCC в проекте неизвестен и его лучше не предполагать.
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
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import RidgeCV
from sklearn.metrics import matthews_corrcoef

from cypsplit import butina_folds

KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
L = np.log10(2.0)
GATE = 4.0 + L


def clipstd(B, trn):
    A = np.nan_to_num(B.astype(np.float64), posinf=0.0, neginf=0.0)
    mu, sd = A[trn].mean(0), A[trn].std(0) + 1e-9
    return np.clip((A - mu) / sd, -5.0, 5.0)


def quantum_residual(Q, base, fold):
    """Q минус его вневыборочная реконструкция из base. Обрезка +-5 обязательна: Ipc
    доходит до 5e14 и без неё гребневая разлетается (R^2 до -2500)."""
    R = np.zeros_like(Q)
    for j in range(Q.shape[1]):
        p = np.zeros(len(Q))
        for f in range(5):
            trn, te = fold != f, fold == f
            Z = clipstd(base, trn)
            p[te] = RidgeCV(alphas=np.logspace(-1, 4, 12)).fit(Z[trn], Q[trn, j]).predict(Z[te])
        R[:, j] = Q[:, j] - p
    return R


def oof(X, y, fold):
    p = np.zeros(len(y))
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0:
            continue
        p[te] = HistGradientBoostingRegressor(**KW).fit(X[trn], y[trn]).predict(X[te])
    return p


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "preds/quantshift.json")
    a = ap.parse_args()

    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    base = np.hstack([z["DESC"], z["MECH"]])
    q = np.load(D + "quantum.npz", allow_pickle=True)
    Q = q["train"].astype(np.float64)
    rows = pd.read_csv(D + "rows.csv")
    tdi = pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv").set_index("Molecule_Name").reindex(rows.Molecule_Name)
    inh = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv").set_index("Molecule_Name").reindex(rows.Molecule_Name)

    recs = []
    for seed in [int(s) for s in a.seeds.split(",")]:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        print(f"\n=== сид {seed} ===", flush=True)
        t0 = time.time()
        QR = quantum_residual(Q, base, fold)
        print(f"  остаток посчитан, {time.time()-t0:.0f} с; доля оставшейся дисперсии по колонкам "
              + " ".join(f"{QR[:,j].var()/Q[:,j].var():.2f}" for j in range(Q.shape[1])), flush=True)
        rng = np.random.default_rng(seed)
        arms = {
            "база": X,
            "+квант": np.hstack([X, Q]),
            "+квант ОСТАТОК": np.hstack([X, QR]),
            "+квант перемешан": np.hstack([X, Q[rng.permutation(len(Q))]]),
            "+10 колонок шума": np.hstack([X, rng.normal(size=(len(X), 10))]),
        }
        for c in ["CYP3A4", "CYP2D6"]:
            lab = tdi[f"{c}_is_TDI"]
            ta = tdi[f"{c}_pIC50_TDI_condition"].to_numpy(float)
            dr = inh[f"{c}_pIC50_direct_inhibition"].to_numpy(float)
            m = (lab.notna() & np.isfinite(ta) & np.isfinite(dr)).to_numpy()
            y = lab[m].astype(bool).to_numpy()
            dl = (ta - dr)[m]
            gate_t = (ta > GATE)[m]
            fi = fold[m]
            print(f"  {c}: n {len(y)}, sd(Delta) {dl.std():.3f}", flush=True)
            for k, XX in arms.items():
                t0 = time.time()
                p = oof(XX[m], dl, fi)
                r2 = 1 - np.mean((p - dl) ** 2) / np.var(dl)
                rho = spearmanr(dl, p).statistic
                mcc = matthews_corrcoef(y, gate_t & (p > L))
                print(f"      {k:18s} R^2(Delta) {r2:+.4f}  rho {rho:+.4f}  "
                      f"MCC(ворота ист.+сдвиг предск.) {mcc:+.4f}  {time.time()-t0:.0f} с",
                      flush=True)
                recs.append({"seed": seed, "cyp": c, "arm": k, "r2": float(r2),
                             "rho": float(rho), "mcc": float(mcc)})
        json.dump(recs, open(a.out, "w"))

    df = pd.DataFrame(recs)
    print("\n\nсреднее по сидам")
    for met in ("r2", "rho", "mcc"):
        print(f"\n{met}:")
        print(df.pivot_table(index="arm", columns="cyp", values=met).round(4).to_string())
    print("\nпротив базы, знак по 8 клеткам:")
    piv = df.pivot_table(index=["seed", "cyp"], columns="arm", values="rho")
    for k in piv.columns:
        if k == "база":
            continue
        d = piv[k] - piv["база"]
        print(f"  {k:18s} rho {d.mean():+.4f}, положительных {int((d>0).sum())}/{len(d)}")
    print(f"\nсохранено: {a.out}", flush=True)


if __name__ == "__main__":
    main()
