"""The pooled arms crossed a threshold in the learner and quietly changed learner.

The finding, which is arithmetic and not a hypothesis. The pinned scikit-learn defaults
`early_stopping='auto'`, and 'auto' means **on above 10000 samples and off below**. The pooled
tables in `src/ablpool.py` sit on both sides of that line:

    рука                всего строк   обучающих (4/5)   ранняя остановка
    независимо (2D6)           1493              1194   нет
    пул                        6525              5220   нет
    пул+TDI                   13063             10450   ДА
    пул+скрининг              18030             14424   ДА

So item 125 -- "pool+TDI adds nothing despite doubling the table" -- compared an arm trained on
all its rows for all 300 iterations against an arm that held out ten per cent of its rows as a
validation set and stopped whenever that set stopped improving. Two learners, one comparison. The
margin is 450 rows, 4.5 per cent over the threshold, which is why nobody saw it.

This is not a claim that item 125's conclusion is wrong. It is a claim that the measurement did
not test what it was written to test, and the difference matters because `пул+TDI` was the arm
that closed "more supervision helps" -- a question `src/ablaux.py` and `src/ablncgc.py` both
reopen with much larger tables, both of which would cross the same line.

What is measured here, on one fold by default because the full re-run is four hours per arm:

    n_iter_ и do_early_stopping_   сработала ли остановка и на какой итерации
    предсказания обоих режимов     насколько они расходятся на отложенном фолде
    ранг и метрика на этом фолде   в какую сторону

If `n_iter_` comes back at 300, early stopping was armed but never fired, the two learners agree
and item 125 stands as written. If it comes back well below 300, the arm was truncated, and the
comparison has to be redone -- `--folds 5` does the whole thing.

Reads data/feats.npz, data/rows.csv and the TDI table. Writes nothing.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--folds", type=int, default=1, help="сколько фолдов считать (5 = полностью)")
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    td = (pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv")
            .set_index("Molecule_Name").reindex(rows.Molecule_Name).reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)
    fold, _ = butina_folds(list(rows.SMILES))

    Y, LO, HI, M, YT, MT = {}, {}, {}, {}, {}, {}
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        M[c] = tr[col].notna().to_numpy()
        Y[c] = tr[col].to_numpy(float)
        LO[c] = tr[col + "_conf_low"].to_numpy(float)
        HI[c] = tr[col + "_conf_high"].to_numpy(float)
        tc = f"{c}_pIC50_TDI_condition"
        MT[c] = td[tc].notna().to_numpy()
        YT[c] = td[tc].to_numpy(float)

    n_dir = int(sum(M[c].sum() for c in CYPS))
    n_tdi = int(sum(MT[c].sum() for c in CYPS))
    print(f"порог early_stopping='auto' в этой версии scikit-learn: 10000 строк\n")
    print(f"{'рука':14s} {'всего':>7s} {'обучающих':>10s}  режим по умолчанию")
    for nm, n in (("пул", n_dir), ("пул+TDI", n_dir + n_tdi)):
        t = int(round(0.8 * n))
        print(f"{nm:14s} {n:7d} {t:10d}  {'ранняя остановка ВКЛ' if t > 10000 else 'выкл'}")
    print()

    def stack(keep, with_tdi):
        Xs, ys = [], []
        for e, c in enumerate(CYPS):
            for cond, (mask, lab) in enumerate(((M[c], Y[c]), (MT[c], YT[c]))):
                if cond and not with_tdi:
                    continue
                sel = mask & keep
                if not sel.any():
                    continue
                ind = np.zeros((int(sel.sum()), 5), np.float32)
                ind[:, e] = 1.0
                ind[:, 4] = cond
                Xs.append(np.hstack([X[sel], ind]))
                ys.append(lab[sel])
        return np.vstack(Xs), np.concatenate(ys)

    for with_tdi, armname in ((False, "пул"), (True, "пул+TDI")):
        for es in ("auto", False):
            t0 = time.time()
            p_all = {c: np.zeros(len(rows)) for c in CYPS}
            info = []
            for f in range(a.folds):
                Xtr, ytr = stack(fold != f, with_tdi)
                mdl = HistGradientBoostingRegressor(**KW, early_stopping=es).fit(Xtr, ytr)
                info.append((int(mdl.n_iter_), bool(mdl.do_early_stopping_), len(ytr)))
                for e, c in enumerate(CYPS):
                    te = M[c] & (fold == f)
                    if not te.any():
                        continue
                    ind = np.zeros((int(te.sum()), 5), np.float32)
                    ind[:, e] = 1.0
                    p_all[c][te] = mdl.predict(np.hstack([X[te], ind]))
            it, doing, ntr = info[0]
            print(f"{armname:8s} early_stopping={str(es):5s}  обучающих {ntr:6d}  "
                  f"do_early_stopping_={doing}  n_iter_={it}/300  ({time.time()-t0:.0f} с)",
                  flush=True)

            sel = np.isin(fold, range(a.folds))
            rr, pp = [], []
            for c in CYPS:
                m = M[c] & sel
                if m.sum() < 20:
                    continue
                y, lo, hi = Y[c][m], LO[c][m], HI[c][m]
                p = p_all[c][m]
                q = fit_apply(p, lo, hi, np.zeros(int(m.sum()), int),
                              np.ones(int(m.sum())) / int(m.sum()))
                rr.append(spearmanr(y, p).statistic)
                pp.append(strae(y, q, y_true_upper=hi, y_true_lower=lo))
            print(f"{'':8s} на этих фолдах: rho {np.mean(rr):.4f}  пара {np.mean(pp):.4f}\n")

    print("""
Как читать. Первое число --- n_iter_. Если у режима 'auto' он равен 300, остановка была
взведена и не сработала, две руки совпадают, и пункт 125 стоит как записан. Если он заметно
меньше 300, рука «пул+TDI» училась меньшим числом итераций и на 90 % своих строк, а
сравнение с «пулом» надо переделывать: --folds 5.

Отдельно стоит сравнить строку 'auto' и строку False у ОДНОЙ руки. Если они дают разный
ранг, то любая будущая рука, переходящая 10000 строк --- а это и ablaux, и ablncgc ---
обязана ставить early_stopping явно, иначе размер таблицы будет менять учителя заодно
с тем, что проверяется.""")


if __name__ == "__main__":
    main()
