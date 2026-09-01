"""The screening table as a training target, on the 11505 rows where no curve exists.

Why this is not item 136 again, which is the first objection and has to be answered before
anything is run. `src/ablcontrast.py` fed the model the level and contrast **predicted from
structure** as extra columns, and measured -0.0015 of rank. Item 136 read that as closing the
idea. It does not, because that design had no branch that could have succeeded:

    непредсказуемо из структуры -> канал есть шум          -> ноль
    предсказуемо из структуры   -> канал есть функция X    -> избыточен -> ноль

Both outcomes give a null, so the null carried no information about the screen. What was never
handed to the model is the screen's **measurements** -- and there are a great many of them:

    фермент   кривых   скрининг без кривой
    CYP1A2      1412           2963
    CYP2C9      1285           3090
    CYP2D6      1493           2882
    CYP3A4      2335           2570
    итого       6525          11505

1.76 times the labelled table, on compounds the labelled table never mentions for that enzyme.

The second reason is the one that makes this the natural next arm rather than one more idea.
Item 132 settled that pooling works by **contrast** -- the enzyme indicator lets the trees
condition splits on which enzyme a row is, and the centred pool (0.4647) does not recover the
blind pool's loss (0.4885). But the curve table carries all four enzymes for **41** molecules.
The screen carries all four for **4375**. If contrast is the mechanism, the place where contrast
is actually measured has never been in the training table at all.

Design. One pooled table with the enzyme indicator, exactly as `src/ablpool.py` builds it, plus
a source column, and two kinds of row:

    кривая   target = pIC50, source = 0
    скрининг target = isotonic(log2fc) fitted out of fold, source = 1

Screen rows are added only where that enzyme has no curve for that molecule, so no cell of the
matrix carries two competing targets and the added rows are new supervision rather than a
reweighting of what is already there.

The calibration is fitted **on training folds only**, on molecules that have both a curve and a
reading. Fitting it on everything would carry held-out labels into the training targets through
the map, which is the quiet version of the leak this file has caught twice.

Arms, and the control is the point:

    независимо                    four models, curve rows only. The reference.
    пул                           pooled curve rows. Must reproduce pooling's gain.
    пул+скрининг                  pooled, curve + 11505 screen rows.
    пул+скрининг, перемешанный    the same rows, the same targets, permuted among molecules
                                  within each enzyme. Row count, feature distribution and target
                                  marginal are identical; only the molecule-to-measurement link
                                  is destroyed. If this scores like the arm above it, the gain is
                                  regularisation from extra rows and not the measurement.
    пул+скрининг w0.3             screen rows down-weighted, since a single point at 49.5 uM is
                                  a far weaker statement than a fitted curve.

Pre-registered readings, per enzyme, so the answer cannot be chosen afterwards. The gain should
be **largest where the screen adds most relative to the curves** -- CYP2C9 (2.4x), CYP1A2 (2.1x),
CYP2D6 (1.9x), CYP3A4 (1.1x, because item 129's 530-compound campaign was never screened). If a
gain appears instead concentrated on CYP3A4, it is not new supervision and something else is
happening. If rank does not move on any enzyme while the permuted control does not move either,
the screen's information is genuinely orthogonal to rank and the 0.354 gap of item 134 is a
property of the bands rather than of what the screen knows.

Two implementation facts that are corrections rather than choices.

**`early_stopping=False` is set explicitly.** The pinned scikit-learn defaults it to `'auto'`,
which turns early stopping ON above 10000 samples. `пул` trains on 5220 rows and `пул+TDI` on
10450 -- so `src/ablpool.py` compared an arm without early stopping against an arm with it, plus
a 10 per cent validation holdout, and item 125's null is confounded. Any pooled arm here would
cross the same threshold at 14424.

**The learner is the plain-tree booster of `src/ablpairloss.py`, not HistGradientBoostingRegressor.**
Not for accuracy: because at 2378 ms per tree HistGB needs about fourteen hours per arm on this
table, and the same trees at max_features=0.3 need forty minutes while scoring 0.5687 against
0.5651 at seed 0. The constants are pinned to that file's; the reference row here therefore lands
near 0.5687, not the journal's 0.5651, and the arms are comparable to each other rather than to
the ablation table.

Reads data/feats.npz, data/rows.csv and the single-concentration screen.
Writes results/preds/oof_aux.json.
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
# Пины скопированы из src/ablpairloss.py и должны совпадать с ним, иначе числа двух
# скриптов несравнимы. Копия, а не общий модуль --- дом строит скрипты самодостаточными.
NTREE, LR, DEPTH, MAXFEAT = 200, 0.06, 5, 0.3


def boost(Xtr, ytr, wtr, Xte, seed, maxfeat=MAXFEAT):
    """Бустинг над обычными деревьями, квадратичная потеря, с весами строк."""
    rng = np.random.default_rng(seed)
    base = float(np.average(ytr, weights=wtr))
    s = np.full(len(ytr), base)
    pred = np.full(len(Xte), base)
    for _ in range(NTREE):
        tree = DecisionTreeRegressor(max_depth=DEPTH, max_features=maxfeat,
                                     random_state=int(rng.integers(1 << 30)))
        tree.fit(Xtr, ytr - s, sample_weight=wtr)
        s = s + LR * tree.predict(Xtr)
        pred = pred + LR * tree.predict(Xte)
    return pred


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--arms", default="независимо|пул|пул+скрининг|"
                                      "пул+скрининг, перемешанный|пул+скрининг w0.3")
    ap.add_argument("--out", default=RES + "preds/oof_aux.json")
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

    # Метки, полосы, маски и показания скрининга --- всё в порядке rows.csv.
    Y, LO, HI, M, L2, S = {}, {}, {}, {}, {}, {}
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        M[c] = tr[col].notna().to_numpy()
        Y[c] = tr[col].to_numpy(float)
        LO[c] = tr[col + "_conf_low"].to_numpy(float)
        HI[c] = tr[col + "_conf_high"].to_numpy(float)
        sub = sc[sc.enzyme == c].set_index("Molecule_Name")["log2fc_estimate"]
        L2[c] = sub.reindex(rows.Molecule_Name).to_numpy(float)
        S[c] = ~np.isnan(L2[c])

    print(f"{'фермент':8s} {'кривых':>8s} {'скрин':>8s} {'скрин без кривой':>18s} {'rho(log2fc, y)':>16s}")
    for c in CYPS:
        both = M[c] & S[c]
        rho = spearmanr(L2[c][both], Y[c][both]).statistic
        print(f"{c:8s} {M[c].sum():8d} {S[c].sum():8d} {(S[c] & ~M[c]).sum():18d} {rho:16.3f}")
    print("\nОтрицательная rho ожидается: сильнее ингибирует --- ниже log2fc. "
          "Изотоника подгоняется убывающей.\n")

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for arm in arms:
            t0 = time.time()
            r = {"seed": seed, "рука": arm}
            w_scr = 0.3 if "w0.3" in arm else 1.0
            use_scr = "скрининг" in arm
            shuffle = "перемешанный" in arm
            pooled = "пул" in arm
            # mf<x> в имени руки переопределяет подвыборку колонок. Нужно потому, что
            # индикатор фермента --- ОДНА колонка из 2300, и при max_features=0.3 он
            # отсутствует в 70 % решений о сплите. Для пулированной модели, весь механизм
            # которой в обусловливании по индикатору, это может быть решающим.
            mf = float(arm.split("mf")[1].split()[0]) if "mf" in arm else MAXFEAT
            p_all = {c: np.zeros(n) for c in CYPS}

            for f in range(5):
                trn_mol, te_mol = fold != f, fold == f
                rng = np.random.default_rng(seed * 100 + f)

                # --- сборка обучающей таблицы ---
                Xs, ys, ws = [], [], []
                for e, c in enumerate(CYPS):
                    sel = M[c] & trn_mol                      # строки с кривой
                    ind = np.zeros((int(sel.sum()), 5), np.float32)
                    ind[:, e] = 1.0                           # source = 0 в столбце 4
                    Xs.append(np.hstack([X[sel], ind]))
                    ys.append(Y[c][sel])
                    ws.append(np.ones(int(sel.sum())))

                    if not use_scr:
                        continue
                    # Калибровка log2fc -> pIC50 ТОЛЬКО на обучающих молекулах с обоими.
                    fitm = M[c] & S[c] & trn_mol
                    addm = S[c] & ~M[c] & trn_mol
                    if fitm.sum() < 50 or addm.sum() == 0:
                        continue
                    iso = IsotonicRegression(increasing=False, out_of_bounds="clip")
                    iso.fit(L2[c][fitm], Y[c][fitm])
                    t = iso.predict(L2[c][addm])
                    if shuffle:
                        # Контроль: те же строки и то же маргинальное распределение
                        # мишени, разорвана только связь «молекула --- измерение».
                        t = t[rng.permutation(len(t))]
                    ind = np.zeros((int(addm.sum()), 5), np.float32)
                    ind[:, e] = 1.0
                    ind[:, 4] = 1.0                           # source = 1: строка скрининга
                    Xs.append(np.hstack([X[addm], ind]))
                    ys.append(t)
                    ws.append(np.full(int(addm.sum()), w_scr))

                if pooled:
                    Xtr = np.vstack(Xs); ytr = np.concatenate(ys); wtr = np.concatenate(ws)
                    for e, c in enumerate(CYPS):
                        te = M[c] & te_mol
                        if not te.any():
                            continue
                        ind = np.zeros((int(te.sum()), 5), np.float32)
                        ind[:, e] = 1.0                       # предсказываем как кривую
                        p_all[c][te] = boost(Xtr, ytr, wtr,
                                             np.hstack([X[te], ind]), seed * 10 + f, mf)
                else:
                    # Поферментная рука. Со скринингом --- строки скрининга ТОГО ЖЕ фермента,
                    # без индикатора фермента (он был бы константой), но с колонкой источника.
                    # Это отделяет вклад скрининга от пулированной конструкции целиком.
                    for e, c in enumerate(CYPS):
                        trn, te = M[c] & trn_mol, M[c] & te_mol
                        if not te.any():
                            continue
                        Xa = [np.hstack([X[trn], np.zeros((int(trn.sum()), 1), np.float32)])]
                        ya = [Y[c][trn]]
                        wa = [np.ones(int(trn.sum()))]
                        if use_scr:
                            fitm = M[c] & S[c] & trn_mol
                            addm = S[c] & ~M[c] & trn_mol
                            if fitm.sum() >= 50 and addm.sum() > 0:
                                iso = IsotonicRegression(increasing=False, out_of_bounds="clip")
                                iso.fit(L2[c][fitm], Y[c][fitm])
                                t = iso.predict(L2[c][addm])
                                if shuffle:
                                    t = t[rng.permutation(len(t))]
                                Xa.append(np.hstack([X[addm],
                                                     np.ones((int(addm.sum()), 1), np.float32)]))
                                ya.append(t)
                                wa.append(np.full(int(addm.sum()), w_scr))
                        p_all[c][te] = boost(
                            np.vstack(Xa), np.concatenate(ya), np.concatenate(wa),
                            np.hstack([X[te], np.zeros((int(te.sum()), 1), np.float32)]),
                            seed * 10 + f, mf)

            for c in CYPS:
                m = M[c]
                y, lo, hi, fi = Y[c][m], LO[c][m], HI[c][m], fold[m]
                p = p_all[c][m]
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{arm}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for tag in ("пара", "rho"):
                r[f"MACRO {tag}"] = round(float(np.mean([r[f"{c} {tag}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {arm:28s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                         + [f"{c} rho" for c in CYPS]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Три сравнения, и порядок обязателен.

«пул» против «независимо» --- воспроизводится ли выигрыш пулирования на этом учителе.
Если нет, дальше читать нечего: механизм контраста здесь не включился.

«пул+скрининг» против «пул» --- добавляют ли 11505 измерений что-нибудь сверх кривых.

«пул+скрининг» против «пул+скрининг, перемешанный» --- ЭТО ГЛАВНОЕ. Строк, признаков и
маргинала мишени в двух руках поровну; разорвана только связь молекулы с её измерением.
Разность между ними --- ровно то, что знает скрининг и не знает структура. Если она нулевая,
а над «пулом» есть выигрыш, то выигрыш дала регуляризация лишними строками, и говорить о
скрининге нельзя.

По ферментам выигрыш предрегистрирован убывающим: 2C9 (скрининга в 2.4 раза больше кривых),
1A2 (2.1), 2D6 (1.9), 3A4 (1.1). Выигрыш, сосредоточенный на 3A4, опровергает прочтение
«новая супервизия».""")


if __name__ == "__main__":
    main()
