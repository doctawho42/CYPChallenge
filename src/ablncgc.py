"""The NCGC panel as extra training rows, with a source indicator and no calibration.

What the panel is and why it is not ChEMBL. Items 60 to 63 measured why the supplied external
set failed: ChEMBL is a sample selected by *publishability*, and selection is not correctable by
a shift. The NCGC panel is a whole-library qHTS campaign -- one laboratory, one protocol, sixteen
concentrations, everything run -- so nothing was filtered by outcome. `src/ncgcmerge.py` folded
the six assays into 54177 rows over 13126 structures.

Two measured facts decide the design, and the second overturned the obvious plan.

**The panel does not touch the blinded test set.** Zero of 750 test molecules appear in it, so
there is no question about using it.

**Calibrating it against our labels is not possible, and must not be attempted.** Only 11 to 48
molecules carry both a fitted NCGC AC50 and one of our curves:

    фермент   общих   rho    сдвиг NCGC-наш   sd разности
    CYP1A2       32  0.72            +0.762         0.577
    CYP2C9       11  0.40            +0.448         0.710
    CYP2D6       48  0.74            +0.442         0.484
    CYP3A4       21  0.49            +0.868         0.843

The offset is *larger* than the ChEMBL offset the project rejected, and eleven molecules at
sd 0.71 give it a standard error of 0.21 -- half the quantity being estimated. So the arm that
subtracts a fitted shift is not built here at all. What is built is item 63's actual lesson: a
**source indicator**, which lets the model learn the offset from all 54177 rows instead of from
eleven, and lets it learn it as an interaction with chemistry rather than as a constant.

The leak that has to be closed by hand. 183 to 319 panel structures are also our molecules. If
one of them sits in the held-out fold, its NCGC measurement is a measurement of a test row, and
letting it into training would be the ordinary kind of leak. Panel rows are therefore mapped to
our folds by canonical SMILES and dropped whenever their molecule is held out. Panel structures
we do not have are unconstrained and always trainable.

Censored rows are kept as a separate arm rather than dropped by default. A qHTS compound with no
fitted AC50 is a measurement -- "not active up to the top concentration" -- and discarding those
7258 CYP3A4 rows would rebuild by hand exactly the selection that killed the ChEMBL merge. Their
target is the lowest observed pAC50 for that enzyme, which is what the assay can say and no more.

Arms:

    база                        pooled curve rows only, 6525. The reference.
    +NCGC активные              + 18711 rows with a fitted AC50 on our four enzymes.
    +NCGC активные, перемешанный  the same rows, targets permuted within enzyme. Identical row
                                count and target marginal, molecule-to-measurement link broken.
                                If this scores the same, the gain was extra rows, not the panel.
    +NCGC активные+2C19         adds 3698 CYP2C19 rows on a fifth indicator column. The challenge
                                does not score CYP2C19, so this asks whether a related isoform we
                                will never be graded on still teaches the shared trunk -- the
                                only form in which a proteochemometric coordinate is testable
                                here, since our own matrix has four enzymes and no fifth.
    +NCGC всё, с цензурой       + the 35466 censored rows.

Pre-registration. The panel is almost entirely chemistry we do not have -- 13126 structures
overlapping ours in about 300 -- so if it helps, it should help most where our own coverage is
thinnest, and the effect should survive the permutation. A gain that does not survive the
permutation is regularisation. A gain on CYP2C19's arm alone, with no movement on the four
scored enzymes, would say the fifth isoform teaches nothing transferable.

The learner is the plain-tree booster pinned to `src/ablpairloss.py`, for the reason given in
`src/ablaux.py`: HistGradientBoostingRegressor needs 2378 ms per tree here, and these tables run
to 60000 rows.

Reads data/ncgc/panel.csv, data/ncgc/feats.npz, data/feats.npz, data/rows.csv.
Writes results/preds/oof_ncgc.json.
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
from sklearn.tree import DecisionTreeRegressor
from rdkit import Chem, RDLogger
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
ALL5 = CYPS + ["CYP2C19"]
NTREE, LR, DEPTH, MAXFEAT = 200, 0.06, 5, 0.3     # пины из src/ablpairloss.py


def canon(s):
    m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
    if m is None:
        return None
    fr = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
    if len(fr) > 1:
        m = max(fr, key=lambda f: f.GetNumHeavyAtoms())
    return Chem.MolToSmiles(m)


def boost(Xtr, ytr, Xte, seed, wtr=None):
    rng = np.random.default_rng(seed)
    base = float(ytr.mean() if wtr is None else np.average(ytr, weights=wtr))
    s = np.full(len(ytr), base)
    pred = np.full(len(Xte), base)
    for _ in range(NTREE):
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MAXFEAT,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(Xtr, ytr - s, sample_weight=wtr)
        s = s + LR * t.predict(Xtr)
        pred = pred + LR * t.predict(Xte)
    return pred


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--arms", default="база|+NCGC активные|+NCGC активные, перемешанный|"
                                      "+NCGC активные+2C19|+NCGC всё, с цензурой")
    ap.add_argument("--out", default=RES + "preds/oof_ncgc.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    arms = a.arms.split("|")

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)

    zn = np.load(D + "ncgc/feats.npz", allow_pickle=True)
    XN = np.hstack([zn["FP"], zn["DESC"], zn["MECH"]]).astype(np.float32)
    nsm = list(zn["SMILES"])
    pos = {s: i for i, s in enumerate(nsm)}
    assert XN.shape[1] == X.shape[1], "ширина матриц панели и обучения не совпала"

    panel = pd.read_csv(D + "ncgc/panel.csv", low_memory=False)
    panel = panel[panel.smiles.isin(pos)].copy()
    # Цензурированной строке ставится наименьшее наблюдённое pAC50 фермента: это то,
    # что прибор может сказать, и ни на десятую больше.
    floor = panel.groupby("enzyme").pAC50.min()
    panel["target"] = panel.pAC50.fillna(panel.enzyme.map(floor))
    panel["row"] = panel.smiles.map(pos)

    # --- отображение панели на наши фолды, чтобы отложенная молекула не попала в обучение ---
    ours = {}
    for i, s in enumerate(rows.SMILES):
        c = canon(s)
        if c is not None:
            ours.setdefault(c, i)
    panel["ourid"] = panel.smiles.map(ours)          # NaN = структуры у нас нет
    print(f"панель: {len(panel)} строк, {panel.smiles.nunique()} структур; "
          f"из них наших — {int(panel.ourid.notna().sum())} строк "
          f"({panel.loc[panel.ourid.notna(), 'smiles'].nunique()} структур)")
    print(f"{'фермент':9s} {'с AC50':>8s} {'цензур.':>8s}   порог pAC50")
    for e in ALL5:
        d = panel[panel.enzyme == e]
        print(f"{e:9s} {int(d.pAC50.notna().sum()):8d} {int(d.pAC50.isna().sum()):8d}   "
              f"{floor.get(e, float('nan')):.2f}")
    print()

    Y, LO, HI, M = {}, {}, {}, {}
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        M[c] = tr[col].notna().to_numpy()
        Y[c] = tr[col].to_numpy(float)
        LO[c] = tr[col + "_conf_low"].to_numpy(float)
        HI[c] = tr[col + "_conf_high"].to_numpy(float)

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for arm in arms:
            t0 = time.time()
            r = {"seed": seed, "рука": arm}
            use = arm != "база"
            shuffle = "перемешанный" in arm
            with_c19 = "2C19" in arm
            with_cens = "цензурой" in arm
            # Вес строк панели. При 1.0 панель --- 74 % обучающей таблицы, и чужой протокол
            # со сдвигом +0.44..+0.87 доминирует над функцией потерь, а индикатор источника
            # вынужден отыгрывать это из меньшинства. Пункт 153.
            w_pan = float(arm.split("w")[1].split()[0]) if " w" in arm else 1.0
            enz_list = ALL5 if with_c19 else CYPS
            p_all = {c: np.zeros(len(rows)) for c in CYPS}

            for f in range(5):
                trn_mol, te_mol = fold != f, fold == f
                rng = np.random.default_rng(seed * 100 + f)
                Xs, ys = [], []
                ws = []
                for e, c in enumerate(CYPS):
                    sel = M[c] & trn_mol
                    ind = np.zeros((int(sel.sum()), 6), np.float32)
                    ind[:, e] = 1.0
                    Xs.append(np.hstack([X[sel], ind]))
                    ys.append(Y[c][sel])
                    ws.append(np.ones(int(sel.sum())))

                if use:
                    for e, c in enumerate(enz_list):
                        d = panel[panel.enzyme == c]
                        if not with_cens:
                            d = d[d.pAC50.notna()]
                        # Наша молекула в отложенном фолде -> её измерение выбрасывается.
                        held = d.ourid.notna() & d.ourid.map(
                            lambda i: te_mol[int(i)] if i == i else False)
                        d = d[~held]
                        if d.empty:
                            continue
                        t = d.target.to_numpy(float)
                        if shuffle:
                            t = t[rng.permutation(len(t))]
                        ind = np.zeros((len(d), 6), np.float32)
                        ind[:, e] = 1.0
                        ind[:, 5] = 1.0                     # источник: панель
                        Xs.append(np.hstack([XN[d.row.to_numpy()], ind]))
                        ys.append(t)
                        ws.append(np.full(len(d), w_pan))

                Xtr, ytr, wtr = np.vstack(Xs), np.concatenate(ys), np.concatenate(ws)
                for e, c in enumerate(CYPS):
                    te = M[c] & te_mol
                    if not te.any():
                        continue
                    ind = np.zeros((int(te.sum()), 6), np.float32)
                    ind[:, e] = 1.0
                    p_all[c][te] = boost(Xtr, ytr, np.hstack([X[te], ind]),
                                         seed * 10 + f, wtr)

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
            print(f"  сид {seed} {arm:30s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                         + [f"{c} rho" for c in CYPS]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Сравнение с перемешанной рукой --- главное и единственное, которое отделяет
панель от лишних строк: там поровну строк, признаков и маргинала мишени, разорвана только
связь структуры с её измерением.

Калибровки здесь нет намеренно. Сдвиг NCGC к нашим меткам +0.44..+0.87 подгонять не на чем
(11..48 общих молекул), поэтому смещение отдано индикатору источника: модель учит его на всех
строках сразу и как взаимодействие с химией, а не как константу из одиннадцати точек.

Рука с CYP2C19 спрашивает то, чего наша матрица спросить не может: помогает ли изоформа,
которую нам никогда не будут оценивать. Движение по четырём зачётным ферментам от её
добавления --- единственная здесь проверяемая форма протеохемометрики.""")


if __name__ == "__main__":
    main()
