"""Is the layer effect real, and does gating the ensemble on similarity pay.

k20 measured each model's rank inside four layers of nearest-neighbour similarity and found
two things that, if real, matter more than anything else currently in the queue:

    layer          < 0.35   0.35-0.45   0.45-0.55   > 0.55
    pool - per-enzyme    -0.0137     +0.0220     +0.0189     +0.0366
    GP   - per-enzyme    +0.0048     -0.0033     -0.0212     -0.0342

The test set sits at a median similarity of 0.587, in the right-hand column. Pooling was
logged at +0.0067 overall; in the regime it will actually be scored in it is five times that.
The Gaussian process runs the other way and is worst exactly where the test lives.

But the right-hand layer holds about 200 compounds per enzyme per seed, and a Spearman on 200
points has a standard error near 0.06. Four seeds and four enzymes are not sixteen independent
estimates -- same compounds, correlated splits. So the differences above are eyeballed, not
established, and this file establishes them or does not.

Two parts.

  1. Paired bootstrap over compounds inside each layer, resampling the compounds and not the
     seeds, so the interval reflects the thing that is actually small.

  2. Whether acting on it pays. Ensemble weights fitted per layer, with the fitting done
     leave-one-fold-out inside the layer so that no compound contributes to the weights used
     to score it. Compared against the equal weights the submission currently uses.

Part 2 is the one that can change the submission, and it is deliberately given the harder
protocol: equal weights cost nothing to justify, fitted weights have to earn the degrees of
freedom they spend.

Reads results/preds/oof_pool_all.json, oof_gp.json, oof_weak.json. Writes nothing. ~8 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from scipy.optimize import nnls
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from cypsplit import butina_folds
from shrinkchoice import fit_apply

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
LAYERS = [(0.0, 0.35), (0.35, 0.45), (0.45, 0.55), (0.55, 1.01)]
LAB = ["< 0.35", "0.35-0.45", "0.45-0.55", "> 0.55"]
SEEDS = (0, 1, 2, 3)
MEMBERS = [("oof_pool_all", "независимо"), ("oof_pool_all", "пул"),
           ("oof_gp", "GP"), ("oof_weak", "гребневая")]
MNAME = ["поферментно", "пул", "GP", "гребневая"]
B = 1000


def _rho(a, b):
    """Spearman without scipy's per-call overhead: Pearson on ranks."""
    ra = np.argsort(np.argsort(a)).astype(np.float64)
    rb = np.argsort(np.argsort(b)).astype(np.float64)
    ra -= ra.mean(); rb -= rb.mean()
    d = np.sqrt((ra * ra).sum() * (rb * rb).sum())
    return float((ra * rb).sum() / d) if d else 0.0


def main():
    rng = np.random.default_rng(0)
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    cache = {f: json.load(open(RES + f"preds/{f}.json"))["preds"]
             for f in {m[0] for m in MEMBERS}}
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]

    # Собираем всё в один список наблюдений, чтобы бутстрап шёл по соединениям.
    pool = {k: [] for k in range(4)}      # слой -> список (y, preds[4], lo, hi, fold)
    for seed in SEEDS:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        nn = np.zeros(len(rows))
        for f in range(5):
            te, trn = np.where(fold == f)[0], np.where(fold != f)[0]
            ref = [fps[j] for j in trn]
            for i in te:
                nn[i] = max(DataStructs.BulkTanimotoSimilarity(fps[i], ref))
        for c in CYPS:
            col = f"{c}_pIC50_direct_inhibition"
            m = tr[col].notna().to_numpy()
            y = tr.loc[m, col].to_numpy()
            lo = tr.loc[m, col + "_conf_low"].to_numpy()
            hi = tr.loc[m, col + "_conf_high"].to_numpy()
            P = np.column_stack([np.asarray(cache[fn][f"{seed}|{arm}|{c}"], float)
                                 for fn, arm in MEMBERS])
            v, fi = nn[m], fold[m]
            for k, (a, b) in enumerate(LAYERS):
                s = (v >= a) & (v < b)
                if s.sum() < 30:
                    continue
                pool[k].append((y[s], P[s], lo[s], hi[s], fi[s]))
        print(f"  сид {seed} собран", flush=True)

    print("\n1. Разрыв рангов внутри слоя, парный бутстрап по соединениям "
          f"({B} пересэмплирований)")
    print(f"{'сравнение':28s}" + "".join(f"{l:>20s}" for l in LAB))
    for i, j, tag in [(1, 0, "пул - поферментно"), (2, 0, "GP - поферментно"),
                      (3, 0, "гребневая - поферментно")]:
        cells = []
        for k in range(4):
            d = []
            for y, P, *_ in pool[k]:
                n = len(y)
                idx = rng.integers(0, n, size=(B, n))
                # Спирмен на пересэмплированных индексах, парно: те же индексы обеим
                # моделям, иначе разность мерила бы шум пересэмплирования, а не модели.
                di = np.empty(B)
                for b in range(B):
                    ii = idx[b]
                    di[b] = (_rho(y[ii], P[ii, i]) - _rho(y[ii], P[ii, j]))
                d.append(di)
            d = np.mean(d, axis=0)          # среднее по ферментам и сидам, внутри реплики
            lo95, hi95 = np.percentile(d, [2.5, 97.5])
            cells.append(f"{d.mean():+.4f} [{lo95:+.3f},{hi95:+.3f}]")
        print(f"{tag:28s}" + "".join(f"{c:>20s}" for c in cells))

    print("\n2. Веса ансамбля, подогнанные внутри слоя (подгонка без выбывшего фолда)")
    print(f"{'слой':>12s} {'равные веса':>14s} {'подогнанные':>14s} {'разница':>10s}   средние веса")
    for k in range(4):
        eq, ft, ws = [], [], []
        for y, P, lo, hi, fi in pool[k]:
            pe = np.zeros(len(y)); pf = np.zeros(len(y))
            for f in np.unique(fi):
                te, trn = fi == f, fi != f
                if trn.sum() < 20 or te.sum() == 0:
                    continue
                pe[te] = P[te].mean(1)
                # Неотрицательные веса, сумма свободна: масштаб всё равно перепишет пара.
                w, _ = nnls(P[trn], y[trn])
                if w.sum() == 0:
                    w = np.ones(P.shape[1])
                pf[te] = P[te] @ w
                ws.append(w / w.sum())
            u = np.ones(len(y)) / len(y)
            eq.append(float(strae(y, fit_apply(pe, lo, hi, fi, u),
                                  y_true_upper=hi, y_true_lower=lo)))
            ft.append(float(strae(y, fit_apply(pf, lo, hi, fi, u),
                                  y_true_upper=hi, y_true_lower=lo)))
        a, b = np.mean(eq), np.mean(ft)
        wm = np.mean(ws, axis=0)
        print(f"{LAB[k]:>12s} {a:14.4f} {b:14.4f} {b-a:+10.4f}   "
              + " ".join(f"{n} {x:.2f}" for n, x in zip(MNAME, wm)))

    print("""
Как читать. Часть 1: интервал, накрывающий ноль, означает, что слой ничего не доказал ---
их справа мало, и это ожидаемо. Смотреть надо, растёт ли оценка слева направо монотонно:
одна значимая клетка слабее, чем согласованный тренд по четырём.

Часть 2: подогнанные веса тратят четыре степени свободы на слой. Отрицательная разница
меньше 0.007 --- это не выигрыш, а шум, и равные веса остаются. Выигрыш только в правом
слое при нулях в остальных был бы самым сильным исходом: он означал бы, что гейт по
близости нужен, а не просто другие веса.""")


if __name__ == "__main__":
    main()
