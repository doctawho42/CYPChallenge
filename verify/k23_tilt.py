"""Reweight the ensemble toward pooling with one parameter, chosen out of sample.

What k22 established. The models' ranking is not the same in every similarity layer, and the
test set sits in the layer where it differs most: pooling beats the per-enzyme model by
+0.0367 of rank there against +0.0067 measured overall, and the Gaussian process is -0.0336.
Both intervals exclude zero. Our equal-weight ensemble is therefore weighted for a regime the
test set is not in.

What k22 also established. Fitting the weights per layer loses in every layer, by 0.004 to
0.019, even though the fitted weights recover exactly the right ordering (0.62 on pooling,
0.04 on the GP in the near layer). Four free parameters estimated on about 160 rows per fold
cost more in variance than the misweighting costs in bias. The effect is real; that way of
spending it is not.

So spend one parameter instead of sixteen. Slide along the segment from equal weights toward
the pooling-heavy corner,

    w(t) = (1 - t) * equal + t * corner,

and choose t by leaving a whole seed out: t is fitted on three seeds and scored on the fourth,
rotated four ways. One number, chosen on data that never sees the data it is scored on.

The affine pair is always fitted on all rows, as the submission does, and only the scoring is
restricted to a layer -- otherwise each layer would receive its own calibration, which the test
set will not get.

Three targets are reported because they answer different questions. The near layer is the
regime the test is in and the reason for doing this at all. The whole out-of-fold set is what
the log's other numbers are measured on, and a gain in the near layer bought by a loss there
is a bet on the layer analysis being right. The far layer is the control: if t helps there too,
it is not the regime that is being exploited but simply better weights.

Reads results/preds/oof_pool_all.json, oof_gp.json, oof_weak.json. Writes nothing. ~3 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json
import numpy as np
import pandas as pd
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from cypsplit import butina_folds
from shrinkchoice import fit_apply

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
SEEDS = (0, 1, 2, 3)
MEMBERS = [("oof_pool_all", "независимо"), ("oof_pool_all", "пул"),
           ("oof_gp", "GP"), ("oof_weak", "гребневая")]
MNAME = ["поферментно", "пул", "GP", "гребневая"]
EQUAL = np.array([0.25, 0.25, 0.25, 0.25])
CORNER = np.array([0.20, 0.62, 0.04, 0.13])     # что NNLS нашла в слое > 0.55 (k22)
TGRID = np.round(np.arange(0.0, 1.001, 0.05), 2)
NEAR, FAR = 0.55, 0.45


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    cache = {f: json.load(open(RES + f"preds/{f}.json"))["preds"]
             for f in {m[0] for m in MEMBERS}}
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]

    data = {}
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
            data[(seed, c)] = (
                tr.loc[m, col].to_numpy(),
                np.column_stack([np.asarray(cache[fn][f"{seed}|{arm}|{c}"], float)
                                 for fn, arm in MEMBERS]),
                tr.loc[m, col + "_conf_low"].to_numpy(),
                tr.loc[m, col + "_conf_high"].to_numpy(),
                fold[m], nn[m])
        print(f"  сид {seed} собран", flush=True)

    # Аффинная пара зависит только от (сид, фермент, t), а не от подмножества, на котором
    # потом считается метрика. Без кэша она пересчитывалась бы по три раза на комбинацию,
    # и это единственное, что тут дорого.
    cache_q = {}

    def macro(seed, t, where):
        w = (1 - t) * EQUAL + t * CORNER
        v = []
        for c in CYPS:
            y, P, lo, hi, fi, s = data[(seed, c)]
            k = (seed, c, round(float(t), 4))
            if k not in cache_q:
                cache_q[k] = fit_apply(P @ w, lo, hi, fi, np.ones(len(y)) / len(y))
            q = cache_q[k]
            sel = {"всё": np.ones(len(y), bool), "близкий": s >= NEAR,
                   "далёкий": s < FAR}[where]
            if sel.sum() < 20:
                continue
            v.append(float(strae(y[sel], q[sel], y_true_upper=hi[sel], y_true_lower=lo[sel])))
        return float(np.mean(v))

    print("\nКривая по t, среднее по всем четырём сидам (для диагностики, не для выбора):")
    print(f"{'t':>5s} {'близкий > 0.55':>15s} {'всё':>10s} {'далёкий < 0.45':>15s}")
    for t in TGRID[::2]:
        print(f"{t:5.2f} {np.mean([macro(s, t, 'близкий') for s in SEEDS]):15.4f} "
              f"{np.mean([macro(s, t, 'всё') for s in SEEDS]):10.4f} "
              f"{np.mean([macro(s, t, 'далёкий') for s in SEEDS]):15.4f}")

    print("\nЧестный выбор: t подобрано на трёх сидах, оценено на четвёртом")
    print(f"{'выбран по':>12s} {'сид':>4s} {'t':>5s} {'близкий':>10s} {'равные':>10s} {'разница':>9s}"
          f" {'всё':>10s} {'равные':>10s} {'разница':>9s}")
    got = {"близкий": [], "всё": []}
    for target in ("близкий", "всё"):
        for held in SEEDS:
            fit = [s for s in SEEDS if s != held]
            t = TGRID[int(np.argmin([np.mean([macro(s, t_, target) for s in fit])
                                     for t_ in TGRID]))]
            n_t, n_e = macro(held, t, "близкий"), macro(held, 0.0, "близкий")
            a_t, a_e = macro(held, t, "всё"), macro(held, 0.0, "всё")
            got[target].append((t, n_t - n_e, a_t - a_e))
            print(f"{target:>12s} {held:4d} {t:5.2f} {n_t:10.4f} {n_e:10.4f} {n_t-n_e:+9.4f}"
                  f" {a_t:10.4f} {a_e:10.4f} {a_t-a_e:+9.4f}")
        g = got[target]
        print(f"{'среднее':>12s} {'':4s} {np.mean([x[0] for x in g]):5.2f} {'':10s} {'':10s} "
              f"{np.mean([x[1] for x in g]):+9.4f} {'':10s} {'':10s} "
              f"{np.mean([x[2] for x in g]):+9.4f}")

    print("""
Как читать. Порог тот же --- 0.007, и отрицательное значит лучше.

Три исхода. Если t уходит в ноль на каждом отложенном сиде, наклона нет и равные веса
остаются --- это самый вероятный исход после k22. Если t стабильно ненулевое и близкий слой
выигрывает, а всё --- проигрывает, то это ставка на то, что слой правильно описывает тест;
её нужно принимать явно, а не молча. Если выигрывает и близкий слой, и всё сразу, то дело
не в режиме, а просто в том, что равные веса были плохими, и слоистая интерпретация тут
ни при чём.""")


if __name__ == "__main__":
    main()
