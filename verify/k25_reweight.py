"""Re-score every saved ablation under weights that match the test set's regime.

The case for doing this is now measured rather than argued. Item 97 found cross-validation runs
at a nearest-neighbour similarity of 0.435 and the test set sits at 0.587. Items 107 and 108
found that the models' ranking is genuinely different there -- pooling +0.0367 against +0.0067
overall, the Gaussian process -0.0336 -- with bootstrap intervals excluding zero. So at least one
conclusion in this log was reached in a regime it will not be scored in, and the only way to know
how many is to re-score them all.

Why this shift is identifiable and the label shift is not. `src/covshift.py` estimates a shift in
the *label*, which nobody observes on the test set, so it carries a posterior and a model-choice
systematic that item 94 measured at 0.11 to 0.51. Here the shifted variable is similarity to the
training set, which is computed from structures alone. The test structures are published. There
is nothing to infer: p_test(s) is counted, not estimated.

The one confound, and the fix. An out-of-fold row's similarity is a maximum over four fifths of
the training set; a test row's is a maximum over all of it. A maximum over a larger reference set
is larger for free, so comparing the two directly would manufacture a shift. Test similarities are
therefore taken against random four-fifths subsets, averaged over draws, matching the reference
size the out-of-fold rows had. This is the same confound as the one caught in item 87's split
control and it is caught the same way.

What is reported. For every arm of every saved ablation: rank and pair, weighted and unweighted,
and whether the ordering of the arms inside that file changes. A file whose ordering is unchanged
is a file whose conclusion is safe. Effective sample size is printed for every enzyme, because
importance weights buy relevance with variance and a re-scoring at an ESS of two hundred is not
evidence of anything.

Reads every results/preds/oof*.json listed in FILES. Writes nothing. ~10 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json
import numpy as np
import pandas as pd
from scipy.stats import rankdata
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from cypsplit import butina_folds
from shrinkchoice import fit_apply

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
SEED = 0
NDRAW = 5            # подвыборок 4/5 для теста
BINS = np.array([0.0, .25, .30, .35, .40, .45, .50, .55, .60, .70, 1.01])
CLIP = 8.0           # потолок веса: без него один редкий бин съедает всё
NBOOT = 300

FILES = {
    "признаки (ablate)":   ("oof", ["FP", "FP+DESC", "FP+DESC+MECH"]),
    "пул (ablpool)":       ("oof_pool_all", ["независимо", "пул"]),
    "GP (gp)":             ("oof_gp", ["GP"]),
    "слабые (ablweak)":    ("oof_weak", ["лес", "гребневая", "kNN"]),
    "форма (ablshape)":    ("oof_shape", ["без формы", "с формой"]),
    "агрегаты (ablagg)":   ("oof_agg", ["база", "флаг", "вес по структуре", "вес со скринингом"]),
    "FCFP (ablfcfp)":      ("oof_fcfp", ["ECFP (как сейчас)", "FCFP вместо ECFP", "оба"]),
}


def wspearman(y, p, w):
    """Weighted Spearman: weighted Pearson on the unweighted ranks."""
    a, b = rankdata(y).astype(float), rankdata(p).astype(float)
    sw = w.sum()
    a -= (w * a).sum() / sw
    b -= (w * b).sum() / sw
    d = np.sqrt((w * a * a).sum() * (w * b * b).sum())
    return float((w * a * b).sum() / d) if d else 0.0


def wstrae(y, q, lo, hi, w, rng):
    """ST-RAE under weights, via weighted resampling so the organisers' code is used as is."""
    pr = w / w.sum()
    n = len(y)
    v = [strae(y[i], q[i], y_true_upper=hi[i], y_true_lower=lo[i])
         for i in (rng.choice(n, size=n, p=pr) for _ in range(NBOOT))]
    return float(np.mean(v))


def main():
    rng = np.random.default_rng(0)
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    te_df = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]
    fte = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in te_df.SMILES
           if Chem.MolFromSmiles(s) is not None]

    fold, _ = butina_folds(list(rows.SMILES), seed=SEED)
    s_oof = np.zeros(len(rows))
    for f in range(5):
        te, trn = np.where(fold == f)[0], np.where(fold != f)[0]
        ref = [fps[j] for j in trn]
        for i in te:
            s_oof[i] = max(DataStructs.BulkTanimotoSimilarity(fps[i], ref))

    # Тестовые сходства против случайных 4/5 --- тот же размер опоры, что был у OOF.
    n_sub = int(round(0.8 * len(fps)))
    s_te = np.zeros(len(fte))
    for _ in range(NDRAW):
        ref = [fps[j] for j in rng.choice(len(fps), n_sub, replace=False)]
        for t, f_ in enumerate(fte):
            s_te[t] += max(DataStructs.BulkTanimotoSimilarity(f_, ref)) / NDRAW

    print(f"сходство: OOF медиана {np.median(s_oof):.3f}, тест медиана {np.median(s_te):.3f} "
          f"(опора {n_sub} молекул в обоих случаях)")
    h_o = np.histogram(s_oof, BINS)[0] / len(s_oof)
    h_t = np.histogram(s_te, BINS)[0] / len(s_te)
    ratio = np.where(h_o > 0, h_t / np.maximum(h_o, 1e-9), 1.0).clip(0, CLIP)
    print(f"\n{'бин':>12s} {'доля OOF':>10s} {'доля тест':>11s} {'вес':>8s}")
    for k in range(len(BINS) - 1):
        print(f"{BINS[k]:5.2f}-{BINS[k+1]:5.2f} {h_o[k]:10.3f} {h_t[k]:11.3f} {ratio[k]:8.2f}")
    wrow = ratio[np.clip(np.digitize(s_oof, BINS) - 1, 0, len(ratio) - 1)]

    print(f"\n{'фермент':8s} {'n':>6s} {'ESS':>7s} {'ESS/n':>8s}")
    W = {}
    for c in CYPS:
        m = tr[f"{c}_pIC50_direct_inhibition"].notna().to_numpy()
        w = wrow[m]
        ess = w.sum() ** 2 / (w * w).sum()
        W[c] = w
        print(f"{c:8s} {m.sum():6d} {ess:7.0f} {100*ess/m.sum():7.1f} %")

    print("\n" + "=" * 78)
    for label, (fn, arms) in FILES.items():
        try:
            preds = json.load(open(RES + f"preds/{fn}.json"))
            preds = preds["preds"] if "preds" in preds else preds
        except FileNotFoundError:
            print(f"\n{label}: файла нет, пропуск"); continue
        res = {}
        for arm in arms:
            got = {}
            for tag, wgt in (("равномерно", None), ("под тест", True)):
                rk, pr = [], []
                for c in CYPS:
                    col = f"{c}_pIC50_direct_inhibition"
                    m = tr[col].notna().to_numpy()
                    y = tr.loc[m, col].to_numpy()
                    lo = tr.loc[m, col + "_conf_low"].to_numpy()
                    hi = tr.loc[m, col + "_conf_high"].to_numpy()
                    for k in (f"{SEED}|{arm}|{c}", f"{arm}|{c}"):
                        if k in preds:
                            p = np.asarray(preds[k], float); break
                    else:
                        p = None
                    if p is None or len(p) != len(y):
                        rk = None; break
                    w = (W[c] if wgt else np.ones(len(y)))
                    w = w / w.sum()
                    q = fit_apply(p, lo, hi, fold[m], w)
                    rk.append(wspearman(y, p, w))
                    pr.append(wstrae(y, q, lo, hi, w, rng) if wgt
                              else float(strae(y, q, y_true_upper=hi, y_true_lower=lo)))
                if rk is None:
                    got = None; break
                got[tag] = (float(np.mean(pr)), float(np.mean(rk)))
            if got:
                res[arm] = got
        if not res:
            print(f"\n{label}: ключи не совпали, пропуск"); continue
        print(f"\n{label}")
        print(f"  {'рука':22s} {'пара равн.':>11s} {'ранг равн.':>11s} "
              f"{'пара тест':>11s} {'ранг тест':>11s}")
        for arm, g in res.items():
            print(f"  {arm:22s} {g['равномерно'][0]:11.4f} {g['равномерно'][1]:11.4f} "
                  f"{g['под тест'][0]:11.4f} {g['под тест'][1]:11.4f}")
        o_u = sorted(res, key=lambda a: -res[a]["равномерно"][1])
        o_w = sorted(res, key=lambda a: -res[a]["под тест"][1])
        print(f"  порядок по рангу: равномерно {o_u}")
        print(f"                    под тест   {o_w}"
              + ("   <-- ПОРЯДОК ИЗМЕНИЛСЯ" if o_u != o_w else "   (тот же)"))

    print("""
Как читать. ESS решает, читать ли вообще: при ESS/n ниже примерно трети веса режут выборку
сильнее, чем стоит любой вывод из неё.

Файл, в котором порядок рук не изменился, --- закрытый вопрос: его вывод вынесен в режиме,
который для этого вывода неважен. Файл, в котором порядок изменился, придётся пересматривать,
и это относится к документу, а не только к коду.

Величины между колонками не сравнивать: взвешенная ST-RAE считается пересэмплированием и у
неё свой шум. Сравнивать надо порядок и знак разностей внутри колонки.""")


if __name__ == "__main__":
    main()
