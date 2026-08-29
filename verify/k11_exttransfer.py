"""Train on the external labels alone, predict ours: which enzyme fails to transfer, and how.

Why. Adding the external ChEMBL rows to the training folds helps three enzymes and hurts
CYP2C9, and it hurts CYP2C9 in the arm where no scale correction is applied at all. That rules
out the first explanation offered - that the paired offset is estimated from only six shared
compounds and is therefore unreliable - because the raw arm applies no offset. It also rules
out a shortage of external labels: CYP2C9 has 2506 of them, second most of the four.

So the fault is either in their labels or in the way we merge them, and one run separates the
two. Train on external rows only and predict our compounds. No merging happens, so anything
that shows up here belongs to their data.

Three readouts, and the third is the one that decides:

  ST-RAE       - the competition metric. Comparable across enzymes because it is a ratio;
                 above 1.0 means worse than predicting our own mean;
  Spearman     - does the ORDER transfer? A model can rank compounds correctly and still miss
                 the metric entirely if it sits on the wrong scale;
  slope of y on p, and the offset - how the scales relate. A slope near 1 with a nonzero
                 offset is a shift and is fixable. A slope far below 1 is regression dilution:
                 their labels carry less signal about our assay than their spread suggests,
                 and no shift repairs that.

If CYP2C9 ranks as well as the others but scores worse, the problem is calibration and the
indicator-column arm should fix it. If CYP2C9 also ranks worst, their CYP2C9 labels are
measuring something other than what we measure, and no merging strategy repairs that. The
chemical candidate is probe dependence: CYP2C9 IC50 depends strongly on the probe substrate
(diclofenac against tolbutamide), and ChEMBL pools assays. That is a hypothesis this file
cannot test - the published file carries SMILES and four label columns and no assay metadata.

Reads data/feats.npz, data/rows.csv, the external CSVs. Prints only.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json

import numpy as np
import pandas as pd
from rdkit import Chem
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

import feats as F

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)


def canon(s):
    try:
        m = Chem.MolFromSmiles(s)
        return Chem.MolToSmiles(m) if m else None
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ext-x", required=True)
    ap.add_argument("--ext-y", required=True)
    ap.add_argument("--cache", default="", help="npz с признаками внешних, чтобы не считать заново")
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])

    ext = pd.concat([pd.read_csv(a.ext_x), pd.read_csv(a.ext_y)], axis=1)
    ext["k"] = [canon(s) for s in ext.OPENADMET_CANONICAL_SMILES]
    ours = set(filter(None, (canon(s) for s in rows.SMILES)))
    ext = ext[~ext.k.isin(ours)].reset_index(drop=True)

    cached = a.cache and _pl.Path(a.cache).exists()
    if cached:
        Xe = np.load(a.cache)["X"]
        if len(Xe) != len(ext):
            raise SystemExit(f"кэш на {len(Xe)} строк, а внешних {len(ext)}")
    else:
        dn = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
        mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
        print(f"считаю признаки для {len(ext)} внешних соединений", flush=True)
        FP, dsc, M, ok = F.build(list(ext.OPENADMET_CANONICAL_SMILES), dn, mn)
        Xe = np.hstack([FP, dsc.to_numpy(np.float32), M.to_numpy(np.float32)])
        ext = ext.iloc[ok].reset_index(drop=True)
        if a.cache:
            np.savez_compressed(a.cache, X=Xe)
    if Xe.shape[1] != X.shape[1]:
        raise SystemExit(f"ширина не совпала: {Xe.shape[1]} против {X.shape[1]}")

    # Наш собственный ранг из сохранённых предсказаний: без него столбец Spearman
    # несравним между ферментами, у них разная предсказуемость на разных соединениях.
    oof = json.load(open(RES + "preds/oof.json"))
    # Эффект от дописывания внешних строк, рука "как есть", сид 0 (src/ablext.py).
    EFF = {"CYP1A2": -0.0260, "CYP2C9": +0.0216, "CYP2D6": -0.0359, "CYP3A4": -0.0107}

    print(f"\n{'фермент':8s} {'внешних':>8s} {'наших':>7s} {'ST-RAE':>8s} {'наш rho':>8s} "
          f"{'внешн.':>7s} {'доля':>6s} {'наклон':>7s} {'сдвиг':>7s} {'эффект':>8s}")
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy()
        lo = tr.loc[m, col + "_conf_low"].to_numpy()
        hi = tr.loc[m, col + "_conf_high"].to_numpy()
        ye = ext[f"OPENADMET_LOGAC50_{c.lower()}"].to_numpy(float)
        em = ~np.isnan(ye)

        p = HistGradientBoostingRegressor(**KW).fit(Xe[em], ye[em]).predict(X[m])
        r = float(strae(y, p, y_true_upper=hi, y_true_lower=lo))
        rho = float(spearmanr(y, p).statistic)
        ours = float(spearmanr(y, np.array(oof[f"FP+DESC+MECH|{c}"])).statistic)
        b, a0 = np.polyfit(p, y, 1)
        print(f"{c:8s} {int(em.sum()):8d} {int(m.sum()):7d} {r:8.4f} {ours:8.3f} "
              f"{rho:7.3f} {rho / ours:6.2f} {b:7.3f} {a0:7.3f} {EFF[c]:+8.4f}", flush=True)

    print("""
Как читать. Абсолютные значения не важны: модель ни разу не видела нашего ассея, и ST-RAE
здесь заведомо хуже нашего собственного. Важен столбец ДОЛЯ --- ранг внешней модели, делённый
на наш собственный на тех же соединениях. Голый Spearman между ферментами несравним, потому
что они предсказуемы по-разному и измерены на разных молекулах; отношение это снимает.

Что получилось. CYP2C9 не выделяется ни рангом, ни метрикой: по голому Spearman он второй из
четырёх. Гипотеза о том, что их метки на 2C9 мешаные, этой проверкой не поддержана. Зато доля
упорядочивает четыре фермента ровно так же, как эффект от дописывания строк: 2D6 переносится
лучше всех и выигрывает больше всех, 2C9 хуже всех и единственный проигрывает. Порядок
совпадает целиком, а точка безубыточности лежит между 0.66 и 0.70.

Чего это НЕ доказывает. Точек четыре. Идеальное совпадение порядка при случайной расстановке
выпадает с вероятностью 1/24. Это гипотеза, а не результат, и держится она на предсказании:
если поднять перенос на CYP2C9, знак эффекта обязан перевернуться. Столбец-индикатор источника
(src/ablsrc.py) --- первая попытка его поднять.""")


if __name__ == "__main__":
    main()
