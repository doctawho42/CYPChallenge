"""Is any of the screening channel's value plate structure rather than chemistry.

Item 13 measured the plate component at 0.9 to 8 per cent of residual variance after subtracting
what the chemistry explains, and concluded no plate term was needed. That was a decomposition
argument. This is the harder version of the same question: block the cross-validation by plate, so
that a held-out compound's plate is never in training, and see whether the screening channel still
pays.

The channel to test is the one that actually moved. `src/ablcontrast.py` measured the screen's
*level* at +0.0040 of rank and its contrast at -0.0015; the level is the arm with something to
lose, so it is the arm blocked here.

Why plate blocking is the right stress and reweighting is not. If part of the level's value is
that compounds sharing a plate share a systematic offset, an ordinary fold split leaks it: the
same plate appears on both sides. Blocking removes exactly that channel and nothing else, so the
difference between the two protocols is the plate component of the signal rather than an estimate
of it.

Folds are blocked per enzyme, because a molecule sits on four different plates -- one per enzyme --
and there is no single plate assignment for a compound. That means these folds are not the
repository's Butina folds and the numbers are not comparable with the rest of the log; both
protocols are recomputed here so the comparison is internal.

Reads data/feats.npz, rows.csv and the screening table. Writes nothing. ~20 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, tutorial
tutorial()

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor

from cypsplit import butina_folds

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
KW_SCR = dict(max_iter=200, learning_rate=0.06, max_leaf_nodes=31,
              l2_regularization=1.0, random_state=0)


def plate_folds(plate, nf=5, seed=0):
    """Каждый планшет целиком в одном фолде. Планшеты раскладываются жадно по размеру,
    чтобы фолды не разъехались по количеству строк."""
    u, cnt = np.unique(plate[plate != ""], return_counts=True)
    order = u[np.argsort(-cnt)]
    load = np.zeros(nf)
    where = {}
    for p_ in order:
        k = int(np.argmin(load))
        where[p_] = k
        load[k] += cnt[list(u).index(p_)]
    return np.array([where.get(p_, 0) for p_ in plate])


def run(X, y, fold, level, use_level):
    p = np.zeros(len(y))
    XX = np.hstack([X, level[:, None].astype(np.float32)]) if use_level else X
    for f in np.unique(fold):
        trn, te = fold != f, fold == f
        if te.sum() == 0 or trn.sum() < 50:
            continue
        p[te] = HistGradientBoostingRegressor(**KW).fit(XX[trn], y[trn]).predict(XX[te])
    return p


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")

    print(f"{'фермент':8s} {'протокол':>12s} {'без уровня':>11s} {'с уровнем':>10s} {'вклад':>9s}")
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy()
        sub = sc[sc.enzyme == c].set_index("Molecule_Name")
        l2 = rows.set_index("Molecule_Name").join(sub["log2fc_estimate"])["log2fc_estimate"]
        pl = rows.set_index("Molecule_Name").join(sub["plate_id"])["plate_id"].fillna("")
        lvl = np.nan_to_num(l2.to_numpy(float))
        Xi, yi = X[m], y

        for tag, fld in (("кластерный", butina_folds(list(rows.SMILES))[0][m]),
                         ("по планшетам", plate_folds(pl.to_numpy()[m]))):
            a = spearmanr(yi, run(Xi, yi, fld, lvl[m], False)).statistic
            b = spearmanr(yi, run(Xi, yi, fld, lvl[m], True)).statistic
            print(f"{c:8s} {tag:>12s} {a:11.4f} {b:10.4f} {b-a:+9.4f}", flush=True)

    print("""
Как читать. Сравнивать надо ВКЛАД уровня между двумя протоколами, а не сами ранги: у
блокированных по планшетам фолдов другая геометрия, и абсолютные числа с кластерными
несопоставимы.

Если вклад держится --- канал несёт химию, и пункт 13 подтверждён более жёсткой проверкой.
Если падает --- часть того, что мы приняли за скрининговый сигнал, была планшетной
систематикой, которую обычный сплит пропускает через себя.""")


if __name__ == "__main__":
    main()
