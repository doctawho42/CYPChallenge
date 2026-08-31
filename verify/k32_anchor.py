"""Validation built to the test set's design, on the one population that permits it.

Why this exists now and not this morning. Earlier today the anchor split was declared
unconstructible: 93.6 per cent of training molecules are Butina singletons, and at the test's own
similarity of 0.587 only 574 of 4905 have any neighbour at all. Item 129 then found why the
remainder is shaped as it is -- the training set is a diversity screen of 4375 singletons glued to
a **CYP3A4-only analog campaign of 530 compounds**, two thirds of which do have a close neighbour.
The split is impossible in general and possible on that campaign, which is also the only
sub-population whose geometry matches the test's.

What the test's design is, and what the ordinary split does instead. The test is 750 compounds
forming about 132 analog series around parents that sit **inside** the training set, at a median
similarity of 0.587. Butina clustering at 0.35 does the exact opposite: it assigns a whole cluster
to one fold, so a held-out compound's analogs are held out with it and it faces the model with no
close relative. Our cross-validation is therefore not a weak version of the test's regime, it is
the mirror image of it, and that is why no reweighting reaches it -- item 123 measured the ceiling
at 26 per cent effective sample size and item 128's chi-square forbids the rest.

The construction here is the mirror of the mirror. Inside the campaign, group compounds by
single-linkage at 0.587; in each group keep one member in training as the **anchor** and hold out
the others as its analogs. A held-out compound then has a labelled near neighbour in training,
exactly as a test compound does. The anchor is rotated so that every member takes the role, which
both removes the arbitrariness of picking one and gives a spread to report.

What it can and cannot say. It reproduces the test's *series* geometry and nothing else: the test
carries four enzymes per compound and the campaign carries one, so this is silent on the
four-enzyme structure and speaks only for CYP3A4. That is one enzyme out of four, and it is one
more than we had.

Two protocols are scored on the same held-out compounds so the comparison is like for like:

    кластерный   the ordinary Butina fold this repository uses -- analogs held out together
    якорный      anchor kept in training, analogs held out

If an intervention's benefit is the same under both, its measurement was regime-independent and
the log's number for it stands. If it differs, we have been measuring it in the mirror image of
the regime it will be scored in.

Reads data/feats.npz and data/rows.csv. Writes nothing. ~10 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from cypsplit import butina_folds

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
THR = 0.587       # медиана сходства тестового соединения с его родителем
ROT = 3           # сколько раз повернуть роль якоря


def groups_at(fps, idx, thr):
    """Одинарная связность внутри idx: группы аналогов."""
    pos = {int(i): k for k, i in enumerate(idx)}
    parent = list(range(len(idx)))

    def find(a):
        while parent[a] != a:
            parent[a] = parent[parent[a]]
            a = parent[a]
        return a

    for k, i in enumerate(idx):
        sims = np.array(DataStructs.BulkTanimotoSimilarity(fps[i], [fps[j] for j in idx]))
        for kk in np.where(sims >= thr)[0]:
            if kk == k:
                continue
            a, b = find(k), find(int(kk))
            if a != b:
                parent[a] = b
    out = {}
    for k in range(len(idx)):
        out.setdefault(find(k), []).append(k)
    return [v for v in out.values() if len(v) >= 2]


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    inscr = rows.Molecule_Name.isin(set(sc.Molecule_Name)).to_numpy()
    z = np.load(D + "feats.npz")
    FP, DESC, MECH = z["FP"], z["DESC"], z["MECH"]
    X_nm = np.hstack([FP, DESC]).astype(np.float32)
    X_full = np.hstack([FP, DESC, MECH]).astype(np.float32)

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]

    col = "CYP3A4_pIC50_direct_inhibition"
    lab = tr[col].notna().to_numpy()
    camp = np.where((~inscr) & lab)[0]
    print(f"кампания: {len(camp)} размеченных соединений вне скрининговой библиотеки")
    grp = groups_at(fps, camp, THR)
    n_in = sum(len(g) for g in grp)
    print(f"групп аналогов при {THR}: {len(grp)}, в них {n_in} соединений "
          f"(медиана размера {int(np.median([len(g) for g in grp]))}, максимум {max(len(g) for g in grp)})")
    print(f"для сравнения, тест: 132 родителя, медиана 6 аналогов\n")

    y = tr[col].to_numpy(float)
    lo = tr[col + "_conf_low"].to_numpy(float)
    hi = tr[col + "_conf_high"].to_numpy(float)
    fold, _ = butina_folds(list(rows.SMILES))

    ARMS = {"FP+DESC": X_nm, "FP+DESC+MECH": X_full}
    acc = {}
    for rot in range(ROT):
        held, anchors = [], []
        for g in grp:
            a = g[rot % len(g)]
            anchors.append(camp[a])
            held += [camp[k] for k in g if k != a]
        held = np.array(held)
        keep = np.ones(len(rows), bool)
        keep[held] = False
        # Подтверждение, что у отложенного действительно есть близкий сосед в обучении.
        simn = []
        for i in held:
            v = np.array(DataStructs.BulkTanimotoSimilarity(
                fps[i], [fps[j] for j in np.where(keep & lab)[0]]))
            simn.append(v.max())
        for nm, XX in ARMS.items():
            trn = keep & lab
            p_anch = HistGradientBoostingRegressor(**KW).fit(XX[trn], y[trn]).predict(XX[held])
            # Кластерный протокол на ТЕХ ЖЕ отложенных соединениях.
            p_clu = np.empty(len(held))
            for f in range(5):
                sel = fold[held] == f
                if sel.sum() == 0:
                    continue
                t2 = lab & (fold != f)
                p_clu[sel] = HistGradientBoostingRegressor(**KW).fit(
                    XX[t2], y[t2]).predict(XX[held[sel]])
            for proto, p in (("якорный", p_anch), ("кластерный", p_clu)):
                acc.setdefault((nm, proto), []).append(
                    (float(strae(y[held], p, y_true_upper=hi[held], y_true_lower=lo[held])),
                     float(spearmanr(y[held], p).statistic)))
        print(f"  поворот {rot}: отложено {len(held)}, якорей {len(anchors)}, "
              f"медиана сходства отложенного с обучением {np.median(simn):.3f}", flush=True)

    print(f"\n{'набор':16s} {'протокол':>12s} {'ST-RAE':>9s} {'ранг':>9s}")
    for nm in ARMS:
        for proto in ("кластерный", "якорный"):
            v = np.array(acc[(nm, proto)])
            print(f"{nm:16s} {proto:>12s} {v[:, 0].mean():9.4f} {v[:, 1].mean():9.4f}")
    print(f"\n{'выигрыш мех. блока':22s} {'кластерный':>12s} {'якорный':>10s}")
    for j, tag in ((0, "ST-RAE"), (1, "ранг")):
        a = np.array(acc[("FP+DESC+MECH", "кластерный")])[:, j].mean() - \
            np.array(acc[("FP+DESC", "кластерный")])[:, j].mean()
        b = np.array(acc[("FP+DESC+MECH", "якорный")])[:, j].mean() - \
            np.array(acc[("FP+DESC", "якорный")])[:, j].mean()
        print(f"{tag:22s} {a:12.4f} {b:10.4f}")

    print("""
Как читать. Строка сходства в поворотах --- проверка самой конструкции: если медиана сходства
отложенного с обучением не поднялась примерно к 0.587, якорь не работает и сравнивать нечего.

Дальше решает разница выигрыша между протоколами. Одинаковый выигрыш --- измерение не зависело
от режима, и число в журнале стоит. Разный --- значит все наши числа получены в зеркальном
отражении того режима, в котором их будут считать, и это первый случай, когда это ВИДНО, а не
выводится.""")


if __name__ == "__main__":
    main()
