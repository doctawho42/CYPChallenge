"""Is the greedy split criterion myopic about the coordinating mode, and is there a mode at all.

The argument this gates, stated so it can fail. Item 156 decoded the fingerprint and found that
every fragment surviving a control for size, lipophilicity and aromaticity is **sp2 ring nitrogen
with an available lone pair** -- pyrimidine C2, N-aryl azole, pyridine nitrogen, azine carbon. That
is one mechanism, type-II coordination to the haem iron. But the model carries it as a
**disjunction**: `n_pyridineN`, `n_imidazoleN` and the fingerprint bits are separate columns, so
"has a coordinating nitrogen of any kind" costs an axis-aligned tree up to four splits, each with a
partial gain, each competing against 2295 columns. One mode indicator is one split.

That is a representational cost, and it is not the strong half of the argument. The strong half is
that a greedy criterion scores a split by its **immediate** gain, and the mode's immediate gain may
be small while the gain from splitting *within* it is large -- because the structure-activity
relationship differs between modes. Coordinators are ranked by lone-pair availability; lipophiles
by bulk. Greedy is short-sighted exactly there.

**The chemistry is drawn where it is sharp, not where it is convenient.** `[nX2]` is a pyridine-type
aromatic nitrogen: two connections, lone pair in the ring plane, available to iron. `[nX3]` is
pyrrole-type: the lone pair is in the pi system and does not coordinate. So the mode is "has at
least one `[nX2]`", which is a mechanism, not a count.

The gate, both halves pre-registered:

    (а) стандартизованная разница средних pIC50 между модами  ->  должна быть МАЛА,
        иначе жадный критерий этот сплит и так бы взял, и близорукости нет;
    (б) обучить ВНУТРИ моды A, предсказать моду B, против B->B  ->  разрыв должен быть ВЕЛИК,
        это и есть расхождение SAR, которого жадность не видит.

    мало и велико   -> близорукость доказана, принудительный корневой сплит оправдан ИЗМЕРЕНИЕМ
    оба малы        -> моды не различаются, закрыто
    оба велики      -> жадность справляется сама, закрыто

The transfer half is measured on the **same held-out rows** in both directions, so the comparison is
not confounded by which compounds are being predicted: for each fold of mode B, train once on all of
mode A and once on mode B minus that fold, and predict that fold with both.

Reads data/feats.npz and data/rows.csv. Writes nothing. ~10 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, tutorial
tutorial()

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.tree import DecisionTreeRegressor
from rdkit import Chem, RDLogger

from cypsplit import butina_folds

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
NTREE, LR, DEPTH, MF = 150, 0.06, 5, 0.3
COORD = Chem.MolFromSmarts("[nX2]")          # пиридинового типа: пара в плоскости, доступна
PYRR = Chem.MolFromSmarts("[nX3;H1]")        # пиррольного типа: пара в pi-системе, НЕ координирует


def boost(Xtr, ytr, Xte, seed):
    rng = np.random.default_rng(seed)
    s = np.full(len(ytr), float(ytr.mean())); p = np.full(len(Xte), s[0])
    for _ in range(NTREE):
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MF,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(Xtr, ytr - s)
        s = s + LR * t.predict(Xtr); p = p + LR * t.predict(Xte)
    return p


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)
    fold, _ = butina_folds(list(rows.SMILES))

    coord, pyrr = [], []
    for s in rows.SMILES:
        m = Chem.MolFromSmiles(s)
        coord.append(bool(m and m.HasSubstructMatch(COORD)))
        pyrr.append(bool(m and m.HasSubstructMatch(PYRR)))
    coord = np.array(coord); pyrr = np.array(pyrr)
    print(f"молекул {len(rows)}; с координирующим [nX2] {coord.sum()} ({100*coord.mean():.1f} %); "
          f"с пиррольным [nX3;H1] {pyrr.sum()}; только пиррольный, без [nX2]: "
          f"{int((pyrr & ~coord).sum())}\n")

    print("=== (а) различаются ли моды ПО УРОВНЮ --- то, что жадность видит сразу ===")
    print(f"{'фермент':8s} {'n коорд':>8s} {'n проч':>7s} {'сред коорд':>11s} {'сред проч':>10s} "
          f"{'станд. разница':>15s}")
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr[col].to_numpy(float)
        a, b = m & coord, m & ~coord
        sd = np.sqrt(((a.sum() - 1) * y[a].var(ddof=1) + (b.sum() - 1) * y[b].var(ddof=1))
                     / max(a.sum() + b.sum() - 2, 1))
        print(f"{c:8s} {a.sum():8d} {b.sum():7d} {y[a].mean():11.3f} {y[b].mean():10.3f} "
              f"{(y[a].mean()-y[b].mean())/max(sd,1e-9):+15.3f}")

    print("\n=== (б) переносится ли SAR между модами --- то, чего жадность не видит ===")
    print("  на ОДНИХ И ТЕХ ЖЕ отложенных строках, и с УРАВНЕННЫМ размером обучения:")
    print("  моды относятся как 71.5 к 28.5, поэтому без уравнивания сравнение идёт про")
    print("  число строк, а не про перенос SAR. Обе обучающие выборки урезаются до общего")
    print("  минимума, и результат усредняется по трём независимым урезаниям.")
    print(f"\n{'фермент':8s} {'мода':>9s} {'n':>5s} {'своя мода':>10s} {'чужая мода':>11s} "
          f"{'разрыв':>8s}")
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr[col].to_numpy(float)
        for nm, sel, other in (("коорд", m & coord, m & ~coord),
                               ("прочие", m & ~coord, m & coord)):
            if sel.sum() < 120 or other.sum() < 120:
                print(f"{c:8s} {nm:>9s} {sel.sum():5d}   мало строк")
                continue
            idx = np.where(sel)[0]
            ros, rcs, sizes = [], [], []
            for rep in range(3):
                g = np.random.default_rng(rep)
                own = np.zeros(len(idx)); cross = np.zeros(len(idx))
                for f in range(5):
                    te = fold[idx] == f
                    if te.sum() < 10:
                        continue
                    to = idx[fold[idx] != f]
                    tc = np.where(other & (fold != f))[0]
                    k = min(len(to), len(tc))            # уравнивание размера
                    to = g.choice(to, k, replace=False)
                    tc = g.choice(tc, k, replace=False)
                    Xte = X[idx[te]]
                    own[te] = boost(X[to], y[to], Xte, 1000 + f * 10 + rep)
                    cross[te] = boost(X[tc], y[tc], Xte, 2000 + f * 10 + rep)
                    sizes.append(k)
                ros.append(spearmanr(y[idx], own).statistic)
                rcs.append(spearmanr(y[idx], cross).statistic)
            ro, rc = float(np.mean(ros)), float(np.mean(rcs))
            print(f"{c:8s} {nm:>9s} {len(idx):5d} {ro:10.3f} {rc:11.3f} {ro-rc:+8.3f}"
                  f"   (обучение по {int(np.mean(sizes))} строк, sd разрыва "
                  f"{np.std(np.array(ros)-np.array(rcs), ddof=1):.3f})")

    print("""
Как читать. Предрегистрация из письма, дословно.

Столбец «станд. разница» МАЛ и столбец «разрыв» ВЕЛИК --- близорукость жадного критерия
доказана: моды почти не отличаются по уровню, поэтому сплит по моде не даёт немедленного
выигрыша и жадность его не берёт, но SAR внутри мод разные, поэтому выигрыш есть на шаг
дальше. Принудительный корневой сплит оправдан измерением.

Оба малы --- мод нет, закрыто. Оба велики --- жадность и так возьмёт этот сплит, закрыто.

Разрыв надо читать против поферментного пола пункта 165 (0.0061 / 0.0071 / 0.0049 / 0.0033),
но здесь он заведомо больше: перенос между популяциями --- эффект другого порядка, чем
добавление признака.""")


if __name__ == "__main__":
    main()
