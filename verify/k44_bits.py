"""Decode the fifty fingerprint bits the boosting actually uses into substructures.

Where this comes from. Item 96 measured the whole 2048-bit block at +0.032 of rank. Item 138
measured that the trees spend 30 to 38 per cent of their splits on 89 per cent of the columns,
so the block is consulted rarely for how much it matters. Item 150 then measured that **fifty
bits carry 80.7 per cent of the block's effect and twenty carry 69.8**. The block is not working
in bulk; it is a short list.

A short list of what, though. A bit index is not a substructure, and this file's rule is that a
number is not a finding until something is counted. Two things have to be counted before any
chemistry may be claimed:

  **коллизии**   Morgan environments are hashed into 2048 buckets, so one bit can collect several
                 chemically unrelated fragments. A bit that is 95 per cent one environment is a
                 substructure with a name. A bit that is five environments at a fifth each is a
                 bucket, and calling it an alert would be inventing chemistry out of a hash
                 collision. Every bit below carries the count of distinct environments mapped to
                 it and the share of the most common one.

  **направление**  a bit the model consults often may still be neutral in the label. For each bit
                 the difference in mean pIC50 between the molecules that carry it and those that
                 do not, per enzyme, with the carrier count. That is what turns "the trees look
                 here" into "this fragment raises or lowers measured inhibition".

The ranking uses the same booster and the same constants as `verify/k40_topk.py`, so the fifty
bits here are the fifty bits measured there. It is fitted on the whole labelled set for each
enzyme rather than per fold: fold-wise ranking exists in k40 because it *scores* the reduced
matrix and selecting on the answer would inflate it, while this script only describes bits and
scores nothing. The two rankings are compared so the difference is visible rather than assumed.

Environments are rendered with MolFragmentToSmiles rooted at the central atom, which is what a
Morgan environment is -- an atom plus everything within the radius -- so the root is chemically
meaningful and not a rendering choice.

Reads data/feats.npz and data/rows.csv. Writes results/bits50.csv. ~3 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES

import argparse
from collections import Counter, defaultdict

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from sklearn.tree import DecisionTreeRegressor

from cypsplit import butina_folds

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
NTREE, LR, DEPTH, MF = 150, 0.06, 5, 0.3      # пины из verify/k40_topk.py
NFP = 2048
TOPK = 50


def bit_counts(X, y, nfp, seed):
    """Сколько раз каждый бит выбран сплитом. Тот же бустинг, что в k40."""
    rng = np.random.default_rng(seed)
    s = np.full(len(y), float(y.mean()))
    cnt = np.zeros(nfp, int)
    for _ in range(NTREE):
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MF,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(X, y - s)
        s = s + LR * t.predict(X)
        f = t.tree_.feature
        f = f[(f >= 0) & (f < nfp)]
        np.add.at(cnt, f, 1)
    return cnt


def env_smiles(mol, atom, rad):
    """SMILES окружения Моргана: центральный атом плюс всё в пределах радиуса."""
    if rad == 0:
        return Chem.MolFragmentToSmiles(mol, atomsToUse=[atom], canonical=True)
    env = Chem.FindAtomEnvironmentOfRadiusN(mol, rad, atom)
    if not env:
        return None
    atoms = {atom}
    for b in env:
        bd = mol.GetBondWithIdx(b)
        atoms.add(bd.GetBeginAtomIdx())
        atoms.add(bd.GetEndAtomIdx())
    try:
        return Chem.MolFragmentToSmiles(mol, atomsToUse=sorted(atoms), bondsToUse=list(env),
                                        rootedAtAtom=atom, canonical=True)
    except Exception:
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--top", type=int, default=TOPK)
    ap.add_argument("--out", default=RES + "bits50.csv")
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    FP = z["FP"].astype(np.float32)
    DM = np.hstack([z["DESC"], z["MECH"]]).astype(np.float32)
    fold, _ = butina_folds(list(rows.SMILES))

    # --- 1. Ранжирование битов, на всей помеченной выборке и по фолдам ---
    cnt_all, cnt_fold = {}, {}
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy()
        Xi = np.hstack([FP[m], DM[m]])
        cnt_all[c] = bit_counts(Xi, y, NFP, 0)
        cf = np.zeros(NFP, int)
        for f in range(5):
            trn = fold[m] != f
            cf += bit_counts(Xi[trn], y[trn], NFP, 10 + f)
        cnt_fold[c] = cf
        print(f"{c}: сплитов по FP {cnt_all[c].sum()}, ненулевых битов "
              f"{int((cnt_all[c] > 0).sum())}", flush=True)

    total = sum(cnt_all[c] for c in CYPS)
    top = np.argsort(-total)[:a.top]
    tf = np.argsort(-sum(cnt_fold[c] for c in CYPS))[:a.top]
    print(f"\nустойчивость ранжирования: из {a.top} битов общими с пофолдовым "
          f"ранжированием {len(set(top.tolist()) & set(tf.tolist()))}")

    # Насколько списки ферментов различны
    tops = {c: set(np.argsort(-cnt_all[c])[:a.top].tolist()) for c in CYPS}
    print(f"\n{'':10s}" + "".join(f"{c[3:]:>8s}" for c in CYPS))
    for c in CYPS:
        print(f"{c:10s}" + "".join(f"{len(tops[c] & tops[d]):8d}" for d in CYPS))
    inter = set.intersection(*tops.values())
    print(f"общих у всех четырёх: {len(inter)} из {a.top}")

    # --- 2. Расшифровка: какие окружения дают эти биты ---
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=NFP)
    want = set(int(b) for b in top)
    envs = defaultdict(Counter)      # бит -> Counter(SMILES окружения)
    example = {}                     # бит -> SMILES молекулы-примера
    print("\nразбираю окружения...", flush=True)
    for i, smi in enumerate(rows.SMILES):
        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            continue
        ao = rdFingerprintGenerator.AdditionalOutput(); ao.AllocateBitInfoMap()
        gen.GetCountFingerprint(mol, additionalOutput=ao)
        for b, occ in ao.GetBitInfoMap().items():
            if b not in want:
                continue
            example.setdefault(b, smi)
            for atom, rad in occ:
                e = env_smiles(mol, atom, rad)
                if e:
                    envs[b][(e, rad)] += 1

    # --- 3. Направление в метке ---
    recs = []
    for b in top:
        b = int(b)
        car = FP[:, b] > 0
        cc = envs[b]
        tot = sum(cc.values())
        (dom, rad), n1 = cc.most_common(1)[0] if cc else (("?", -1), 0)
        r = {"бит": b, "сплитов": int(total[b]), "молекул": int(car.sum()),
             "разных окружений": len(cc), "доля главного": round(n1 / tot, 3) if tot else 0.0,
             "радиус": rad, "подструктура": dom, "пример": example.get(b, "")}
        for c in CYPS:
            col = f"{c}_pIC50_direct_inhibition"
            m = tr[col].notna().to_numpy()
            y = tr[col].to_numpy(float)
            g, h = m & car, m & ~car
            r[f"n_{c[3:]}"] = int(g.sum())
            r[f"d_{c[3:]}"] = (round(float(y[g].mean() - y[h].mean()), 3)
                               if g.sum() >= 10 and h.sum() >= 10 else np.nan)
        recs.append(r)
    df = pd.DataFrame(recs)
    df.to_csv(a.out, index=False)

    clean = df[df["доля главного"] >= 0.8]
    print(f"\nиз {len(df)} битов НЕ смешанных (главное окружение >= 80 %): {len(clean)}")
    print(f"медиана доли главного окружения: {df['доля главного'].median():.2f}")
    print(f"медиана числа разных окружений на бит: {df['разных окружений'].median():.0f}\n")

    cols = ["бит", "сплитов", "молекул", "доля главного", "радиус", "подструктура"] + \
           [f"d_{c[3:]}" for c in CYPS]
    print("ВСЕ ОТОБРАННЫЕ БИТЫ, по числу сплитов")
    print(df[cols].to_string(index=False))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. «доля главного» --- первое, на что смотреть. Бит с долей ниже 0.8 --- это
коллизия хеша, а не подструктура, и называть его алертом значит выдумывать химию из
совпадения адресов. Такие строки надо выбросить ДО любого химического чтения.

«d_XXX» --- разность средней pIC50 у носителей фрагмента и остальных, по ферментам.
Положительное значит «носители ингибируют сильнее». NaN --- носителей меньше десяти,
разность не считается.

Радиус говорит, что это: 0 --- одиночный атом с окружением по свойствам, 1 --- атом с
соседями, 2 --- полноценный фрагмент. Биты радиуса 0 несут не подструктуру, а тип атома, и
их химическое прочтение другое.""")


if __name__ == "__main__":
    main()
