"""Before building the screening head: is the screening readout informative, and about whom.

Item 5 of the outside reading proposes a second head on the shared trunk, predicting the
single-point screen at 49.5 uM alongside the curve-derived pIC50. The argument is instrumental
-- under a Hill model the screen is a monotone squashed version of the curve, so transfer is
forced -- and the payoff is claimed to be 11509 extra labels, most of them on compounds the
curve set does not contain, which would correct the organisers' selection rather than merely
add rows.

Every part of that is checkable in an hour, and the same discipline closed items 1, 6 and 7 in
less time than building them would have taken. Three questions, in the order that can kill it
fastest:

  1. Does the screen rank compounds the way the curve does, on the compounds that have both?
     If the correlation is weak, the instrumental argument is wrong about this instrument and
     nothing downstream matters.

  2. How many compounds have a screen and no curve? That is the actual extra information; a
     screen on a compound we already have a curve for adds a noisy copy of a label we hold.

  3. Are those compounds different from the ones with curves? The proposal's real claim is
     that they correct a selection, which requires them to sit somewhere else -- in potency,
     in structure, or both. If they are a random subset, the head is just more rows.

Reads data/cyp-challenge-single-concentration-TRAIN.csv and the inhibition table. ~3 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]


def main():
    tr = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    print(f"скрининг: {len(sc)} строк, {sc.Molecule_Name.nunique()} молекул, "
          f"{sc.enzyme.nunique()} ферментов")
    print(f"кривые:   {tr.Molecule_Name.nunique()} молекул\n")

    print("1. Ранжирует ли скрининг так же, как кривая (на пересечении):")
    print(f"{'фермент':8s} {'общих':>7s} {'rho(log2fc, pIC50)':>20s} {'|rho|':>8s}")
    for c in CYPS:
        s = sc[sc.enzyme == c][["Molecule_Name", "log2fc_estimate"]]
        s = s.groupby("Molecule_Name").log2fc_estimate.mean()
        col = f"{c}_pIC50_direct_inhibition"
        t = tr.loc[tr[col].notna(), ["Molecule_Name", col]].set_index("Molecule_Name")[col]
        j = pd.concat([s, t], axis=1, join="inner").dropna()
        r = spearmanr(j.log2fc_estimate, j[col]).statistic if len(j) > 10 else float("nan")
        print(f"{c:8s} {len(j):7d} {r:20.3f} {abs(r):8.3f}")

    print("\n2. Сколько скрининга приходится на молекулы БЕЗ кривой:")
    print(f"{'фермент':8s} {'всего скрин.':>13s} {'есть кривая':>12s} {'нет кривой':>11s} {'доля новых':>12s}")
    new = {}
    for c in CYPS:
        s = set(sc[sc.enzyme == c].Molecule_Name)
        have = set(tr.loc[tr[f"{c}_pIC50_direct_inhibition"].notna(), "Molecule_Name"])
        n = s - have
        new[c] = n
        print(f"{c:8s} {len(s):13d} {len(s & have):12d} {len(n):11d} {100*len(n)/max(len(s),1):11.1f} %")
    allnew = set().union(*new.values())
    print(f"  молекул без единой кривой, но со скринингом: {len(allnew)}")

    print("\n3. Отличаются ли они по показанию скрининга от тех, у кого кривая есть:")
    print(f"{'фермент':8s} {'медиана |log2fc| с кривой':>26s} {'без кривой':>12s} {'rho(есть кривая, |log2fc|)':>28s}")
    for c in CYPS:
        s = sc[sc.enzyme == c].groupby("Molecule_Name").log2fc_estimate.mean()
        have = set(tr.loc[tr[f"{c}_pIC50_direct_inhibition"].notna(), "Molecule_Name"])
        h = np.abs(s[s.index.isin(have)]); n = np.abs(s[~s.index.isin(have)])
        flag = np.concatenate([np.ones(len(h)), np.zeros(len(n))])
        val = np.concatenate([h.to_numpy(), n.to_numpy()])
        r = spearmanr(flag, val).statistic
        print(f"{c:8s} {np.median(h):26.3f} {np.median(n):12.3f} {r:28.3f}")

    print("\n4. Структурно ли они где-то ещё:")
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    have = sorted(set(tr.Molecule_Name))
    smi = dict(zip(tr.Molecule_Name, tr.SMILES))
    smi.update(dict(zip(sc.Molecule_Name, sc.SMILES)))
    rng = np.random.default_rng(0)
    ref_names = list(rng.choice(have, min(800, len(have)), replace=False))
    ref = [gen.GetFingerprint(Chem.MolFromSmiles(smi[n])) for n in ref_names
           if Chem.MolFromSmiles(smi.get(n, "")) is not None]
    for tag, names in (("с кривой", ref_names[:300]),
                       ("без кривой", list(rng.choice(sorted(allnew),
                                                      min(300, len(allnew)), replace=False)))):
        v = []
        for n in names:
            m = Chem.MolFromSmiles(smi.get(n, ""))
            if m is None:
                continue
            f = gen.GetFingerprint(m)
            sims = DataStructs.BulkTanimotoSimilarity(f, ref)
            sims = sorted(sims)[:-1]          # выбросить самого себя
            v.append(max(sims) if sims else 0.0)
        print(f"  {tag:12s} медиана сходства с набором кривых: {np.median(v):.3f}  (n={len(v)})")

    print("""
Как читать. Пункт 1 --- это выключатель: слабая корреляция означает, что инструментальный
довод к этому прибору не относится, и голова учила бы шум.

Пункт 2 отделяет объём от новизны. Пункты 3 и 4 проверяют то, ради чего предложение и стоит
дороже простого добавления строк: если молекулы без кривой лежат там же и показывают то же,
это не коррекция отбора, а просто ещё данные --- и тогда 2--3 дня работы покупают то же, что
покупает любая регуляризация.""")


if __name__ == "__main__":
    main()
