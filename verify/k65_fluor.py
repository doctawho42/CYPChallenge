"""Is the fluorescence artefact visible in the LABELS? The one confounder with a built-in control.

Section 7 of the document carries a "sample artefact channel" and marks it as never measured, and
the journal confirms: no item mentions fluorescence, chromophores or quenching. It is the strongest
unrun idea in the document, because unlike everything else here it comes with a control the design
supplies for free.

    1A2, 2C9, 3A4   флуоресцентное считывание   -> собственная флуоресценция и тушение
                                                   дают ЛОЖНОЕ «ингибирование»
    2D6             масс-спектрометрия          -> артефакта нет ПО ПОСТРОЕНИЮ

**This file does not build the channel. It asks whether there is anything to correct**, which is a
day cheaper and decides the question either way. If a chromophore score predicts the label on the
three fluorescent enzymes and not on the mass-spectrometry one, the artefact is in the benchmark's
labels -- a statement about the data rather than about our model, and one the organisers would want.
If it predicts everywhere, it is a shape descriptor and the channel closes; if nowhere, likewise.

**The direction is pre-registered and it is not symmetric.** Autofluorescence and quenching both
make the read look like less product, that is like MORE inhibition, so a chromophore should appear
FALSELY POTENT. The coefficient must be **positive** on 1A2, 2C9, 3A4.

**The confounder that makes the naive test invalid, and it is specific.** CYP2D6 binds through a
salt bridge to a protonated basic nitrogen (item 81), so nitrogen-bearing aromatics are genuinely
more active there for reasons having nothing to do with the readout. Aromatic N-heterocycles are
simultaneously chromophores and bases, so a raw correlation on 2D6 confounds the two and the control
arm stops being a control. Every association below is therefore reported twice: raw, and partial
after removing molecular weight, lipophilicity **and the basic-nitrogen count**.

Item 156 is why this could be large rather than another hundredth. It decoded the fifty fingerprint
bits carrying most of the model's signal and found them to be ring nitrogen able to coordinate the
haem iron -- pyrimidine, N-aryl-azole, pyridine nitrogen, azine carbon, methoxyaryl. Those are also
chromophore classes. Part of what the model reads as a pharmacophore on three enzymes may be a
signature of readout interference. Which of the two it is, is exactly what the 2D6 arm decides.

Reads data/rows.csv and the training labels. RDKit only. Minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D

import argparse

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
FLUOR = {"CYP1A2": True, "CYP2C9": True, "CYP3A4": True, "CYP2D6": False}

# Сильные хромофоры и тушители, по которым спрашивают отдельно.
PATT = {
    "нитро":   "[N+](=O)[O-]",
    "азо":     "[#6]N=N[#6]",
    "хинон":   "O=C1C=CC(=O)C=C1",
    "стильбен": "[#6]=[#6]c1ccccc1",
}


def chromo(smi):
    """Признаки хромофорности. Физически определяющая величина --- размер сопряжённой системы:
    именно она задаёт длину волны поглощения, а совпадение с длиной волны пробы и есть помеха."""
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    ri = m.GetRingInfo()
    arom_atoms = sum(1 for a in m.GetAtoms() if a.GetIsAromatic())
    conj = sum(1 for b in m.GetBonds() if b.GetIsConjugated())
    # крупнейшая слитая ароматическая система
    ars = [set(r) for r in ri.AtomRings()
           if all(m.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    merged = []
    for r in ars:
        hit = [g for g in merged if g & r]
        for g in hit:
            merged.remove(g)
            r = r | g
        merged.append(r)
    biggest = max((len(g) for g in merged), default=0)
    d = {"аром_атомов": arom_atoms, "сопряж_связей": conj,
         "крупнейшая_система": biggest, "аром_колец": rdMolDescriptors.CalcNumAromaticRings(m)}
    for k, sm in PATT.items():
        q = Chem.MolFromSmarts(sm)
        d[k] = int(m.HasSubstructMatch(q)) if q is not None else 0
    d["MolWt"] = Descriptors.MolWt(m)
    d["logP"] = Descriptors.MolLogP(m)
    # основные азоты: то, чем живёт CYP2D6, и главный конфаундер этого теста
    d["основн_N"] = sum(1 for a in m.GetAtoms() if a.GetSymbol() == "N"
                        and not a.GetIsAromatic() and a.GetTotalNumHs() >= 0
                        and a.GetFormalCharge() == 0 and a.GetDegree() <= 3)
    return d


def partial(x, y, Z):
    """Спирменова частная связь x с y после снятия столбцов Z (регрессия по рангам)."""
    R = lambda v: pd.Series(v).rank().to_numpy()
    Zr = np.column_stack([R(Z[:, j]) for j in range(Z.shape[1])] + [np.ones(len(y))])
    rx = R(x) - Zr @ np.linalg.lstsq(Zr, R(x), rcond=None)[0]
    ry = R(y) - Zr @ np.linalg.lstsq(Zr, R(y), rcond=None)[0]
    return float(spearmanr(rx, ry).statistic)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    F = pd.DataFrame([chromo(s) for s in rows.SMILES])
    print(f"признаки хромофорности посчитаны для {F.notna().all(1).sum()} из {len(rows)}\n")
    for k in PATT:
        print(f"  {k:10s} {int(F[k].sum()):5d} молекул")
    print()

    keys = ["крупнейшая_система", "сопряж_связей", "аром_колец"]
    ctrl = F[["MolWt", "logP", "основн_N"]].to_numpy(float)
    rng = np.random.default_rng(a.seed)

    for k in keys:
        print(f"=== {k} ===")
        print(f"{'фермент':8s} {'считыв.':>9s} {'сырая rho':>10s} {'частная':>9s} "
              f"{'переставл.':>11s}")
        for c in CYPS:
            col = f"{c}_pIC50_direct_inhibition"
            m = tr[col].notna().to_numpy()
            y = tr.loc[m, col].to_numpy()
            x = F.loc[m, k].to_numpy(float)
            raw = float(spearmanr(x, y).statistic)
            par = partial(x, y, ctrl[m])
            perm = float(np.mean([partial(rng.permutation(x), y, ctrl[m]) for _ in range(5)]))
            print(f"{c:8s} {('флуор' if FLUOR[c] else 'МАСС-СП'):>9s} "
                  f"{raw:10.3f} {par:9.3f} {perm:11.3f}")
        print()

    print("""Как читать. Предрегистрация: артефакт даёт ЛОЖНУЮ потентность, поэтому связь должна
быть ПОЛОЖИТЕЛЬНОЙ на 1A2, 2C9, 3A4 и около нуля на 2D6, который читается масс-спектрометрией.

Решает ЧАСТНАЯ связь, а не сырая. Сырая смешивает хромофорность с размером и липофильностью, а
на 2D6 ещё и с основным азотом --- у этого фермента солевой мостик, и азотистые ароматические
соединения там активнее по настоящей химии. Без снятия основного азота контрольная рука
перестаёт быть контролем, и именно это делает наивную версию теста недействительной.

Равномерная связь по всем четырём --- это дескриптор формы, а не артефакт считывания, и канал
закрывается. Связь на трёх флуоресцентных при нуле на 2D6 --- артефакт в МЕТКАХ бенчмарка, и это
утверждение о данных, а не о нашей модели.""")


if __name__ == "__main__":
    main()
