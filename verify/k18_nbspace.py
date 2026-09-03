"""In which space is a neighbour informative? Measured with no model at all.

Why this comes before the methods that use neighbours. Two things in the queue - neighbour
labels as features, and a Gaussian process with a similarity kernel - both need a definition of
"close", and everything in this repository has silently used one: Tanimoto over Morgan counts.
That choice has never been tested, and it may weigh more than the model built on top of it.

The measurement needs no training. For every compound, take the label of its nearest neighbour
among the OTHER folds and call that the prediction. Rank-correlate with the truth. Whichever
space wins is the space in which a neighbour actually carries information about activity.

Four spaces, chosen because they encode different things:

  Morgan / Tanimoto        substructure. What the molecule is made of;
  pharmacophore pairs      the arrangement of donors, acceptors, charges and rings, without
                           regard to scaffold. Two molecules with different skeletons and the
                           same basic-nitrogen-to-aromatic geometry are close here and far in
                           Morgan space, which for a binding site is the relevant sense;
  chemprop embedding       what a pretrained D-MPNN considers close. Item 68 showed it carries
                           something our features do not;
  descriptors, cosine      bulk properties - size, lipophilicity, polarity.

The prediction worth recording: for a binding site, pharmacophore similarity should beat
substructure similarity, most clearly on CYP2D6 where the pharmacophore is explicit and known.
If Morgan wins everywhere the current choice is vindicated and there is nothing to change.

A one-nearest-neighbour rule is a deliberately weak predictor and its absolute rank will be far
below the boosting's. The comparison here is between spaces, not against the model.

Reads data/feats.npz, data/emb_chemprop_rdkit2d.npz, data/rows.csv. Prints only.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D

import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Pharm2D import Generate, Gobbi_Pharm2D

from cypsplit import butina_folds

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]


def nn_predict(sim_fn, y, fold):
    """Метка ближайшего соседа из ДРУГИХ фолдов, для каждого соединения."""
    p = np.zeros(len(y))
    for f in np.unique(fold):
        te = np.where(fold == f)[0]
        trn = np.where(fold != f)[0]
        if len(te) == 0 or len(trn) == 0:
            continue
        S = sim_fn(te, trn)
        p[te] = y[trn[np.argmax(S, axis=1)]]
    return p


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    smiles = list(rows.SMILES)
    mols = [Chem.MolFromSmiles(s) for s in smiles]

    print("строю представления", flush=True)
    t0 = time.time()
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    FPm = [gen.GetFingerprint(m) for m in mols]
    print(f"  Морган {time.time()-t0:.0f} с", flush=True)
    t0 = time.time()
    FPp = [Generate.Gen2DFingerprint(m, Gobbi_Pharm2D.factory) for m in mols]
    print(f"  фармакофорные пары {time.time()-t0:.0f} с", flush=True)

    E = np.load(D + "emb_chemprop_rdkit2d.npz")["train"].astype(np.float64)
    E = E / (np.linalg.norm(E, axis=1, keepdims=True) + 1e-12)
    Dsc = z["DESC"].astype(np.float64)
    Dsc = np.nan_to_num(Dsc, posinf=0.0, neginf=0.0)
    Dsc = (Dsc - Dsc.mean(0)) / (Dsc.std(0) + 1e-9)
    Dsc = Dsc / (np.linalg.norm(Dsc, axis=1, keepdims=True) + 1e-12)

    def bulk(fps):
        return lambda te, trn: np.array([DataStructs.BulkTanimotoSimilarity(
            fps[i], [fps[j] for j in trn]) for i in te])

    SPACES = {
        "Морган": lambda idx: bulk([FPm[i] for i in idx]),
        "фармакофор": lambda idx: bulk([FPp[i] for i in idx]),
        "эмбеддинг": lambda idx: (lambda te, trn: E[idx][te] @ E[idx][trn].T),
        "дескрипторы": lambda idx: (lambda te, trn: Dsc[idx][te] @ Dsc[idx][trn].T),
    }

    print(f"\n{'пространство':14s} " + " ".join(f"{c:>9s}" for c in CYPS) + f" {'среднее':>9s}")
    for name, mk in SPACES.items():
        vals = []
        for c in CYPS:
            col = f"{c}_pIC50_direct_inhibition"
            m = tr[col].notna().to_numpy()
            idx = np.where(m)[0]
            y = tr.loc[m, col].to_numpy()
            fold, _ = butina_folds(smiles, seed=0)
            fi = fold[m]
            p = nn_predict(mk(idx), y, fi)
            vals.append(float(spearmanr(y, p).statistic))
        print(f"{name:14s} " + " ".join(f"{v:9.4f}" for v in vals)
              + f" {np.mean(vals):9.4f}", flush=True)

    print("""
Как читать. Числа заведомо ниже бустинга --- один сосед это нарочно слабый предсказатель.
Сравниваются ПРОСТРАНСТВА между собой, а не с моделью.

Если побеждает фармакофор, то и признаки по соседям, и ядро гауссова процесса надо строить в
нём, а не в морганском, --- и это меняет два метода из очереди до того, как они написаны. Если
побеждает Морган, текущий выбор оправдан и менять нечего.""")


if __name__ == "__main__":
    main()
