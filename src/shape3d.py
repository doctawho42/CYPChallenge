"""Three-dimensional shape descriptors, the one axis the feature set has never had.

Why this gap is specific rather than general. Every feature in this repository is
two-dimensional: Morgan counts, RDKit's 2D descriptors, SMARTS counts. The mechanistic block
is entirely about basicity, which is to say entirely about CYP2D6. Shape is the other half of
isoform selectivity and nothing encodes it: CYP1A2 is a narrow planar slot, CYP3A4 an enormous
flexible cavity, and no column in the matrix distinguishes a flat polyaromatic from a globular
one of the same weight and logP.

One conformer is enough for what is being asked. These are whole-molecule shape indices - plane
of best fit, the normalised principal moments that place a molecule on the rod-disc-sphere
triangle, asphericity, radius of gyration - and they are stable across reasonable conformers in
a way a docking score is not. ETKDG plus a short MMFF minimisation, one conformer per molecule.

The prediction, per enzyme, and it is falsifiable in the useful way because it differs by
enzyme: planarity should matter on CYP1A2, volume and flexibility on CYP3A4, and neither much
on CYP2D6, where charge decides. Scored the way item 80 requires - the change in rank
correlation, not the raw metric.

One extra column that is not a shape index and is here for another reason. The mechanistic block
measures the distance from a basic nitrogen to an aromatic ring in BONDS, and the CYP2D6
pharmacophore is defined in ANGSTROMS, 5 to 7. A conformer gives the real distance instead of a
proxy that misreports on rings and branches, and item 82 showed that geometry group is what
carries CYP2D6.

Writes data/shape3d.npz with `train` and `test` and data/shape_names.csv.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D

import argparse
import time

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, Descriptors3D, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")

NAMES = ["PBF", "NPR1", "NPR2", "Asphericity", "Eccentricity", "InertialShapeFactor",
         "RadiusOfGyration", "SpherocityIndex", "PMI1", "PMI2", "PMI3",
         "vol_vdw", "sasa", "n_rot_frac", "bN_arom_ang_min", "bN_arom_ang_mean"]
BASIC_N = Chem.MolFromSmarts("[NX3;H2,H1,H0;!$(N=*);!$(N-[!#6;!#1]);!$(N#*);!$([N-]);"
                             "!$(N[a]);!$(N*=[O,S,N])]")


def one(smi):
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None
    mh = Chem.AddHs(m)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 0xC0FFEE
    if AllChem.EmbedMolecule(mh, ps) != 0:
        return None
    try:
        AllChem.MMFFOptimizeMolecule(mh, maxIters=400)
    except Exception:
        pass
    try:
        v = [rdMolDescriptors.CalcPBF(mh), Descriptors3D.NPR1(mh), Descriptors3D.NPR2(mh),
             Descriptors3D.Asphericity(mh), Descriptors3D.Eccentricity(mh),
             Descriptors3D.InertialShapeFactor(mh), Descriptors3D.RadiusOfGyration(mh),
             Descriptors3D.SpherocityIndex(mh), Descriptors3D.PMI1(mh),
             Descriptors3D.PMI2(mh), Descriptors3D.PMI3(mh),
             AllChem.ComputeMolVolume(mh), rdMolDescriptors.CalcTPSA(m)]
    except Exception:
        return None
    nb = m.GetNumBonds()
    v.append(rdMolDescriptors.CalcNumRotatableBonds(m) / max(nb, 1))

    # Настоящее расстояние от основного азота до центра ароматического кольца, в ангстремах.
    conf = mh.GetConformer()
    ns = [a[0] for a in mh.GetSubstructMatches(BASIC_N)]
    rings = [r for r in mh.GetRingInfo().AtomRings()
             if all(mh.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    ds = []
    for n in ns:
        pn = np.array(conf.GetAtomPosition(n))
        for r in rings:
            cen = np.mean([np.array(conf.GetAtomPosition(i)) for i in r], axis=0)
            ds.append(float(np.linalg.norm(pn - cen)))
    v += [min(ds) if ds else 0.0, float(np.mean(ds)) if ds else 0.0]
    return v


def build(smiles, tag):
    out, bad = [], 0
    t0 = time.time()
    for i, s in enumerate(smiles):
        v = one(s)
        if v is None:
            v = [0.0] * len(NAMES)
            bad += 1
        out.append(v)
        if (i + 1) % 500 == 0:
            print(f"  {tag}: {i+1}/{len(smiles)}, {time.time()-t0:.0f} с, сбоев {bad}", flush=True)
    print(f"  {tag}: готово {len(out)}, конформер не построен у {bad}", flush=True)
    return np.array(out, np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=D + "shape3d.npz")
    a = ap.parse_args()
    rows = pd.read_csv(D + "rows.csv")
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    print(f"{len(NAMES)} признаков формы на {len(rows)} + {len(te)} молекул", flush=True)
    tr_v = build(list(rows.SMILES), "обучение")
    te_v = build(list(te.SMILES), "тест")
    np.savez_compressed(a.out, train=tr_v, test=te_v)
    pd.Series(NAMES).to_csv(D + "shape_names.csv", index=False, header=False)
    print(f"\nсохранено: {a.out}")
    fin = np.isfinite(tr_v).all(0)
    print(f"столбцов с нечисловым: {int((~fin).sum())}; "
          f"мёртвых (нулевой разброс): {int((tr_v.std(0) < 1e-9).sum())}")


if __name__ == "__main__":
    main()
