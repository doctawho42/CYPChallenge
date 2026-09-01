"""Shape and pharmacophore overlay on the four co-crystal ligands: a pair feature, not a ligand one.

Why this and not another descriptor block. Item 189 collected six measurements: physics enters
through the measurement model and returns zero through the feature matrix -- SMARTCyp, hand-built
site blocks, the quantum block, the Hill residual as a weight. Item 168 named the reason and the
exception. Every one of those is a **function of the ligand**, and the ligand block already has 2295
columns of those. A quantity computed from the pair **(ligand, cavity)** is not, and overlay onto a
cavity's own co-crystal ligand is the cheapest such quantity: hours rather than the night twenty
thousand docking runs would take.

**The bound conformation is the whole point and is not optional.** The co-crystal ligand's pose is a
cast of the cavity; a freely generated conformer of alpha-naphthoflavone is just another ligand, and
overlaying molecules on that measures ligand-to-ligand similarity, which is exactly the class that
returned zero six times. So the references are the deposited coordinates, downloaded from RCSB:

    CYP1A2   2HI4  BHF   альфа-нафтофлавон      21 атом
    CYP2C9   1R9O  FLP   флурбипрофен           18
    CYP2D6   4WNV  QI9   хинин                  24
    CYP3A4   3NXU  RIT   ритонавир              50

Two of those corrected a wrong recollection when the PDB was asked instead: 4WNV is quinine, not
thioridazine, and 3NXU is ritonavir, not ketoconazole.

**The feature is differences between cavities, not the scores.** The common part of the four scores
is how large and how greasy the molecule is -- which the ligand block already carries in full -- and
it cancels in the contrast, leaving complementarity to a particular cavity. That is also why it fits
this repository: item 132 established that contrast is the quantity pooling exploits.

Alignment is Crippen O3A rather than MMFF O3A: it types atoms by their logP and molar-refractivity
contributions, so shape and "colour" enter together, and it does not fail on molecules MMFF cannot
type. Several conformers are generated per molecule and the best alignment against each reference is
kept, because a single ETKDG conformer is an arbitrary draw and the question is whether the molecule
*can* present the cavity's shape.

Writes data/overlay.npz with the raw scores; centring and the wrong-cavity control belong to the
ablation, not here.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D

import argparse
import glob
import time

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem, rdMolAlign, rdMolDescriptors

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]


def load_refs():
    refs = {}
    for f in sorted(glob.glob(D + "cocrystal/*.sdf")):
        enz = f.split("/")[-1].split("_")[0]
        m = next(Chem.ForwardSDMolSupplier(f, removeHs=False), None)
        if m is None:
            raise SystemExit(f"не прочитан {f}")
        m = Chem.AddHs(m, addCoords=True)
        refs[enz] = m
    miss = [c for c in CYPS if c not in refs]
    if miss:
        raise SystemExit(f"нет со-кристальных лигандов для {miss}")
    return refs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nconf", type=int, default=3)
    ap.add_argument("--out", default=D + "overlay.npz")
    a = ap.parse_args()

    refs = load_refs()
    print("эталоны (связанные позы из PDB):")
    for c in CYPS:
        print(f"  {c}: {refs[c].GetNumAtoms()} атомов с водородами")
    ref_par = {c: rdMolDescriptors._CalcCrippenContribs(refs[c]) for c in CYPS}

    rows = pd.read_csv(D + "rows.csv")
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    out = {}
    for tag, smiles in (("train", list(rows.SMILES)), ("test", list(te.SMILES))):
        t0 = time.time()
        S = np.full((len(smiles), 4), np.nan, np.float32)
        bad = 0
        for i, smi in enumerate(smiles):
            m = Chem.MolFromSmiles(smi)
            if m is None:
                bad += 1
                continue
            m = Chem.AddHs(m)
            p = AllChem.ETKDGv3()
            p.randomSeed = 0xBEEF
            p.pruneRmsThresh = 0.5
            cids = AllChem.EmbedMultipleConfs(m, numConfs=a.nconf, params=p)
            if not len(cids):
                bad += 1
                continue
            try:
                AllChem.MMFFOptimizeMoleculeConfs(m, maxIters=200)
            except Exception:
                pass
            par = rdMolDescriptors._CalcCrippenContribs(m)
            for e, c in enumerate(CYPS):
                best = -1e9
                for cid in cids:
                    try:
                        o3 = rdMolAlign.GetCrippenO3A(m, refs[c], par, ref_par[c],
                                                      prbCid=int(cid), refCid=-1)
                        o3.Align()
                        best = max(best, o3.Score())
                    except Exception:
                        continue
                if best > -1e8:
                    S[i, e] = best
            if (i + 1) % 500 == 0:
                print(f"    {tag}: {i+1}/{len(smiles)}  ({time.time()-t0:.0f} с)", flush=True)
        out[tag] = S
        print(f"  {tag}: {len(smiles)} молекул, не получилось {bad}, "
              f"NaN ячеек {int(np.isnan(S).sum())}  ({time.time()-t0:.0f} с)", flush=True)

    out["names"] = np.array([f"o3a_{c}" for c in CYPS], dtype=object)
    np.savez_compressed(a.out, **out)
    S = out["train"]
    print(f"\nсохранено: {a.out}")
    print(f"\n{'фермент':8s} {'медиана':>9s} {'sd':>8s} {'корр. с числом атомов':>22s}")
    n_heavy = np.array([Chem.MolFromSmiles(s).GetNumHeavyAtoms() if Chem.MolFromSmiles(s)
                        else np.nan for s in rows.SMILES], float)
    ok = np.isfinite(S).all(1) & np.isfinite(n_heavy)
    for e, c in enumerate(CYPS):
        r = np.corrcoef(S[ok, e], n_heavy[ok])[0, 1]
        print(f"{c:8s} {np.nanmedian(S[:, e]):9.2f} {np.nanstd(S[:, e]):8.2f} {r:22.3f}")
    C = S[ok] - S[ok].mean(1, keepdims=True)
    print(f"\nпосле центрирования (контраст) корреляция с числом атомов:")
    for e, c in enumerate(CYPS):
        print(f"  {c}: {np.corrcoef(C[:, e], n_heavy[ok])[0, 1]:+.3f}")
    print("""
Последние две таблицы --- предварительная проверка замысла. Сырые оценки обязаны сильно
коррелировать с размером молекулы: это и есть общая часть, которую блок лиганда уже несёт.
Если после центрирования корреляция падает почти до нуля, контраст действительно вычел
размер, и остаток --- про полость. Если не падает, признак переодетый объём, и абляцию можно
не запускать.""")


if __name__ == "__main__":
    main()
