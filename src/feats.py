"""Mechanistic feature block for CYP2D6 + generic feature builders. Caches to npz."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import re, os, time, numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import rdFingerprintGenerator, Descriptors, Crippen, rdMolDescriptors
RDLogger.DisableLog('rdApp.*')

from pka import most_basic_pka, frac_protonated
# Ship-along SMARTS table of dimorphite's protonation sites. Taken from the installed
# package rather than a hard-coded site-packages path, which only existed in the
# original Linux image.
import dimorphite_dl
SM = str(_pl.Path(dimorphite_dl.__file__).parent / "smarts" / "site_substructures.smarts")

def load_pka_table():
    """name -> (mol, [pka_mean,...]) for every protonation site dimorphite knows."""
    out = []
    for line in open(SM):
        line = line.strip()
        if not line or line.startswith("#"): continue
        f = line.split("\t")
        name, smarts = f[0], f[1]
        pkas = [float(f[i+1]) for i in range(2, len(f), 3) if abs(float(f[i+1])) < 900]
        m = Chem.MolFromSmarts(smarts)
        if m is not None and pkas: out.append((name, m, pkas))
    return out
PKA = load_pka_table()
# sites whose protonated form is a CATION (bases). Everything else is an acid.
# Exactly the sites whose protonated form carries the positive charge (Brønsted bases).
# Everything else in dimorphite's table is an acid (loses H+ to give an anion), and its
# pKa must NOT enter a "most basic centre" statistic.
BASIC = {"AmidineGuanidine1", "AmidineGuanidine2", "Anilines_primary", "Anilines_secondary",
         "Anilines_tertiary", "Aromatic_nitrogen_unprotonated", "*Aromatic_nitrogen_protonated",
         "Amines_primary_secondary_tertiary", "Primary_hydroxyl_amine"}

SMARTS = {
 "basicN_ali":  "[NX3;H2,H1,H0;!$(N[#6]=[O,N,S]);!$(N[a]);!$(N[SX4](=O)=O);!$(N#*);!$([N+]);!$(N=*)]",
 "amine_prim":  "[NX3;H2;!$(N[#6]=[O,N,S]);!$(N[a]);!$(N[SX4](=O)=O)]",
 "amine_sec":   "[NX3;H1;!$(N[#6]=[O,N,S]);!$(N[a]);!$(N[SX4](=O)=O);!$(N=*)]",
 "amine_tert":  "[NX3;H0;!$(N[#6]=[O,N,S]);!$(N[a]);!$(N[SX4](=O)=O);!$(N=*);!$([N+])]",
 "amidine":     "[NX3][CX3]=[NX2]",
 "guanidine":   "[NX3][CX3](=[NX2])[NX3]",
 "pyridineN":   "[nX2;$(n1ccccc1),$(n1ccncc1)]",
 "imidazoleN":  "[nX2;$(n1ccnc1)]",
 "anilineN":    "[NX3;!$(N[#6]=[O,N,S])][a]",
 "amideN":      "[NX3][CX3]=[OX1]",
 "morpholine":  "C1COCCN1",
 "piperazine":  "C1CNCCN1",
 "piperidine":  "[NX3;R]1[CX4][CX4][CX4][CX4][CX4]1",
 "acid":        "[CX3](=O)[OX2H1]",
 "tetrazole":   "c1nnn[nH]1",
 "sulfonamide": "[SX4](=O)(=O)[NX3]",
}
SMARTS = {k: Chem.MolFromSmarts(v) for k, v in SMARTS.items()}

def basic_pka_profile(mol):
    """Highest pKa among basic (cation-forming) sites, and how many such sites."""
    best, n = -99.0, 0
    for name, patt, pkas in PKA:
        if name not in BASIC: continue
        hits = mol.GetSubstructMatches(patt)
        if not hits: continue
        n += len(hits)
        best = max(best, max(pkas))
    return best, n

def n_aromatic_centroid_topo(mol):
    """Topological bond distance from every candidate basic N to the nearest aromatic ring atom."""
    ri = mol.GetRingInfo()
    arom_atoms = [a.GetIdx() for a in mol.GetAtoms() if a.GetIsAromatic()]
    bN = [i for (i,) in mol.GetSubstructMatches(SMARTS["basicN_ali"])] if SMARTS["basicN_ali"] else []
    bN += [m[0] for m in mol.GetSubstructMatches(SMARTS["amidine"])]
    bN = sorted(set(bN))
    if not bN or not arom_atoms: return (-1.0, -1.0, 0)
    dm = Chem.GetDistanceMatrix(mol)
    ds = [min(dm[i][j] for j in arom_atoms) for i in bN]
    return (float(min(ds)), float(np.mean(ds)), len(bN))

def mech_block(mols):
    rows = []
    for m in mols:
        d = {}
        for k, p in SMARTS.items():
            d["n_" + k] = len(m.GetSubstructMatches(p))
        pka, idx, cls = most_basic_pka(m)
        d["pka_max_basic"] = pka if pka > -50 else 0.0
        d["frac_prot_74"] = frac_protonated(pka)
        d["is_base_74"] = 1.0 if pka > 7.4 else 0.0
        _, nsite = basic_pka_profile(m)
        d["n_basic_sites"] = nsite
        dmin, dmean, nb = n_aromatic_centroid_topo(m)
        d["topo_bN_to_arom_min"] = dmin
        d["topo_bN_to_arom_mean"] = dmean
        d["n_bN_geom"] = nb
        # distance from the MOST BASIC centre to the nearest aromatic ring atom
        if idx >= 0:
            ar = [a.GetIdx() for a in m.GetAtoms() if a.GetIsAromatic()]
            dmx = Chem.GetDistanceMatrix(m)
            d["topo_mbc_to_arom"] = float(min((dmx[idx][j] for j in ar), default=-1.0))
        else:
            d["topo_mbc_to_arom"] = -1.0
        # CYP2D6 pharmacophore: cation at pH 7.4 sitting 2-5 bonds (~5-7 A) from an aryl
        d["pharm_2d6"] = 1.0 if (2 <= d["topo_mbc_to_arom"] <= 5) else 0.0
        d["pharm_2d6_x_prot"] = d["pharm_2d6"] * d["frac_prot_74"]
        d["formal_charge"] = float(Chem.GetFormalCharge(m))
        rows.append(d)
    return pd.DataFrame(rows)

def dimorphite_charge(smiles_list, ph=7.4):
    from dimorphite_dl import protonate_smiles
    net, npos, nneg = [], [], []
    for s in smiles_list:
        try:
            out = protonate_smiles(s, ph_min=ph, ph_max=ph, precision=0.0)
            m = Chem.MolFromSmiles(out[0]) if out else None
        except Exception:
            m = None
        if m is None:
            net.append(np.nan); npos.append(np.nan); nneg.append(np.nan); continue
        chs = [a.GetFormalCharge() for a in m.GetAtoms()]
        net.append(sum(chs)); npos.append(sum(1 for c in chs if c > 0)); nneg.append(sum(1 for c in chs if c < 0))
    return pd.DataFrame({"ph74_net_charge": net, "ph74_n_cation": npos, "ph74_n_anion": nneg})

def build(smiles, desc_names=None, mech_names=None):
    """FP, DESC and MECH for an arbitrary SMILES list, plus the kept molecule indices.

    Pass desc_names / mech_names to reindex the descriptor and mechanistic blocks onto a
    fixed column set. That is what the test set needs: the training build selects
    descriptor columns by `isna().mean() < 0.05`, a filter fitted on the training
    molecules, and recomputing it on 750 test molecules would silently produce a
    different column set. Checking the column *count* would not catch a permutation
    either, so the selection has to be by name.
    """
    mols = [Chem.MolFromSmiles(s) for s in smiles]
    ok = [i for i, m in enumerate(mols) if m is not None]
    mols = [mols[i] for i in ok]
    kept = [smiles[i] for i in ok]

    M = pd.concat([mech_block(mols), dimorphite_charge(kept)], axis=1).fillna(0.0)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    FP = np.array([gen.GetCountFingerprintAsNumPy(m) for m in mols], dtype=np.float32)
    dsc = pd.DataFrame([Descriptors.CalcMolDescriptors(m) for m in mols]).replace([np.inf, -np.inf], np.nan)

    if desc_names is None:
        dsc = dsc.loc[:, dsc.isna().mean() < 0.05]
    else:
        dsc = dsc.reindex(columns=list(desc_names))
    if mech_names is not None:
        M = M.reindex(columns=list(mech_names))
    dsc = _guard(dsc)
    return FP, dsc.fillna(0.0), M.fillna(0.0), ok


# float32 хранит до 3.4e38, а стандартизация возводит в квадрат, поэтому опасный порог
# для дисперсии --- корень из него. 1e12 --- порог другого рода: выше него колонка для
# дерева уже вырождена по рангу, и это стоит знать, но ронять сборку из-за этого нельзя.
VAR_OVERFLOW = float(np.sqrt(np.finfo(np.float32).max))
RANK_DEGENERATE = 1e12


def _guard(dsc):
    """Catch descriptor columns that overflow float32, and say which.

    The existing `replace([inf, -inf], nan)` runs on float64 values and therefore misses the
    case that actually occurs: a value finite in float64 that becomes infinite when cast to
    float32 downstream. Ipc does this - it grows exponentially with molecule size, reaches
    1.5e36 on public ChEMBL compounds and overflows outright on one of them, while our own
    training set tops out at 5.1e14 and never triggers it. The cleaning has to happen after
    the cast that creates the problem, not before it.

    Columns are cleared to NaN rather than clipped, because a value this large carries no
    usable information either way and clipping would invent a number. RDKit ships `AvgIpc`
    for exactly this reason; swapping it in would be the real fix and would move every
    published number, so it is a deliberate change and not this function's business.
    """
    a = dsc.to_numpy(np.float64)
    bad = np.isfinite(a) & ~np.isfinite(a.astype(np.float32))
    if bad.any():
        for j in np.where(bad.any(0))[0]:
            print(f"  ВНИМАНИЕ: {dsc.columns[j]} переполняет float32 в "
                  f"{int(bad[:, j].sum())} строках, эти ячейки обнулены")
        dsc = dsc.mask(pd.DataFrame(bad, index=dsc.index, columns=dsc.columns))
    mx = np.nanmax(np.abs(np.where(np.isfinite(a), a, np.nan)), axis=0)
    for j in np.where(mx > RANK_DEGENERATE)[0]:
        lvl = "переполнит дисперсию" if mx[j] > VAR_OVERFLOW else "вырождена по рангу"
        print(f"  ВНИМАНИЕ: {dsc.columns[j]} доходит до {mx[j]:.2e}, {lvl}")
    return dsc


if __name__ == "__main__":

    tr = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
    t = time.time()
    FP, dsc, M, ok = build(list(tr.SMILES))
    tr = tr.loc[ok].reset_index(drop=True)
    print("features", round(time.time() - t, 1), "s")
    np.savez_compressed(D + "feats.npz", FP=FP, DESC=dsc.to_numpy(np.float32), MECH=M.to_numpy(np.float32))
    pd.Series(list(dsc.columns)).to_csv(D + "desc_names.csv", index=False, header=False)
    pd.Series(list(M.columns)).to_csv(D + "mech_names.csv", index=False, header=False)
    tr[["Molecule_Name", "SMILES"]].to_csv(D + "rows.csv", index=False)
    print("FP", FP.shape, "DESC", dsc.shape, "MECH", M.shape)
    print(M.describe().T[["mean", "std", "min", "max"]].round(2).to_string())
