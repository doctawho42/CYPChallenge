"""Rule-based pKa of the most basic centre, with a Taft-style inductive correction.

Not a research-grade pKa predictor: the point is a *chemically honest* number for the
one thing CYP2D6 cares about - is there a nitrogen that carries a positive charge at
pH 7.4, and how strongly. Validated against literature pKa for known drugs below.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')

# (name, SMARTS matching the basic atom as the FIRST mapped atom, intrinsic pKa of the
#  unsubstituted parent).  Order matters: first match wins for a given atom.
CLASSES = [
    ("guanidine",   "[NX3;!$(N[CX3]=[OX1])][CX3](=[NX2;!$(N[CX3]=[OX1])])[NX3]", 13.0),
    ("amidine",     "[NX2;!$(N[CX3]=[OX1])]=[CX3][NX3]",                          11.5),
    ("amine_2ary",  "[NX3;H1;!$(N[a]);!$(N[CX3]=[OX1,SX1,NX2]);!$(N[SX4]);!$(N[PX4]);!$(N#*);!$(N=*);!$([N+])]([CX4])[CX4]", 10.9),
    ("amine_1ary",  "[NX3;H2;!$(N[a]);!$(N[CX3]=[OX1,SX1,NX2]);!$(N[SX4]);!$(N[PX4])][CX4]", 10.6),
    ("amine_3ary",  "[NX3;H0;!$(N[a]);!$(N[CX3]=[OX1,SX1,NX2]);!$(N[SX4]);!$(N[PX4]);!$(N#*);!$(N=*);!$([N+])]([CX4])([CX4])[CX4]", 9.8),
    ("imidazole",   "[nX2;H0;$(n1cncc1),$(n1cnc[nH]1),$(n1cncn1)]",                7.0),
    ("pyridine",    "[nX2;H0]",                                                    5.2),
    ("aniline",     "[NX3;!$(N[CX3]=[OX1]);!$(N[SX4]);!$(N#*);!$(N=*)][a]",        4.6),
]
CLASSES = [(n, Chem.MolFromSmarts(s), p) for n, s, p in CLASSES]

# Inductive electron withdrawal, damped by topological distance from the basic nitrogen.
_W = {2: 1.55, 3: 0.62, 4: 0.20}

def _parent_ring_atoms(mol, n_idx, cls):
    """Atoms already accounted for by the class's intrinsic pKa - must not be counted twice.

    For pyridine/imidazole the intrinsic value is that of the parent heteroarene, so the
    whole aromatic ring system containing N is already in the number. For an aniline the
    parent is aniline itself, so the aryl ring the N hangs off is already in the number.
    """
    if cls not in ("pyridine", "imidazole", "aniline"): return set()
    out = set()
    ri = mol.GetRingInfo()
    if cls in ("pyridine", "imidazole"):
        for ring in ri.AtomRings():
            if n_idx in ring: out.update(ring)
    else:
        for nb in mol.GetAtomWithIdx(n_idx).GetNeighbors():
            if nb.GetIsAromatic():
                for ring in ri.AtomRings():
                    if nb.GetIdx() in ring: out.update(ring)
    return out

def _ewg_penalty(mol, n_idx, dm, skip=frozenset()):
    pen = 0.0
    for a in mol.GetAtoms():
        i = a.GetIdx()
        if i == n_idx or i in skip: continue
        d = int(dm[n_idx][i])
        if d not in _W: continue
        sym, arom = a.GetSymbol(), a.GetIsAromatic()
        s = 0.0
        if sym in ("F",):                      s = 1.0
        elif sym in ("Cl", "Br", "I"):         s = 0.75
        elif sym == "O":                       s = 1.15
        elif sym == "N" and not arom:          s = 0.55
        elif sym == "N" and arom:              s = 0.5
        elif sym == "S":                       s = 0.5
        elif sym == "C" and arom:              s = 0.45
        elif sym == "C" and any(b.GetBondTypeAsDouble() >= 2 for b in a.GetBonds()): s = 0.35
        pen += _W[d] * s
    return pen

def most_basic_pka(mol):
    """(pKa of the most basic centre, its atom index, class name). pKa = -99 if none."""
    dm = Chem.GetDistanceMatrix(mol)
    best = (-99.0, -1, "none")
    seen = set()
    for name, patt, base in CLASSES:
        for match in mol.GetSubstructMatches(patt):
            i = match[0]
            if i in seen: continue
            seen.add(i)
            pk = base - _ewg_penalty(mol, i, dm, _parent_ring_atoms(mol, i, name))
            if pk > best[0]: best = (pk, i, name)
    return best

def frac_protonated(pka, ph=7.4):
    if pka < -50: return 0.0
    return 1.0 / (1.0 + 10 ** (ph - pka))

if __name__ == "__main__":
    # literature pKa of the most basic centre
    REF = [
        ("этиламин",        "CCN",                                                        10.7),
        ("бензиламин",      "NCc1ccccc1",                                                  9.3),
        ("анилин",          "Nc1ccccc1",                                                   4.6),
        ("пиридин",         "c1ccncc1",                                                    5.2),
        ("имидазол",        "c1c[nH]cn1",                                                  7.0),
        ("дифенгидрамин",   "CN(C)CCOC(c1ccccc1)c1ccccc1",                                 9.0),
        ("пропранолол",     "CC(C)NCC(O)COc1cccc2ccccc12",                                 9.5),
        ("амитриптилин",    "CN(C)CCC=C1c2ccccc2CCc2ccccc21",                              9.4),
        ("флуоксетин",      "CNCCC(Oc1ccc(C(F)(F)F)cc1)c1ccccc1",                         10.1),
        ("верапамил",       "COc1ccc(CCN(C)CCCC(C#N)(C(C)C)c2ccc(OC)c(OC)c2)cc1OC",        8.9),
        ("морфолин",        "C1COCCN1",                                                    8.4),
        ("кофеин",          "Cn1c(=O)c2c(ncn2C)n(C)c1=O",                                  0.6),
        ("парацетамол",     "CC(=O)Nc1ccc(O)cc1",                                         -1.0),
        ("хинидин",         "C=CC1CN2CCC1CC2C(O)c1ccnc2ccc(OC)cc12",                       8.7),
        ("пиперидин",       "C1CCNCC1",                                                   11.1),
        ("N-метилпиперазин","CN1CCNCC1",                                                   9.1),
        ("хинуклидин",      "C1CN2CCC1CC2",                                               11.0),
        ("лидокаин",        "CCN(CC)CC(=O)Nc1c(C)cccc1C",                                  7.9),
        ("никотин",         "CN1CCCC1c1cccnc1",                                            8.0),
        ("прокаин",         "CCN(CC)CCOC(=O)c1ccc(N)cc1",                                  9.0),
    ]
    print(f"{'соединение':<16}{'класс':<12}{'оценка':>8}{'лит.':>8}{'ошибка':>8}")
    errs = []
    for name, smi, lit in REF:
        m = Chem.MolFromSmiles(smi)
        pk, idx, cls = most_basic_pka(m)
        pk_show = pk if pk > -50 else -1.0
        e = pk_show - lit; errs.append(abs(e))
        print(f"{name:<16}{cls:<12}{pk_show:>8.1f}{lit:>8.1f}{e:>+8.1f}")
    print(f"\nсредняя абсолютная ошибка: {np.mean(errs):.2f} ед. pKa (n={len(errs)})")
    print("доля правильно классифицированных как протонированные при 7.4:",
          f"{np.mean([(most_basic_pka(Chem.MolFromSmiles(s))[0] > 7.4) == (l > 7.4) for _, s, l in REF]):.2f}")
