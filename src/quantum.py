"""A semi-empirical block aimed at type II heme coordination, not at more columns.

The mechanism this is for. CYP inhibition has more than one mode, and the one our matrix has no
handle on is coordination of the heme iron by a lone pair on an sp2 nitrogen. It is why the
azoles -- ketoconazole, fluconazole -- are the reference inhibitors of CYP3A4 and CYP1A2. Our
mechanistic block carries pKa as a proxy for basicity, which is about *protonation*, and
protonation is the opposite of what coordination needs: a protonated nitrogen has no lone pair to
donate. Nothing in the matrix says whether the lone pair is available, how easily it ionises, or
whether anything can reach it.

What is computed, and why each one rather than a descriptor dump:

  HOMO, LUMO, gap        the frontier orbitals. Donation to an empty iron d orbital is a
                         HOMO-controlled interaction, so its energy is the first-order term.
  заряд на самом основном азоте
                         local electron density at the atom that would coordinate.
  минимальный заряд по ароматическим N
                         the most electron-rich aromatic nitrogen in the molecule, which is the
                         candidate donor whether or not it is the most basic by pKa.
  f- по Фукуи на нём     the local softness: how much density that atom actually loses on
                         ionisation, computed as q_cation - q_neutral. This is the part pKa
                         cannot proxy, and it costs a second single-point.
  диполь                 orientation term, cheap, and it comes out of the same call.
  свободный конус        geometric: the fraction of directions around the candidate nitrogen not
                         blocked by the molecule's own atoms within 4 A. A perfectly available
                         lone pair on a buried nitrogen coordinates nothing, and this is the
                         steric half of the mechanism that no electronic descriptor carries.

Pre-registered reading, and it separates mechanism from padding. The gain must appear on CYP3A4
and CYP1A2, where the azole story is the reference chemistry, and must be **absent on CYP2D6**,
whose pharmacophore is a protonated basic amine forming a salt bridge -- the opposite
requirement. A uniform gain across four enzymes means forty extra columns helped a boosting, not
that coordination chemistry was the missing axis, and should be read as a negative.

The prior is honestly low. Items 96, 117 and 119 each added a feature family and each failed, and
this project's record on "more features" is poor. What distinguishes this one is that a mechanism
predicts *where* it should work and where it should not, so it can fail informatively.

Runs in a separate environment, like src/embed.py and src/embed2.py, and for the same reason:

    uv pip install --python /tmp/enc2/bin/python tblite
    /tmp/enc2/bin/python src/quantum.py

Writes data/quantum.npz with `train` (n_train, d) in rows.csv order, `test` (750, d) in the
blinded file's order, and `names`.
"""
import argparse
import pathlib
import time

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

RDLogger.DisableLog("rdApp.*")
ROOT = pathlib.Path(__file__).resolve().parents[1]
D = str(ROOT / "data") + "/"

NAMES = ["homo", "lumo", "gap", "dipole", "q_basicN", "q_aromN_min",
         "fukui_minus", "cone_free", "n_arom_N", "has_donor"]
BASIC_N = Chem.MolFromSmarts("[nX2,NX2;!$(N=O)]")     # sp2 азот с неподелённой парой
HARTREE = 27.211386


def geometry(smi):
    """Один низкоэнергетический конформер: тот же рецепт, что в src/shape3d.py.

    Настройки те же (ETKDGv3, сид 0xC0FFEE, MMFF 400 итераций), но rdkit здесь другой ---
    отдельное окружение поставило свою версию. Значит конформеры этого блока НЕ обязаны
    совпадать с конформерами shape3d покоординатно, и признак свободного конуса нельзя
    сравнивать с геометрическими признаками того блока как посчитанные на одной геометрии.
    Внутри блока всё согласовано, между блоками --- нет, и это надо помнить, если они
    когда-нибудь окажутся в одной таблице.
    """
    m = Chem.MolFromSmiles(smi)
    if m is None:
        return None, None
    mh = Chem.AddHs(m)
    ps = AllChem.ETKDGv3()
    ps.randomSeed = 0xC0FFEE
    if AllChem.EmbedMolecule(mh, ps) != 0:
        return None, None
    try:
        AllChem.MMFFOptimizeMolecule(mh, maxIters=400)
    except Exception:
        pass
    return mh, mh.GetConformer()


def cone_free(pos, idx, cutoff=4.0, ndir=200):
    """Доля направлений вокруг атома, не перекрытых собственными атомами молекулы.

    Грубая, но нужная величина: электронный дескриптор ничего не говорит о том, дотянется ли
    железо до пары. Направления берутся с равномерной сферы, направление считается занятым,
    если в конусе 30 градусов вокруг него есть атом ближе cutoff.
    """
    c = pos[idx]
    d = pos - c
    r = np.linalg.norm(d, axis=1)
    keep = (r > 0.1) & (r < cutoff)
    if not keep.any():
        return 1.0
    u = d[keep] / r[keep, None]
    g = np.random.default_rng(0).standard_normal((ndir, 3))
    g /= np.linalg.norm(g, axis=1, keepdims=True)
    blocked = (g @ u.T > np.cos(np.deg2rad(30))).any(axis=1)
    return float(1.0 - blocked.mean())


def one(smi, calc_cls):
    mh, conf = geometry(smi)
    if mh is None:
        return None
    pos_ang = np.array([list(conf.GetAtomPosition(i)) for i in range(mh.GetNumAtoms())])
    nums = np.array([a.GetAtomicNum() for a in mh.GetAtoms()])
    pos_bohr = pos_ang / 0.52917721092

    def run(charge):
        c = calc_cls("GFN2-xTB", nums, pos_bohr, charge=float(charge),
                     uhf=int(abs(charge) % 2))
        c.set("verbosity", 0)
        return c.singlepoint()

    try:
        res = run(0)
    except Exception:
        return None
    eps = np.asarray(res.get("orbital-energies"), float)
    occ = np.asarray(res.get("orbital-occupations"), float)
    q0 = np.asarray(res.get("charges"), float)
    dip = float(np.linalg.norm(np.asarray(res.get("dipole"), float)))
    filled = occ > 0.5
    homo = float(eps[filled].max() * HARTREE) if filled.any() else 0.0
    lumo = float(eps[~filled].min() * HARTREE) if (~filled).any() else 0.0

    # Кандидат-донор: самый электронно-богатый sp2 азот.
    cand = [a[0] for a in mh.GetSubstructMatches(BASIC_N)]
    arom_n = [i for i in range(mh.GetNumAtoms())
              if nums[i] == 7 and mh.GetAtomWithIdx(i).GetIsAromatic()]
    if cand:
        j = int(min(cand, key=lambda i: q0[i]))
    elif arom_n:
        j = int(min(arom_n, key=lambda i: q0[i]))
    else:
        j = -1

    fukui = 0.0
    if j >= 0:
        try:
            qc = np.asarray(run(1).get("charges"), float)
            fukui = float(qc[j] - q0[j])
        except Exception:
            fukui = 0.0

    return [homo, lumo, lumo - homo, dip,
            float(q0[j]) if j >= 0 else 0.0,
            float(min((q0[i] for i in arom_n), default=0.0)),
            fukui,
            cone_free(pos_ang, j) if j >= 0 else 0.0,
            float(len(arom_n)), float(j >= 0)]


def build(smiles, tag, calc_cls):
    out, bad, t0 = [], 0, time.time()
    for i, s in enumerate(smiles):
        v = one(s, calc_cls)
        if v is None:
            v = [0.0] * len(NAMES)
            bad += 1
        out.append(v)
        if (i + 1) % 250 == 0:
            print(f"  {tag} {i+1}/{len(smiles)}, не сошлось {bad}, "
                  f"{time.time()-t0:.0f} с", flush=True)
    return np.array(out, np.float32), bad


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=D + "quantum.npz")
    ap.add_argument("--limit", type=int, default=0)
    a = ap.parse_args()
    from tblite.interface import Calculator

    rows = pd.read_csv(D + "rows.csv")
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    tr_smi, te_smi = list(rows.SMILES), list(te.SMILES)
    if a.limit:
        tr_smi, te_smi = tr_smi[:a.limit], te_smi[:a.limit]

    A, bad_a = build(tr_smi, "обучение", Calculator)
    B, bad_b = build(te_smi, "тест", Calculator)
    print(f"обучение {A.shape}, не сошлось {bad_a}; тест {B.shape}, не сошлось {bad_b}")

    sd = A.std(0)
    dead = int((sd < 1e-9).sum())
    print("колонки: " + " ".join(f"{n}({s:.3f})" for n, s in zip(NAMES, sd)))
    if dead > len(NAMES) // 2:
        raise SystemExit("больше половины колонок стоит --- блок вырожден")
    np.savez_compressed(a.out, train=A, test=B, names=np.array(NAMES))
    print(f"сохранено: {a.out}")


if __name__ == "__main__":
    main()
