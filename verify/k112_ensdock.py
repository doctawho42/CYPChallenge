"""Ensemble docking into experimental CYP3A4 receptors, and whether contact consensus beats the
docking score at choosing a pose. The decisive first measurement, not the campaign.

THE IDEA, and why it is not a foundation model. The structure track asks for twenty CYP3A4-ligand
complexes from ligand SMILES. A co-folding model re-predicts the entire 503-residue fold for every
ligand; but the protein is THE SAME PROTEIN in all twenty tasks, and the PDB already holds 129
experimental answers to that half (item 333). So the protein is taken, not predicted, and what
remains is placing a ligand in a fixed pocket -- which smina does in about 32 seconds per run
(item 306's campaign: 22608 runs, 204 CPU-hours).

WHAT REPLACES THE MISSING HALF. One rigid receptor biases every pose toward the shape its own
ligand carved, which is the likely reason a serious rigid-docking entry (ADMETForge, Matcha) scores
BELOW the stock Boltz-2 baseline. The fix needs no protein prediction: dock into MANY experimental
conformations and choose among the results. The choosing is the part with a claim to being ours.

WHY NOT THE DOCKING SCORE. Scoring functions rank poses badly -- that is textbook -- and the metric
here is LDDT-PLI, which counts PRESERVED CONTACTS. So the criterion should be contact agreement, not
predicted affinity: among all (receptor, pose) candidates, prefer the one whose protein-ligand
contact set is most reproduced by the others. The ensemble serves as its own uncertainty estimate,
and item 333 measured why that matters -- for a third of the same-ligand pairs in the PDB the SAME
ligand adopts a drastically different pose, so the modal pose is a better bet than the best-scoring
one.

WHAT THIS FILE MEASURES, and it is deliberately small. Two target ligands chosen to straddle item
333's bimodality -- MYT, which reproduces at 0.443 A between its own two crystals, and PK9, which
differs by 6.247 A between its own two -- docked LEAVE-ONE-OUT into receptors that are not their own,
then every pose scored against the held-out crystal. Three selection rules are compared on the same
pose set:

    по скору smina     the baseline anyone would use
    по согласию контактов   ours
    оракул             the best pose present, which bounds both

If consensus does not beat the score here, the idea is wrong and the campaign is not run.

HONEST LIMITS, stated before the numbers. (1) Two ligands is not a sample; this is a go/no-go, and
the campaign that follows is what would carry a claim. (2) The metrics are k111's local
implementations, whose absolute agreement with OpenStructure is unverified -- they rank, they do not
predict a board value. (3) Iridium complexes in the PDB's CYP3A4 set (X2Q, WZN) are excluded: they
are photo-probes, not drug-like, and the twenty targets are not like them.

    uv run python verify/k112_ensdock.py --targets MYT,PK9

Writes results/logs/k112_ensdock.json. Exit: 0 ran, 2 inputs missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES

import argparse
import json
import os
import subprocess
import tempfile
import urllib.parse
import urllib.request

import numpy as np
from rdkit import Chem, RDLogger
from rdkit.Chem import AllChem

from k111_structscore import (INCLUSION, SITE_RADIUS, THRESH, fetch, kabsch, parse,
                              split_complex)

RDLogger.DisableLog("rdApp.*")
SMINA = os.path.expanduser("~/smina")

# catalytic-site holo entries from item 333's listing, metal complexes excluded
ENSEMBLE = {"RIT": ["3NXU", "5VC0", "5VCE", "38FV"], "MYT": ["1W0G", "6MA6"],
            "PK9": ["4D6Z", "4D75"], "08Y": ["3UA1", "5VCG"], "FZV": ["6DAC", "6DAJ"],
            "G0D": ["6DA2", "6DA3"], "CFF": ["8SO1", "8SO2"], "MWY": ["38FX", "6OOA"],
            "KLN": ["2V0M", "38FW"], "TPF": ["6MA7"]}
EXCLUDE_LIG = {"X2Q", "WZN"}          # iridium photo-probes, not drug-like


def comp_smiles(ids):
    q = ('{chem_comps(comp_ids:[%s]){chem_comp{id} '
         'rcsb_chem_comp_descriptor{SMILES_stereo}}}' % ",".join('"%s"' % i for i in ids))
    d = json.load(urllib.request.urlopen(
        "https://data.rcsb.org/graphql?query=" + urllib.parse.quote(q)))["data"]["chem_comps"]
    return {c["chem_comp"]["id"]: (c.get("rcsb_chem_comp_descriptor") or {}).get("SMILES_stereo")
            for c in d}


def pdb_line(i, a, hetero):
    rec = "HETATM" if hetero else "ATOM  "
    nm = a["name"] if len(a["name"]) >= 4 else " " + a["name"].ljust(3)
    x, y, z = a["xyz"]
    return (f"{rec}{i:5d} {nm[:4]} {a['res']:>3} A{int(a['seq']) % 10000:4d}    "
            f"{x:8.3f}{y:8.3f}{z:8.3f}  1.00  0.00          {a['el']:>2}")


def heme_box(prot, size=18.0, lift=5.0):
    """A search box anchored on the HEME, identical for every receptor.

    The first version used `--autobox_ligand` on each receptor's own co-crystal ligand, which makes
    the search volume depend on which drug happened to be deposited: ritonavir's box is
    11 x 13 x 11 A and caffeine's 1.3 x 5.7 x 5.8, so metyrapone was searched in a volume cut for a
    721-dalton ligand and the pose set came out mediocre everywhere (oracle 3.04 A against a 1.45 A
    ceiling). The catalytic site is defined by the heme, not by a guest: centre is the iron lifted
    along the porphyrin normal AWAY from the proximal cysteine thiolate, which is the distal face
    where catalysis happens.
    """
    fe = [a["xyz"] for a in prot if a["res"] == "HEM" and a["el"] == "FE"]
    ns = [a["xyz"] for a in prot if a["res"] == "HEM" and a["el"] == "N"]
    if not fe or len(ns) < 3:
        return None
    fe = fe[0]
    P = np.array(ns) - fe
    _, _, Vt = np.linalg.svd(P - P.mean(0))
    nrm = Vt[2] / np.linalg.norm(Vt[2])
    sg = [a["xyz"] for a in prot if a["name"] == "SG"]
    if sg:
        near = min(sg, key=lambda x: np.linalg.norm(x - fe))
        if np.dot(near - fe, nrm) > 0:      # normal points at the thiolate: flip to the distal side
            nrm = -nrm
    return fe + lift * nrm, size


def receptor_files(pdb_id, lig_id, work):
    """Protein chain + its heme as PDB (ligand stripped), and the heme-anchored box."""
    prot, lig, d_fe = split_complex(parse(fetch(pdb_id)), lig_id)
    if not lig:
        return None
    rec = work / f"{pdb_id}_rec.pdb"
    rec.write_text("\n".join(pdb_line(i, a, a["res"] == "HEM")
                             for i, a in enumerate(prot, 1)) + "\nEND\n")
    hb = heme_box(prot)
    if hb is None:
        return None
    return rec, hb, prot, lig, d_fe


def lig_mol(lig_atoms, smiles):
    """RDKit mol of a crystal ligand, bond orders from its SMILES template."""
    blk = "\n".join(pdb_line(i, a, True) for i, a in enumerate(lig_atoms, 1)) + "\nEND\n"
    m = Chem.MolFromPDBBlock(blk, removeHs=False, sanitize=False)
    if m is None:
        return None
    ref = Chem.MolFromSmiles(smiles)
    if ref is None:
        return None
    try:
        return AllChem.AssignBondOrdersFromTemplate(ref, m)
    except Exception:
        return m


def embed3d(smiles, path):
    m = Chem.AddHs(Chem.MolFromSmiles(smiles))
    ps = AllChem.ETKDGv3(); ps.randomSeed = 42
    if AllChem.EmbedMolecule(m, ps) != 0:
        return False
    AllChem.MMFFOptimizeMolecule(m, maxIters=500)
    w = Chem.SDWriter(str(path)); w.write(m); w.close()
    return True


def dock(rec, box, lig_sdf, out, modes=9, exh=8):
    (cx, cy, cz), sz = box
    cmd = [SMINA, "-r", str(rec), "-l", str(lig_sdf),
           "--center_x", f"{cx:.3f}", "--center_y", f"{cy:.3f}", "--center_z", f"{cz:.3f}",
           "--size_x", f"{sz:.1f}", "--size_y", f"{sz:.1f}", "--size_z", f"{sz:.1f}",
           "--exhaustiveness", str(exh), "--num_modes", str(modes),
           "--cpu", "1", "--seed", "42", "-o", str(out)]
    r = subprocess.run(cmd, capture_output=True, text=True, timeout=900)
    return r.returncode == 0 and _pl.Path(out).exists()


def poses_from(sdf):
    """Every pose with its smina affinity. sanitize=False is load-bearing (item 304): smina writes
    poses back without charges and RDKit returns None on a strict read, which silently drops them."""
    out = []
    for m in Chem.SDMolSupplier(str(sdf), sanitize=False, removeHs=False):
        if m is None:
            continue
        try:
            aff = float(m.GetProp("minimizedAffinity"))
        except Exception:
            aff = float("nan")
        out.append((m, aff))
    return out


def correspondence(tmpl, ref_mol, pose_mol):
    """Pose-index -> reference-index maps, with symmetry, independent of atom order.

    The first version assumed smina preserves the input SDF's heavy-atom order. It does for
    metyrapone and NOT for PK9 or ritonavir: smina keeps polar hydrogens (25 atoms out of a
    23-heavy input) and reorders. The assertion caught it -- 171 of 171 poses refused rather than
    mismatched -- but the assumption was wrong, so the order is no longer used. Instead the pose
    gets its bond orders from the SMILES template, exactly as the organisers' own
    `posebusters_utils.assign_bond_orders_from_template` does and for the same reason (RDKit's PDB
    and smina's SDF both leave bonds unresolved), and every template embedding is then one
    candidate correspondence. Measured match counts: MYT 2, PK9 6, RIT 16 -- that IS the symmetry.
    """
    try:
        noh = Chem.RemoveHs(pose_mol, sanitize=False)
    except Exception:
        noh = pose_mol
    try:
        fixed = AllChem.AssignBondOrdersFromTemplate(tmpl, noh)
    except Exception:
        return None
    to_pose = fixed.GetSubstructMatches(tmpl, uniquify=False, useChirality=False, maxMatches=200)
    to_ref = ref_mol.GetSubstructMatches(tmpl, uniquify=False, useChirality=False, maxMatches=200)
    if not to_pose or not to_ref:
        return None
    # One template->pose embedding is enough; varying the template->reference side covers every
    # distinct pose-to-reference correspondence up to the automorphism group.
    maps = [list(zip(to_pose[0], r)) for r in to_ref]
    return maps, fixed.GetConformer().GetPositions()


def score_pose(ref_prot, ref_mol, mod_prot, pose_mol, tmpl):
    """BiSyRMSD and LDDT-PLI of one pose against the held-out crystal complex."""
    got = correspondence(tmpl, ref_mol, pose_mol)
    if got is None:
        return None, None
    maps, PL_all = got
    RL = ref_mol.GetConformer().GetPositions()
    heavy = sorted({r for _, r in maps[0]})
    # --- superpose the model receptor's binding site onto the reference's, by residue number ---
    site = [a for a in ref_prot if a["name"] == "CA"
            and np.min(np.linalg.norm(RL[heavy] - a["xyz"], axis=1)) <= SITE_RADIUS]
    mod_ca = {a["seq"]: a["xyz"] for a in mod_prot if a["name"] == "CA"}
    pairs = [(a["xyz"], mod_ca[a["seq"]]) for a in site if a["seq"] in mod_ca]
    if len(pairs) < 3:
        return None, None
    Rm, t = kabsch(np.array([p[1] for p in pairs]), np.array([p[0] for p in pairs]))
    best_r, best_p = None, None
    # --- LDDT-PLI needs no superposition, but the model protein must be in the reference frame ---
    mod_map = {(a["seq"], a["name"]): a["xyz"] for a in mod_prot}
    keep = [j for j, a in enumerate(ref_prot) if (a["seq"], a["name"]) in mod_map]
    if len(keep) < 20:
        return None, None
    MP = np.array([Rm @ mod_map[(ref_prot[j]["seq"], ref_prot[j]["name"])] + t for j in keep])
    RP = np.array([ref_prot[j]["xyz"] for j in keep])
    Dref = np.linalg.norm(RL[heavy][:, None, :] - RP[None, :, :], axis=2)
    mask = Dref <= INCLUSION
    for mp in maps:
        order = dict(mp)                      # ref index -> pose index
        PL = np.array([Rm @ PL_all[order[i]] + t for i in heavy])
        r = float(np.sqrt(np.mean(np.sum((RL[heavy] - PL) ** 2, axis=1))))
        Dmod = np.linalg.norm(PL[:, None, :] - MP[None, :, :], axis=2)
        dif = np.abs(Dmod[mask] - Dref[mask])
        p = float(np.mean([(dif <= th).mean() for th in THRESH])) if mask.any() else 0.0
        if best_r is None or r < best_r:
            best_r = r
        if best_p is None or p > best_p:
            best_p = p
    return best_r, best_p


def contact_set(prot, pos, cutoff=4.5):
    """Residue-level contacts of a pose: the set of residue numbers within cutoff of any atom."""
    P = np.array([a["xyz"] for a in prot])
    seq = [a["seq"] for a in prot]
    d = np.linalg.norm(pos[:, None, :] - P[None, :, :], axis=2).min(0)
    return {seq[j] for j in range(len(seq)) if d[j] <= cutoff}


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="MYT,PK9")
    ap.add_argument("--modes", type=int, default=9)
    ap.add_argument("--out", default=RES + "logs/k112_ensdock.json")
    a = ap.parse_args(argv)
    targets = [t.strip().upper() for t in a.targets.split(",") if t.strip()]

    smi = comp_smiles(sorted(set(list(ENSEMBLE) + targets)))
    missing = [t for t in targets if not smi.get(t)]
    if missing:
        print("нет SMILES для:", missing)
        return 2

    work = _pl.Path(tempfile.mkdtemp(prefix="k112_"))
    print(f"рабочий каталог: {work}\n")
    results = {}

    for tgt in targets:
        own = set(ENSEMBLE.get(tgt, []))
        pool = [(l, p) for l, ids in ENSEMBLE.items() if l not in EXCLUDE_LIG
                for p in ids if p not in own]
        print(f"=== цель {tgt}: рецепторов {len(pool)} (свои записи {sorted(own)} ИСКЛЮЧЕНЫ) ===")
        lig_sdf = work / f"{tgt}.sdf"
        if not embed3d(smi[tgt], lig_sdf):
            print("  3D не построилась, пропуск")
            continue
        tmpl = Chem.MolFromSmiles(smi[tgt])
        # held-out truth: the target's own crystal complex
        ref_id = sorted(own)[0]
        rp, rl, _ = split_complex(parse(fetch(ref_id)), tgt)
        ref_mol = lig_mol(rl, smi[tgt])
        if ref_mol is None:
            print("  эталонный лиганд не собрался, пропуск")
            continue

        cands, unscored = [], 0
        for lig_id, pdb_id in pool:
            rf = receptor_files(pdb_id, lig_id, work)
            if rf is None:
                continue
            rec, box, prot, _, _ = rf
            out = work / f"{tgt}_{pdb_id}.sdf"
            if not dock(rec, box, lig_sdf, out, a.modes):
                print(f"  {pdb_id}: smina не отработала")
                continue
            for k, (m, aff) in enumerate(poses_from(out)):
                r, p = score_pose(rp, ref_mol, prot, m, tmpl)
                if r is None:
                    unscored += 1
                    continue
                pos = np.array([m.GetConformer().GetPositions()[i]
                                for i, at in enumerate(m.GetAtoms()) if at.GetSymbol() != "H"])
                cands.append({"rec": pdb_id, "k": k, "aff": aff, "rmsd": r, "pli": p,
                              "contacts": contact_set(prot, pos)})
            print(f"  {pdb_id}: поз {len(poses_from(out))}")

        print(f"  поз оценено {len(cands)}, НЕ оценено {unscored}")
        if unscored and not cands:
            print("  ВСЕ позы остались неоценёнными --- это дефект сопоставления, не результат\n")
            continue
        if not cands:
            print("  ни одной позы\n")
            continue

        # --- the three selection rules, on one pose set ---
        freq = {}
        for c in cands:
            for s in c["contacts"]:
                freq[s] = freq.get(s, 0) + 1
        n = len(cands)
        for c in cands:
            c["cons"] = (sum(freq[s] for s in c["contacts"]) / (len(c["contacts"]) * n)
                         if c["contacts"] else 0.0)
        by_aff = min(cands, key=lambda c: c["aff"] if c["aff"] == c["aff"] else 1e9)
        by_cons = max(cands, key=lambda c: c["cons"])
        oracle_r = min(cands, key=lambda c: c["rmsd"])
        oracle_p = max(cands, key=lambda c: c["pli"])
        print(f"\n  поз всего {n}, рецепторов {len({c['rec'] for c in cands})}")
        print("  %-26s %-7s %9s %9s" % ("правило отбора", "рецептор", "RMSD", "LDDT-PLI"))
        for nm, c in (("по скору smina", by_aff), ("ПО СОГЛАСИЮ КОНТАКТОВ", by_cons),
                      ("оракул (лучший RMSD)", oracle_r), ("оракул (лучший PLI)", oracle_p)):
            print("  %-26s %-7s %9.3f %9.4f" % (nm, c["rec"], c["rmsd"], c["pli"]))
        print("  для сравнения: потолок пункта 333 --- RMSD 1.454, PLI 0.7515;"
              " лидер доски PLI 0.4252\n")
        results[tgt] = {"n_poses": n, "по_скору": {k: by_aff[k] for k in ("rec", "rmsd", "pli")},
                        "по_согласию": {k: by_cons[k] for k in ("rec", "rmsd", "pli")},
                        "оракул_rmsd": {k: oracle_r[k] for k in ("rec", "rmsd", "pli")},
                        "оракул_pli": {k: oracle_p[k] for k in ("rec", "rmsd", "pli")}}

    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump(results, open(a.out, "w"), ensure_ascii=False, indent=1, default=float)
    print(f"сохранено: {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
