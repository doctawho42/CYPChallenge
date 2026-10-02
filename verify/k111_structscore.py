"""Local LDDT-PLI and BiSyRMSD, and the intrinsic ceiling of the CYP3A4 structure track.

WHY A LOCAL SCORER AT ALL. The structure track gives twenty blinded targets and one board
reading per upload every twelve hours. Nothing can be chosen that way. But the PDB holds 129
experimental CYP3A4 entries (UniProt P08684), 119 of them with a drug-like ligand, so the same
question -- does this pipeline place a ligand correctly in this protein -- has 119 answers already
in hand. This file implements the two metrics the organisers score on, so every design choice can
be measured before an upload is spent.

WHAT IT DOES NOT CLAIM, said first. OpenStructure is not installed here (conda or docker), so the
ABSOLUTE agreement of these numbers with the organisers' `SCRMSDScorer` and `LDDTPLIScorer` is
UNVERIFIED and this file does not assert it. What is implemented is each metric's stated definition,
with controls that fix its behaviour at the ends. That is enough for its purpose -- RANKING our own
options against each other on real complexes -- and not enough to predict a board value, which is
why no board value is predicted from it anywhere.

    BiSyRMSD  superpose on the reference binding site (CA within 8 A of any reference ligand atom,
              matched by residue), then ligand RMSD minimised over the ligand's graph automorphisms.
    LDDT-PLI  no superposition. Over protein-ligand atom pairs whose REFERENCE distance is inside
              the inclusion radius, the fraction of pairs whose distance is preserved, averaged over
              the four standard tolerances 0.5, 1, 2, 4 A, maximised over automorphisms.

THE CEILING, which is this file's own measurement rather than an instrument reading. Twelve ligands
appear in more than one CYP3A4 entry -- STR in three, RIT in four, and nine others in two. Scoring
those pairs against each other asks what the SAME ligand does to the SAME protein in two independent
experiments. Whatever that number is, no prediction method can beat it: it is the target's own
irreproducibility plus crystallographic error. The organisers say CYP3A4 is notoriously hard because
similar ligands adopt drastically different poses; this measures the stronger version of that claim,
for identical ligands, and nobody has published it for this target.

CONTROLS, and the file refuses its own output without them. A reference scored against itself must
give BiSyRMSD 0 and LDDT-PLI 1. A reference scored against copies of itself displaced by a growing
offset must degrade MONOTONICALLY in both metrics. A scorer that cannot fail at the ends cannot be
trusted in the middle, and "all pairs look reasonable" is this project's most frequent defect.

    uv run python verify/k111_structscore.py --controls
    uv run python verify/k111_structscore.py --ceiling

Writes results/logs/k111_structscore.json. PDB files are cached under data/pdb3a4/.
Exit: 0 controls passed, 1 a control failed (any measurement printed is void), 2 inputs missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES

import argparse
import itertools
import json
import urllib.request

import numpy as np

CACHE = _pl.Path(D) / "pdb3a4"
THRESH = (0.5, 1.0, 2.0, 4.0)
INCLUSION = 6.0          # contact radius for PLI pairs
SITE_RADIUS = 8.0        # binding-site CA radius, per the organisers' notebook
SALT = {"HOH", "GOL", "DMS", "EDO", "PEG", "SO4", "PO4", "CL", "NA", "K", "MG",
        "CA", "ZN", "FE", "ACT", "TRS", "MES", "IMD", "NO3", "1PE", "P6G", "PGE",
        "MPD", "FMT", "CIT", "BME", "EPE", "NH4", "BR", "IOD", "PG4"}


def fetch(pdb_id):
    """One PDB file, cached. Same URL shape src/prep_receptors.py uses."""
    CACHE.mkdir(parents=True, exist_ok=True)
    path = CACHE / f"{pdb_id.lower()}.pdb"
    if not path.exists():
        url = f"https://files.rcsb.org/download/{pdb_id.upper()}.pdb"
        with urllib.request.urlopen(url, timeout=60) as r:
            path.write_bytes(r.read())
    return path


def parse(path):
    """Atoms by COLUMN, never by split() -- item 302: a five-digit serial leaves no space and
    `line.split()[0]` stops being 'HETATM', which is how a heme once went missing."""
    out = []
    for ln in _pl.Path(path).read_text().splitlines():
        rec = ln[0:6]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        if ln[16] not in (" ", "A"):          # altloc: keep the first only
            continue
        out.append({"rec": rec.strip(), "name": ln[12:16].strip(), "res": ln[17:20].strip(),
                    "chain": ln[21:22], "seq": ln[22:27].strip(),
                    "xyz": np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]),
                    "el": (ln[76:78].strip() or ln[12:14].strip()).upper()})
    return out


def split_complex(atoms, lig_id):
    """Protein (+heme) and the CATALYTIC copy of the named ligand, from one chain.

    The copy is chosen as the one nearest the heme iron, and its distance is returned, because
    "the first copy by residue number" is not a site. CYP3A4 has a peripheral surface site as
    well as the catalytic one: caffeine appears in three and six copies in 8SO1 and 8SO2, from
    4.6 to 23.7 A from the iron, and all three steroid entries (1W0F, 5A1P, 5A1R) place STR at
    ~22 A, i.e. not in the catalytic pocket at all. Picking by residue number got the right site
    in those cases by luck; the first version of this file did exactly that and its ceiling
    numbers were inflated by the peripheral triplet.
    """
    fe = [a["xyz"] for a in atoms if a["res"] == "HEM" and a["el"] == "FE"]
    copies = {}
    for a in atoms:
        if a["res"] == lig_id and a["el"] != "H":
            copies.setdefault((a["chain"], a["seq"]), []).append(a)
    if not copies:
        return None, None, None
    def dist(v):
        c = np.mean([x["xyz"] for x in v], axis=0)
        return min(float(np.linalg.norm(c - f)) for f in fe) if fe else float("inf")
    key = min(copies, key=lambda k: dist(copies[k]))
    lig, d_fe = copies[key], dist(copies[key])
    ch = key[0]
    prot = [a for a in atoms if a["chain"] == ch and a["res"] != lig_id
            and a["res"] not in SALT and a["el"] != "H"]
    return prot, lig, d_fe


def automorphisms(lig):
    """Index permutations of the ligand that preserve element and intramolecular distances.

    Built from geometry rather than from a chemical graph on purpose: both sides of every
    comparison here are the SAME chem-comp, so atom NAMES already correspond, and the only
    symmetry that matters is the one that leaves the distance matrix invariant. Capped, because
    a large symmetric ligand can have many and the cap is reported rather than silent.
    """
    n = len(lig)
    el = [a["el"] for a in lig]
    Dm = np.linalg.norm(np.array([a["xyz"] for a in lig])[:, None, :]
                        - np.array([a["xyz"] for a in lig])[None, :, :], axis=2)
    # candidate images per atom: same element and same sorted distance profile
    prof = [tuple(np.round(np.sort(Dm[i]), 2)) for i in range(n)]
    cand = [[j for j in range(n) if el[j] == el[i] and prof[j] == prof[i]] for i in range(n)]
    if all(len(c) == 1 for c in cand):
        return [list(range(n))], False
    perms, capped = [], False
    groups = {}
    for i, c in enumerate(cand):
        groups.setdefault(tuple(c), []).append(i)
    # permute within each equivalence class, product across classes
    blocks = []
    for key, idxs in groups.items():
        if len(idxs) > 6:
            capped = True
            blocks.append([tuple(idxs)])
        else:
            blocks.append(list(itertools.permutations(idxs)))
    for combo in itertools.islice(itertools.product(*blocks), 5000):
        p = list(range(n))
        for key_idxs, perm in zip(groups.values(), combo):
            for src, dst in zip(key_idxs, perm):
                p[src] = dst
        if len(set(p)) == n:
            perms.append(p)
    if not perms:
        perms = [list(range(n))]
    return perms, capped


def kabsch(P, Q):
    """Rotation+translation taking P onto Q."""
    pc, qc = P.mean(0), Q.mean(0)
    H = (P - pc).T @ (Q - qc)
    U, _, Vt = np.linalg.svd(H)
    d = np.sign(np.linalg.det(Vt.T @ U.T))
    Rm = Vt.T @ np.diag([1, 1, d]) @ U.T
    return Rm, qc - Rm @ pc


def bisyrmsd(ref_prot, ref_lig, mod_prot, mod_lig):
    """Binding-site superposed, symmetry-corrected ligand RMSD."""
    L = np.array([a["xyz"] for a in ref_lig])
    site = [a for a in ref_prot if a["name"] == "CA"
            and np.min(np.linalg.norm(L - a["xyz"], axis=1)) <= SITE_RADIUS]
    # Keyed by residue number and atom name ONLY. Chain letters are arbitrary between PDB
    # entries -- 38FV calls its copies A/B/C and 3NXU calls its A/B -- and keying on the chain
    # silently emptied the intersection, which raised deep inside numpy and dropped four of
    # nineteen pairs without a word. split_complex already restricts each side to one chain,
    # so (seq, name) is unambiguous here.
    mod_ca = {a["seq"]: a["xyz"] for a in mod_prot if a["name"] == "CA"}
    pairs = [(a["xyz"], mod_ca[a["seq"]]) for a in site if a["seq"] in mod_ca]
    if len(pairs) < 3:
        return None, len(pairs)
    Rm, t = kabsch(np.array([p[1] for p in pairs]), np.array([p[0] for p in pairs]))
    M = np.array([Rm @ a["xyz"] + t for a in mod_lig])
    names = {a["name"]: i for i, a in enumerate(mod_lig)}
    order = [names.get(a["name"]) for a in ref_lig]
    if any(o is None for o in order):
        return None, len(pairs)
    M = M[order]
    perms, _ = automorphisms(ref_lig)
    best = min(float(np.sqrt(np.mean(np.sum((L - M[p]) ** 2, axis=1)))) for p in perms)
    return best, len(pairs)


def lddt_pli(ref_prot, ref_lig, mod_prot, mod_lig):
    """Contact preservation, no superposition."""
    names = {a["name"]: i for i, a in enumerate(mod_lig)}
    order = [names.get(a["name"]) for a in ref_lig]
    if any(o is None for o in order):
        return None, 0
    RP = np.array([a["xyz"] for a in ref_prot])
    RL = np.array([a["xyz"] for a in ref_lig])
    Dref = np.linalg.norm(RL[:, None, :] - RP[None, :, :], axis=2)
    mask = Dref <= INCLUSION
    if not mask.any():
        return None, 0
    mod_map = {(a["seq"], a["name"]): a["xyz"] for a in mod_prot}
    keep = [j for j, a in enumerate(ref_prot) if (a["seq"], a["name"]) in mod_map]
    if len(keep) < 20:
        return None, len(keep)          # loudly, rather than failing inside numpy
    MP = np.array([mod_map[(ref_prot[j]["seq"], ref_prot[j]["name"])] for j in keep])
    Dref = Dref[:, keep]
    mask = mask[:, keep]
    ML_all = np.array([a["xyz"] for a in mod_lig])[order]
    perms, _ = automorphisms(ref_lig)
    best = 0.0
    for p in perms:
        Dmod = np.linalg.norm(ML_all[p][:, None, :] - MP[None, :, :], axis=2)
        dif = np.abs(Dmod[mask] - Dref[mask])
        sc = float(np.mean([(dif <= t).mean() for t in THRESH]))
        best = max(best, sc)
    return best, int(mask.sum())


def score_pair(ref_id, mod_id, lig_id):
    rp, rl, dr = split_complex(parse(fetch(ref_id)), lig_id)
    mp, ml, dm = split_complex(parse(fetch(mod_id)), lig_id)
    if not rl or not ml:
        return None
    # A pair is only about the catalytic pocket if BOTH copies are in it. Peripheral copies and
    # mismatched sites are reported and excluded rather than averaged in.
    site = "каталитический" if max(dr, dm) <= 15.0 else "ПЕРИФЕРИЧЕСКИЙ"
    if abs(dr - dm) > 5.0:
        site = "РАЗНЫЕ САЙТЫ"
    r, npair = bisyrmsd(rp, rl, mp, ml)
    p, ncon = lddt_pli(rp, rl, mp, ml)
    return {"ref": ref_id, "mod": mod_id, "lig": lig_id, "bisyrmsd": r, "lddt_pli": p,
            "site_ca": npair, "contacts": ncon, "lig_atoms": len(rl),
            "fe_ref": round(dr, 1), "fe_mod": round(dm, 1), "сайт": site}


def controls():
    """Both metrics must be exact at zero displacement and monotone away from it."""
    ok = True
    print("КОНТРОЛЬ 1 --- структура против СЕБЯ: RMSD 0, LDDT-PLI 1\n")
    for pid, lig in (("3NXU", "RIT"), ("6MA7", "TPF"), ("8SO1", "CFF")):
        s = score_pair(pid, pid, lig)
        good = s and abs(s["bisyrmsd"]) < 1e-9 and abs(s["lddt_pli"] - 1.0) < 1e-9
        ok &= bool(good)
        print("  %-5s %-5s BiSyRMSD %.2e  LDDT-PLI %.6f  контактов %5d  %s" % (
            pid, lig, s["bisyrmsd"], s["lddt_pli"], s["contacts"], "ок" if good else "ПРОВАЛ"))

    print("\nКОНТРОЛЬ 2 --- сдвиг лиганда: обе метрики обязаны деградировать МОНОТОННО\n")
    rp, rl, _ = split_complex(parse(fetch("3NXU")), "RIT")
    prev_r, prev_p, mono = -1.0, 2.0, True
    print("  сдвиг, A   BiSyRMSD   LDDT-PLI")
    for off in (0.0, 0.25, 0.5, 1.0, 2.0, 4.0):
        ml = [dict(a, xyz=a["xyz"] + np.array([off, 0.0, 0.0])) for a in rl]
        r, _ = bisyrmsd(rp, rl, rp, ml)
        p, _ = lddt_pli(rp, rl, rp, ml)
        bad = (r < prev_r - 1e-9) or (p > prev_p + 1e-9)
        mono &= not bad
        print("  %8.2f %10.4f %10.4f%s" % (off, r, p, "   <-НЕ МОНОТОННО" if bad else ""))
        prev_r, prev_p = r, p
    ok &= mono
    print("\n  монотонность:", "ок" if mono else "ПРОВАЛ")
    print("\nКОНТРОЛЬ 3 --- симметрия: перестановка эквивалентных атомов не должна менять счёт")
    perms, capped = automorphisms(rl)
    print("  автоморфизмов у RIT: %d%s" % (len(perms), " (ограничено)" if capped else ""))
    return ok


# ligands appearing in more than one CYP3A4 entry, from the P08684 PDB listing
MULTI = {"STR": ["1W0F", "5A1P", "5A1R"], "RIT": ["3NXU", "5VC0", "5VCE", "38FV"],
         "MYT": ["1W0G", "6MA6"], "PK9": ["4D6Z", "4D75"], "08Y": ["3UA1", "5VCG"],
         "FZV": ["6DAC", "6DAJ"], "G0D": ["6DA2", "6DA3"], "CFF": ["8SO1", "8SO2"],
         "MWY": ["38FX", "6OOA"], "WZN": ["8EWD", "8EWR"], "X2Q": ["8EWP", "8EWS"],
         "KLN": ["2V0M", "38FW"]}


def ceiling():
    print("ВНУТРЕННИЙ ПОТОЛОК --- один лиганд, два независимых кристалла\n")
    print("%-6s %-6s %-6s %10s %10s %7s %6s %7s %s" % (
        "лиганд", "эталон", "модель", "BiSyRMSD", "LDDT-PLI", "контакт", "атомов", "Fe A", "сайт"))
    rows, dropped = [], []
    for lig, ids in sorted(MULTI.items()):
        for a, b in itertools.combinations(ids, 2):
            try:
                s = score_pair(a, b, lig)
            except Exception as e:
                print("  %-6s %-6s %-6s  ОШИБКА: %s" % (lig, a, b, str(e)[:60]))
                dropped.append((lig, a, b, str(e)[:60]))
                continue
            if not s or s["bisyrmsd"] is None or s["lddt_pli"] is None:
                print("  %-6s %-6s %-6s  НЕ СОПОСТАВИЛОСЬ (мало общих атомов)" % (lig, a, b))
                dropped.append((lig, a, b, "мало общих атомов"))
                continue
            print("%-6s %-6s %-6s %10.3f %10.4f %7d %6d %3.1f/%3.1f %s" % (
                lig, a, b, s["bisyrmsd"], s["lddt_pli"], s["contacts"], s["lig_atoms"],
                s["fe_ref"], s["fe_mod"], s["сайт"]))
            if s["сайт"] == "каталитический":
                rows.append(s)
    if rows:
        r = np.array([x["bisyrmsd"] for x in rows]); p = np.array([x["lddt_pli"] for x in rows])
        print("\n  пар в КАТАЛИТИЧЕСКОМ сайте (только они идут в сводку): %d" % len(rows))
        print("  BiSyRMSD  медиана %.3f  среднее %.3f  разброс %.3f..%.3f" % (
            np.median(r), r.mean(), r.min(), r.max()))
        print("  LDDT-PLI  медиана %.4f  среднее %.4f  разброс %.4f..%.4f" % (
            np.median(p), p.mean(), p.min(), p.max()))
    print("\n  КОНТРОЛЬ ПОЛНОТЫ: пар потеряно %d%s" % (
        len(dropped), "" if not dropped else " --- " + "; ".join(
            "%s %s/%s (%s)" % d for d in dropped)))
    if dropped:
        print("  потерянная пара --- это не ноль, а неизмеренное; сводка выше НЕПОЛНА")
    return rows


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--controls", action="store_true")
    ap.add_argument("--ceiling", action="store_true")
    ap.add_argument("--out", default=RES + "logs/k111_structscore.json")
    a = ap.parse_args(argv)
    if not (a.controls or a.ceiling):
        a.controls = a.ceiling = True

    # Merge into whatever is already logged. A --controls-only run used to overwrite the file
    # and silently erase the ceiling measurement; a log a partial run can truncate is not a log.
    out = {}
    if _pl.Path(a.out).exists():
        try:
            out = json.load(open(a.out))
        except Exception:
            out = {}
    if a.controls:
        ok = controls()
        out["контроли_прошли"] = bool(ok)
        if not ok:
            print("\nКОНТРОЛЬ ПРОВАЛЕН --- любое измерение ниже НЕДЕЙСТВИТЕЛЬНО")
            return 1
        print()
    if a.ceiling:
        out["потолок"] = ceiling()
    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump(out, open(a.out, "w"), ensure_ascii=False, indent=1, default=float)
    print(f"\nсохранено: {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
