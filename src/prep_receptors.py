"""Build the four docking receptors from the PDB, by COLUMNS rather than whitespace.

Why this file exists. Item 298 documented the recipe in prose -- chain A plus its own heme,
waters, glycerol, DMSO and ions stripped -- but nothing performed it reproducibly, and
`data/receptors/` is gitignored, so the inputs the docking campaign runs on existed on
exactly one machine with no way to rebuild them. That became blocking the moment a second
person offered to run docking. Closing that hole is the whole of this file's value.

**A false alarm is recorded here because it nearly cost the campaign, and it was mine.**
Counting heme atoms across the four prepared receptors with `awk '$1=="HETATM"'` returned
129, which is 43 x 3 rather than 43 x 4, and I read that as 4WNV carrying no heme -- serious,
since 4WNV is CYP2D6 and a cavity missing its iron is a void. It was a parsing bug in the
diagnostic, not a defect in the data. PDB is fixed-column: columns 1-6 are the record name
and 7-11 the right-justified serial, so a five-digit serial gives "HETATM14450" with no
space and `line.split()[0]` is no longer "HETATM". 4WNV has 14449 ATOM records across
chains A-D, so all 615 of its HETATM lines are invisible to a whitespace split, while 2HI4
(3845), 1R9O (3650) and 3NXU (7356) stay under 10000 and parse fine. The heme was present
in all four the whole time: each old file's length equals its coordinate count plus one
terminator (3889, 3694, 3722, 3709), which only holds if 43 heme atoms were in it.

Two lessons, both worth more than the recipe below. **Never `split()` a PDB line** -- every
field here is sliced by column index. And **verify before writing**: the check that would
have caught this ("the three intact files must come out atom-identical") was stated and then
skipped in favour of one convenient call, which overwrote the originals in a gitignored
directory with no backup, so the retraction had to be argued from line counts instead of a
direct comparison. Hence `--dry-run` below, which is the default-safe way to use this.

Idempotent: writes only when the content actually changes, and does so by atomic rename, so
a docking run holding the old file open is unaffected (its descriptor keeps the old inode).
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D

import argparse
import os
import urllib.request

# enzyme -> PDB id. Chain A in every case; item 298 verified that all four co-crystal
# ligands coincide with chain A at 0.00 A while the other copies sit 50-89 A away.
CAVITY = {
    "CYP1A2": "2HI4",
    "CYP2C9": "1R9O",
    "CYP2D6": "4WNV",
    "CYP3A4": "3NXU",
}
KEEP_HET = {"HEM"}
N_HEM_ATOMS = 43          # identical in all four entries, checked
RAW = D + "receptors/raw/"
OUT = D + "receptors/"


def fetch(pdb_id):
    """Cache the raw entry under data/receptors/raw/ (the whole tree is gitignored)."""
    _pl.Path(RAW).mkdir(parents=True, exist_ok=True)
    path = RAW + pdb_id + ".pdb"
    if not os.path.exists(path) or os.path.getsize(path) < 10000:
        url = "https://files.rcsb.org/download/%s.pdb" % pdb_id
        with urllib.request.urlopen(url, timeout=60) as r, open(path, "wb") as f:
            f.write(r.read())
    return path


def prepare(pdb_id):
    """Chain A protein plus chain A heme, as raw PDB lines. Columns, not fields."""
    kept, stats = [], {"atom": 0, "hem": 0, "fe": 0, "dropped_het": {}}
    for ln in open(fetch(pdb_id), errors="ignore"):
        rec = ln[0:6]
        if rec not in ("ATOM  ", "HETATM"):
            continue
        res = ln[17:20].strip()
        chain = ln[21]
        if rec == "ATOM  ":
            if chain != "A":
                continue
            kept.append(ln)
            stats["atom"] += 1
        else:
            if chain != "A" or res not in KEEP_HET:
                stats["dropped_het"][res] = stats["dropped_het"].get(res, 0) + 1
                continue
            kept.append(ln)
            stats["hem"] += 1
            if ln[76:78].strip().upper() == "FE":
                stats["fe"] += 1

    # Assertions, because the defect this file fixes produced a receptor that looked fine.
    assert stats["atom"] > 3000, "%s: only %d chain-A protein atoms" % (pdb_id, stats["atom"])
    assert stats["hem"] == N_HEM_ATOMS, \
        "%s: %d heme atoms, expected %d -- the parse dropped the heme again" \
        % (pdb_id, stats["hem"], N_HEM_ATOMS)
    assert stats["fe"] == 1, "%s: %d Fe in the heme, expected 1" % (pdb_id, stats["fe"])
    chains = {ln[21] for ln in kept}
    assert chains == {"A"}, "%s: chains %s leaked through" % (pdb_id, sorted(chains))
    return kept, stats


def write_if_changed(path, lines, dry=False):
    """Atomic replace, and only when the bytes differ. `dry` compares and writes nothing.

    Atomic because a docking run may hold the old file open: `os.replace` keeps that
    descriptor pointing at the old inode, so a chunk in flight cannot read a half-written
    receptor.
    """
    body = "".join(lines) + "TER\nEND\n"
    if os.path.exists(path):
        with open(path) as f:
            if f.read() == body:
                return "без изменений"
        if dry:
            return "ИЗМЕНИЛСЯ БЫ (сухой прогон, ничего не записано)"
    elif dry:
        return "БЫЛ БЫ СОЗДАН (сухой прогон)"
    tmp = path + ".tmp"
    with open(tmp, "w") as f:
        f.write(body)
    os.replace(tmp, path)
    return "ПЕРЕЗАПИСАН"


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--dry-run", action="store_true",
                    help="сравнить и НЕ записывать; так и надо смотреть первый раз")
    a = ap.parse_args()

    print("Сборка рецепторов докинга из PDB, разбор ПО КОЛОНКАМ")
    print("(зачем файл существует и какая ложная тревога тут записана -- в докстринге)")
    if a.dry_run:
        print("РЕЖИМ: сухой прогон, файлы не записываются")
    print()
    print("    %-8s %-6s %7s %6s %5s %-28s %s"
          % ("фермент", "PDB", "белок", "гем", "Fe", "выброшено HETATM", "файл"))
    for enz in sorted(CAVITY):
        pdb_id = CAVITY[enz]
        lines, st = prepare(pdb_id)
        path = OUT + pdb_id + "_rec.pdb"
        before = None
        if os.path.exists(path):
            with open(path) as f:
                before = sum(1 for x in f if x[0:6] == "HETATM")
        verdict = write_if_changed(path, lines, dry=a.dry_run)
        dropped = ", ".join("%s:%d" % kv for kv in sorted(st["dropped_het"].items()))
        print("    %-8s %-6s %7d %6d %5d %-28s %s"
              % (enz, pdb_id, st["atom"], st["hem"], st["fe"], dropped[:28], verdict))
        if before is not None and before != st["hem"]:
            print("        ^ было %d записей HETATM, стало %d" % (before, st["hem"]))

    print("\n    Проверки внутри prepare(): 43 атома гема, ровно один Fe, единственная")
    print("    цепь A, ни одного не-гемового HETATM. Дефект, который этот файл чинит,")
    print("    давал рецептор, выглядевший исправным, поэтому проверки в коде, а не в прозе.")


if __name__ == "__main__":
    main()
