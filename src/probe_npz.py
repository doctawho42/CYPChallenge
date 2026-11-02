"""Re-key an embedding file into the two-key contract `src/submit.py:_probe_load` reads.

Why this exists rather than a flag on the writer. `_probe_load` takes the array from the key
`train` and the row identities from `molecule_name`, and it realigns against the table instead of
trusting either -- including a shifted negative control, because a comparison that cannot fail is
not a control. `src/embed_cypmodel.py` writes `pred`, `emb`, `ffn` and `smiles`, since one run of a
multitask model produces three candidate feature sets and the choice between them is a measurement
(item 362's device, `verify/k139_oadmetcyp.py`), not something to fix at write time. This script
makes whichever one won loadable by the submission builder, and re-runs the order control on the
way through so a re-key cannot silently desync the rows.

    uv run python src/probe_npz.py --src data/emb_oadmetcyp_train.npz --key emb \\
        --table data/rows.csv --out data/emb_oadmetcyp_probe.npz

Row identity. `rows.csv` and the blinded test file both carry `Molecule_Name` and `SMILES`; the
embedding files carry `smiles` only, so the match is on SMILES and the names are taken from the
table in that same order. The control is that the embedding's own SMILES column must equal the
table's exactly and must NOT equal it rolled by one.
Exit: 0 ran, 2 a control failed, 3 an input is missing.
"""
import argparse
import pathlib
import sys

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--src", required=True)
    ap.add_argument("--key", default="emb", help="emb | ffn | pred, or any array in the file")
    ap.add_argument("--table", required=True, help="rows.csv or the blinded test csv")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    for f in (a.src, a.table):
        if not pathlib.Path(f).exists():
            print("нет %s" % f, file=sys.stderr)
            return 3
    z = np.load(a.src, allow_pickle=True)
    if a.key not in z:
        print("в %s нет ключа %r; есть %r" % (a.src, a.key, list(z)), file=sys.stderr)
        return 3
    E = np.asarray(z[a.key], np.float64)
    smi = [str(x) for x in z["smiles"]]
    t = pd.read_csv(a.table)
    if len(t) != len(E):
        print("КОНТРОЛЬ: %d строк в %s против %d в таблице" % (len(E), a.src, len(t)),
              file=sys.stderr)
        return 2
    ok = smi == list(t.SMILES)
    rolled = (smi[-1:] + smi[:-1]) == list(t.SMILES)
    if not ok or rolled:
        print("КОНТРОЛЬ ПРОВАЛЕН: порядок (совпал=%s, сдвинутый=%s)" % (ok, rolled),
              file=sys.stderr)
        return 2
    if not np.isfinite(E).all():
        print("КОНТРОЛЬ: не-конечные значения в %s[%r]" % (a.src, a.key), file=sys.stderr)
        return 2
    np.savez_compressed(a.out, train=E,
                        molecule_name=t.Molecule_Name.to_numpy(object))
    print("%s[%s] %s -> %s; порядок совпал: %s; сдвинутый контроль (ожидается False): %s"
          % (pathlib.Path(a.src).name, a.key, E.shape, a.out, ok, rolled))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
