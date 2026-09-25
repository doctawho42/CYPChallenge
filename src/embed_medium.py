"""Frozen chemprop_medium embeddings for every row of data/rows.csv.

Follows src/embed.py's contract exactly: this runs OUT OF PROCESS in a separate
interpreter, touches none of the repository's measured pins, and leaves a .npz on
disk. Nothing downstream imports chemprop.

    uv venv <scratch>/encenv --python 3.12
    uv pip install --python <scratch>/encenv/bin/python 'chemprop>=2.1'
    <scratch>/encenv/bin/python embed_medium.py --ckpt <...>/chemprop_medium.pt

Which route. The task described the extraction as
``model.agg(model.message_passing(bmg, V_d), bmg.batch)``. The committed competitor
code does NOT do that: cyp_submission/chemprop_transfer.py calls
``model.fingerprint(bmg, V_d, X_d)``, which is message_passing -> agg -> self.bn.
Those two differ whenever ``bn`` is not Identity, so BOTH are computed here and
compared, and the one their ridge probes actually saw (fingerprint) is what is
written out.

Unparseable SMILES are handled the way src/feats.py handles them -- RDKit parses
first, failures are collected and reported by name, never zero-filled. rows.csv is
already feats.py's surviving set, so the expected count is zero and a non-zero count
is a defect to report rather than to paper over.
"""
import argparse
import hashlib
import pathlib
import sys

import numpy as np
import pandas as pd
import torch
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


def load_model(ckpt):
    from chemprop.models.model import MPNN
    # The published checkpoint was saved on a CUDA device; without map_location
    # torch refuses to deserialise it on this machine.
    m = MPNN.load_from_file(str(ckpt), map_location="cpu")
    m.eval()
    return m


def embed(model, smiles, batch=256, route="fingerprint"):
    """Graph-level embedding per molecule, in the order given.

    drop_last=False is not a default worth trusting: the competitor's own file
    records that chemprop silently drops a size-1 remainder batch, which desyncs
    every index-aligned array downstream. It is passed explicitly, and the row
    count is asserted by the caller regardless.
    """
    from chemprop import data, featurizers

    feat = featurizers.SimpleMoleculeMolGraphFeaturizer()
    pts = [data.MoleculeDatapoint.from_smi(s, [0.0]) for s in smiles]
    dset = data.MoleculeDataset(pts, feat)
    loader = data.build_dataloader(dset, num_workers=0, batch_size=batch,
                                   shuffle=False, drop_last=False)
    out = []
    with torch.inference_mode():
        for b in loader:
            bmg, V_d, X_d, *_ = b
            if route == "fingerprint":
                h = model.fingerprint(bmg, V_d, X_d)
            elif route == "bare":
                h = model.agg(model.message_passing(bmg, V_d), bmg.batch)
            else:
                raise ValueError(route)
            out.append(h.numpy())
    E = np.concatenate(out, axis=0)
    if len(E) != len(smiles):
        raise SystemExit(f"embedded {len(E)} of {len(smiles)} molecules -- rows desynced")
    return E.astype(np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", required=True)
    ap.add_argument("--rows", required=True)
    ap.add_argument("--out", required=True)
    a = ap.parse_args()

    rows = pd.read_csv(a.rows)
    smi = list(rows.SMILES)
    print(f"rows.csv: {len(rows)} molecules", flush=True)

    # --- parse control, matching src/feats.py's rule -------------------------
    bad = [(i, rows.Molecule_Name[i]) for i, s in enumerate(smi)
           if Chem.MolFromSmiles(s) is None]
    print(f"RDKit-unparseable rows: {len(bad)}"
          + (f" -> {bad[:10]}" if bad else " (expected 0: rows.csv is feats.py's survivors)"),
          flush=True)
    if bad:
        raise SystemExit("unparseable rows present; feats.py would have dropped these, "
                         "so rows.csv and the checkpoint disagree -- stopping rather than zero-filling")

    model = load_model(a.ckpt)
    print(f"encoder: mp={type(model.message_passing).__name__} "
          f"output_dim={model.message_passing.output_dim} "
          f"agg={type(model.agg).__name__} bn={type(model.bn).__name__}", flush=True)

    E = embed(model, smi, route="fingerprint")
    print(f"embedding (fingerprint route): {E.shape} {E.dtype}", flush=True)

    Eb = embed(model, smi[:512], route="bare")
    same = np.allclose(Eb, E[:512], atol=1e-6)
    print(f"bare agg(message_passing) route vs fingerprint on first 512: "
          f"{'IDENTICAL' if same else 'DIFFERENT'} "
          f"(max abs diff {np.abs(Eb - E[:512]).max():.3e})", flush=True)

    # --- control 1: determinism on a repeated SMILES -------------------------
    probe = [smi[0], smi[1], smi[0], smi[2], smi[1]]
    P = embed(model, probe, route="fingerprint")
    id_ok = np.array_equal(P[0], P[2]) and np.array_equal(P[1], P[4])
    diff_ok = not np.array_equal(P[0], P[1]) and not np.array_equal(P[0], P[3])
    print(f"control identical-SMILES-twice: {'PASS' if id_ok else 'FAIL'} "
          f"(row0 vs row2 max diff {np.abs(P[0]-P[2]).max():.3e})", flush=True)
    print(f"control different-molecules-differ: {'PASS' if diff_ok else 'FAIL'} "
          f"(row0 vs row1 max diff {np.abs(P[0]-P[1]).max():.3e})", flush=True)

    # --- control 2: order is preserved, not merely assumed -------------------
    # Embed a PERMUTED copy and check it equals the permuted embedding. A loader
    # that reordered or dropped a batch fails here; a loader that merely "looked
    # fine" cannot pass this by accident.
    rng = np.random.default_rng(0)
    perm = rng.permutation(len(smi))[:1024]
    Ep = embed(model, [smi[i] for i in perm], route="fingerprint")
    ord_ok = np.allclose(Ep, E[perm], atol=1e-5)
    print(f"control order-preserved under permutation (n=1024): "
          f"{'PASS' if ord_ok else 'FAIL'} (max abs diff {np.abs(Ep - E[perm]).max():.3e})",
          flush=True)

    # --- control 3: the unparseable path actually raises ---------------------
    try:
        embed(model, ["this-is-not-a-smiles"], route="fingerprint")
        raised = False
    except Exception as e:
        raised = True
        why = type(e).__name__
    print(f"control unparseable-SMILES-raises: {'PASS (' + why + ')' if raised else 'FAIL (returned silently)'}",
          flush=True)

    if not (id_ok and diff_ok and ord_ok and raised):
        raise SystemExit("a control failed; not writing the array")

    np.savez_compressed(a.out, train=E,
                        molecule_name=np.array(rows.Molecule_Name, dtype=object))
    h = hashlib.sha256(pathlib.Path(a.out).read_bytes()).hexdigest()
    print(f"\nwrote {a.out}")
    print(f"  shape {E.shape} dtype {E.dtype} sha256 {h}")
    print(f"  array sha256 (raw bytes) {hashlib.sha256(E.tobytes()).hexdigest()}")
    print(f"  column spread: median sd {np.median(E.std(0)):.4f}, "
          f"dead columns {int((E.std(0) < 1e-6).sum())} of {E.shape[1]}")


if __name__ == "__main__":
    main()
