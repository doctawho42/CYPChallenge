"""Embed every molecule with the organisers' pretrained chemprop encoder.

Why this exists. The ablation grid compares six feature sets - FP, DESC, MECH and their
unions - and every one of them is Morgan counts plus RDKit descriptors plus the mechanistic
block. So when the document says the choice of model is worth 0.005 while the post-processing
is worth 0.051, that is true *inside one representation on one dataset*, and it was written as
if it were true in general. The search was local and got read as global.

Section 7 has specified a pretrained encoder since the first draft and names the weights the
organisers published. It was never built. This script builds the cheapest honest version of
it: take the pretrained encoder, embed the molecules, and hand the embedding to exactly the
learner and the folds the ablation already uses. One row of the same table, directly
comparable, with nothing else changed.

Which weights, and why this one. The organisers published sixteen chemprop pretraining
checkpoints, each a D-MPNN trained to predict a different descriptor family. `rdkit2d` is
chosen deliberately: it is pretrained to predict the very RDKit descriptors that make up our
DESC block, so the comparison is not "learned features against hand features" in general but
the sharper question - does a *learned* representation of that information beat handing the
information over directly? The CYP-fine-tuned baseline in the same organisation was NOT used,
because a model already trained on CYP pIC50 data would make the comparison circular.

Runs in a separate environment. chemprop pulls in forty-one dependencies and this repository's
pins are measured rather than cautious (see CLAUDE.md), so the encoder must not become a
runtime dependency of the pipeline. The contract is a file: this script writes
data/emb_chemprop_rdkit2d.npz and everything downstream just reads an array.

    uv venv /tmp/enc --python 3.12
    uv pip install --python /tmp/enc/bin/python chemprop
    curl -L -o /tmp/rdkit2d__s42.pt https://huggingface.co/openadmet/\\
chemprop-foundation-pretraining-weights/resolve/main/rdkit2d__s42.pt
    /tmp/enc/bin/python src/embed.py --ckpt /tmp/rdkit2d__s42.pt

Writes data/emb_chemprop_rdkit2d.npz with `train` (n_train, 2048) in rows.csv order and
`test` (750, 2048) in the blinded file's order.
"""
import argparse
import pathlib

import numpy as np
import pandas as pd
import torch

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = str(ROOT / "data") + "/"


def build_encoder(ckpt_path):
    """Reconstruct the message-passing block and load the published weights into it.

    The checkpoint carries its own hyper-parameters, so nothing here is guessed: d_v 72,
    d_e 14, d_h 2048, depth 6, LeakyReLU, no bias. Loading is strict - a silently partial
    load would produce embeddings that look reasonable and mean nothing.
    """
    from chemprop.nn import BondMessagePassing

    blob = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    hp = blob["hyper_parameters"]
    mp = BondMessagePassing(
        d_v=hp["d_v"], d_e=hp["d_e"], d_h=hp["d_h"], bias=hp["bias"],
        depth=hp["depth"], dropout=hp["dropout"], activation=hp["activation"],
        undirected=hp["undirected"],
    )
    missing, unexpected = mp.load_state_dict(blob["state_dict"], strict=False)
    if missing or unexpected:
        raise SystemExit(f"веса легли не полностью: пропущено {missing}, лишнее {unexpected}")
    mp.eval()
    return mp, hp


def embed(smiles, mp, batch=64):
    """Mean-pooled node representation per molecule, in the order given."""
    from chemprop.data import MoleculeDatapoint, MoleculeDataset, build_dataloader

    ds = MoleculeDataset([MoleculeDatapoint.from_smi(s) for s in smiles])
    dl = build_dataloader(ds, batch_size=batch, shuffle=False)
    out = []
    with torch.no_grad():
        for bmb in dl:
            H = mp(bmb.bmg)                      # (атомы всего батча, d_h)
            idx = torch.as_tensor(bmb.bmg.batch, dtype=torch.long)
            n = int(idx.max()) + 1
            s = torch.zeros(n, H.shape[1], dtype=H.dtype)
            s.index_add_(0, idx, H)
            cnt = torch.zeros(n, dtype=H.dtype).index_add_(
                0, idx, torch.ones(len(idx), dtype=H.dtype))
            out.append((s / cnt[:, None]).numpy())
    return np.vstack(out).astype(np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", required=True)
    ap.add_argument("--out", default=D + "emb_chemprop_rdkit2d.npz")
    ap.add_argument("--limit", type=int, default=0, help="только первые N, для проверки")
    a = ap.parse_args()

    mp, hp = build_encoder(a.ckpt)
    print(f"энкодер собран: d_h {hp['d_h']}, глубина {hp['depth']}, "
          f"{sum(p.numel() for p in mp.parameters()):,} параметров", flush=True)

    rows = pd.read_csv(D + "rows.csv")
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    tr_smi = list(rows.SMILES)
    te_smi = list(te.SMILES)
    if a.limit:
        tr_smi, te_smi = tr_smi[:a.limit], te_smi[:a.limit]

    print(f"обучение: {len(tr_smi)} молекул", flush=True)
    E_tr = embed(tr_smi, mp)
    print(f"  готово {E_tr.shape}", flush=True)
    print(f"тест: {len(te_smi)} молекул", flush=True)
    E_te = embed(te_smi, mp)
    print(f"  готово {E_te.shape}", flush=True)

    # Строки должны совпасть с rows.csv элемент в элемент, иначе признаки и метки
    # разъедутся молча - та же ловушка, о которой предупреждает CLAUDE.md.
    if not a.limit and (len(E_tr) != len(rows) or len(E_te) != len(te)):
        raise SystemExit("число строк не совпало с rows.csv или тестовым файлом")

    np.savez_compressed(a.out, train=E_tr, test=E_te)
    print(f"\nсохранено: {a.out}")
    print(f"  разброс по столбцам: медиана {np.median(E_tr.std(0)):.4f}, "
          f"мёртвых столбцов {int((E_tr.std(0) < 1e-6).sum())} из {E_tr.shape[1]}")


if __name__ == "__main__":
    main()
