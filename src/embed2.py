"""Embed every molecule with a pretrained SMILES transformer.

Why another encoder script. `src/embed.py` already runs a pretrained encoder, and item 68
read its result as closing the question. It does not, and the reason is in that script's own
docstring: it deliberately chose chemprop `rdkit2d`, a checkpoint pretrained to predict the
very RDKit descriptors that make up our DESC block. That is the sharpest possible choice for
the question "does a *learned* representation beat handing the same information over
directly", and the **least** favourable one for the question that actually matters here --
does a pretrained representation carry anything we do not already have. An encoder trained to
reproduce our own features has little room to add to them.

So this script takes the other kind: transformers pretrained on raw structure at scale, with
no descriptor supervision anywhere in their objective. What they know that our matrix does not
is whatever survives from a hundred million to a billion molecules of chemistry, and that is
information we have never put in front of the learner.

Three checkpoints, chosen to differ in pretraining objective rather than in size:

  ibm/MoLFormer-XL-both-10pct   linear-attention transformer, masked language modelling over
                                1.1 billion molecules from ZINC and PubChem. The largest
                                pretraining corpus available for SMILES.
  DeepChem/ChemBERTa-77M-MTR    RoBERTa over 77M PubChem SMILES, but pretrained by multitask
                                *regression* onto computed properties rather than by masking.
                                A different objective on a smaller corpus -- the contrast with
                                MoLFormer is the point, not the ranking.
  seyonec/PubChem10M_SMILES_BPE_450k
                                RoBERTa, masked language modelling over 10M PubChem SMILES.
                                The small-corpus control for MoLFormer: same objective, two
                                orders of magnitude less data.

Runs in a separate environment, exactly as `src/embed.py` does and for the same reason -- the
pins in `pyproject.toml` are measured and `transformers` must not become a runtime dependency
of the pipeline. The contract is a file; everything downstream reads an array.

    uv venv /tmp/enc2 --python 3.12
    uv pip install --python /tmp/enc2/bin/python torch transformers sentencepiece pandas numpy
    /tmp/enc2/bin/python src/embed2.py --model ibm/MoLFormer-XL-both-10pct

Writes data/emb_<slug>.npz with `train` (n_train, d) in rows.csv order and `test` (750, d) in
the blinded file's order -- the same contract as data/emb_chemprop_rdkit2d.npz.
"""
import argparse
import pathlib
import re
import time

import numpy as np
import pandas as pd
import torch

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = str(ROOT / "data") + "/"


def load(model_id, device):
    """Tokeniser and encoder, in eval mode.

    `trust_remote_code` is required by MoLFormer, whose linear-attention block lives in the
    repository rather than in transformers. It is off for everything else.
    """
    from transformers import AutoModel, AutoTokenizer

    remote = "molformer" in model_id.lower()
    tok = AutoTokenizer.from_pretrained(model_id, trust_remote_code=remote)
    mod = AutoModel.from_pretrained(model_id, trust_remote_code=remote,
                                    deterministic_eval=True) if remote else \
        AutoModel.from_pretrained(model_id)
    mod.eval().to(device)
    return tok, mod


def embed(smiles, tok, mod, device, pool="mean", batch=64, maxlen=202):
    """One vector per molecule, in the order given.

    Mean pooling over non-padding tokens rather than the CLS vector, and this is a correctness
    requirement rather than a preference. Loading ChemBERTa-77M-MTR prints

        pooler.dense.weight | MISSING | newly initialized

    because the published checkpoint carries a regression head and no pooler. Anything that
    reads `pooler_output` -- or CLS, which reaches the same untrained layer -- would be running
    the molecule through **randomly initialised weights** and returning a plausible-looking
    vector that means nothing. Mean pooling reads `last_hidden_state` directly and never
    touches it.

    The mask is applied explicitly because padding tokens otherwise drag every short molecule
    toward the same vector, which is the second way this silently degrades.
    """
    out = []
    with torch.no_grad():
        for i in range(0, len(smiles), batch):
            chunk = smiles[i:i + batch]
            enc = tok(chunk, padding=True, truncation=True, max_length=maxlen,
                      return_tensors="pt").to(device)
            h = mod(**enc).last_hidden_state
            if pool == "cls":
                v = h[:, 0]
            else:
                m = enc["attention_mask"].unsqueeze(-1).to(h.dtype)
                v = (h * m).sum(1) / m.sum(1).clamp(min=1)
            out.append(v.float().cpu().numpy())
    return np.vstack(out).astype(np.float32)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--out", default="")
    ap.add_argument("--pool", default="mean", choices=["mean", "cls"])
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--device", default="mps" if torch.backends.mps.is_available() else "cpu")
    ap.add_argument("--limit", type=int, default=0, help="только первые N, для проверки")
    a = ap.parse_args()

    slug = re.sub(r"[^a-z0-9]+", "_", a.model.lower()).strip("_")
    out = a.out or D + f"emb_{slug}.npz"

    tok, mod = load(a.model, a.device)
    npar = sum(p.numel() for p in mod.parameters())
    print(f"{a.model}: {npar:,} параметров, устройство {a.device}, пулинг {a.pool}", flush=True)

    rows = pd.read_csv(D + "rows.csv")
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    tr_smi, te_smi = list(rows.SMILES), list(te.SMILES)
    if a.limit:
        tr_smi, te_smi = tr_smi[:a.limit], te_smi[:a.limit]

    t0 = time.time()
    E_tr = embed(tr_smi, tok, mod, a.device, a.pool, a.batch)
    print(f"  обучение: {E_tr.shape} за {time.time()-t0:.0f} с", flush=True)
    t0 = time.time()
    E_te = embed(te_smi, tok, mod, a.device, a.pool, a.batch)
    print(f"  тест:     {E_te.shape} за {time.time()-t0:.0f} с", flush=True)

    # Вырожденный эмбеддинг выглядит как обычный и ничего не значит. Две проверки:
    # сколько колонок вообще шевелится и не совпали ли две разные молекулы.
    sd = E_tr.std(0)
    dead = int((sd < 1e-6).sum())
    u = np.unique(np.round(E_tr, 4), axis=0)
    print(f"  мёртвых колонок {dead} из {E_tr.shape[1]}, "
          f"различных векторов {len(u)} из {len(E_tr)}", flush=True)
    if dead > E_tr.shape[1] // 2 or len(u) < len(E_tr) // 2:
        raise SystemExit("эмбеддинг вырожден: половина колонок стоит или молекулы слиплись")

    np.savez_compressed(out, train=E_tr, test=E_te,
                        model=np.array(a.model), pool=np.array(a.pool))
    print(f"сохранено: {out}")


if __name__ == "__main__":
    main()
