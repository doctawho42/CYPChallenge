"""Embed and PREDICT with the organisers' published CYP model --- the one checkpoint in their
collection that was trained on the endpoint we are scored on, and the one this project refused.

WHY IT WAS REFUSED AND WHY THAT REASON DOES NOT APPLY HERE. `src/embed.py` records the decision
verbatim: "The CYP-fine-tuned baseline in the same organisation was NOT used, because a model
already trained on CYP pIC50 data would make the comparison circular." That is correct about the
comparison `embed.py` was built to make -- does a LEARNED representation of the descriptors beat
handing the descriptors over directly -- and it is irrelevant to whether the checkpoint is a member.
A model trained on public ChEMBL CYP data is external prior information of exactly the kind item
350 named as the remaining form of the transfer axis, and it costs no GPU because the training has
already been done by someone else.

WHAT THE MODEL IS. `openadmet/cyp1a2-cyp2d6-cyp3a4-cyp2c9-chemeleon-v1`: a multitask chemprop
D-MPNN warm-started from CheMeleon, 2048-wide message passing into a 1024x3 head with four
outputs, trained to convergence on 8068 ChEMBL 37 structures with no validation split. Its own
model card calls the data public and curated.

WHAT IS MEASURED BEFORE IT IS USED, and both numbers are gates rather than notes:

  * **Overlap with the blinded test set is ZERO** of 750, by canonical SMILES of the largest
    fragment, with the matcher checked against the test set itself first. The overlap with our
    TRAINING set is 64, which reproduces item 290's recorded count exactly.
  * **The checkpoint's unspecified hyper-parameters are identified, not guessed.** `model.json`
    pins the head and the dropout but not the message-passing depth, the activation or the
    aggregation, and a bare state dict loads for any of them because W_h is shared across depth
    steps -- so a wrong choice gives plausible numbers silently. The model is a no-split model
    trained to convergence, so the right configuration must reproduce its OWN training labels and
    the wrong ones must not. Over depth 2-6 x {relu, leakyrelu} x {mean, norm, sum}:

        depth 6, relu, mean   ро 0.9987   R2  0.9977   MAE 0.0292
        depth 6, relu, sum    ро 0.9329   R2 -25.52    MAE 2.2331   <- the runner-up
        depth 5, relu, mean   ро 0.7319   R2  0.4804   MAE 0.4457

    A gap of 0.0658 of rank and three orders of magnitude of MAE to the next candidate. The check
    is re-run inside this script on every invocation and it STOPS if the fit falls below 0.99, so
    the configuration cannot drift silently.

WHAT IS WRITTEN. Three arrays per split, because the project's probe device wants an embedding and
the endpoint model also has something an encoder does not --- four direct estimates of our targets:

    pred  (n, 4)     the model's pIC50 for CYP3A4, CYP2D6, CYP2C9, CYP1A2, in ITS column order
    emb   (n, 2048)  agg(message_passing(.)), the same graph embedding the shipped probe reads
    ffn   (n, 1024)  the penultimate head layer, which the 2048 does not contain

The ChEMBL offset that killed the merge (items 60-63, +0.22 to +0.60 against our labels) is a
statement about LEVEL. Spearman is invariant to it exactly, and the affine pair absorbs what is
left for ST-RAE, so the quantity this project is graded on is the one least exposed to the defect
that closed the external-rows axis.

    /tmp/encenv/bin/python src/embed_cypmodel.py \\
        --model <dir>/anvil_training --rows data/rows.csv --test data/cyp-challenge-TEST-BLINDED.csv \\
        --out data/emb_oadmetcyp

Writes <out>_train.npz and <out>_test.npz, each with `pred`, `emb`, `ffn` and `smiles`, in the row
order of the file given. Runs out of process for the reason `src/embed.py` gives: chemprop must not
become a runtime dependency of the pipeline, and the contract between the two environments is a file.
Exit: 0 ran, 2 a control failed, 3 an input is missing.
"""
import argparse
import pathlib
import sys

import numpy as np
import pandas as pd
import torch
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

DEPTH, ACT, AGG = 6, "relu", "mean"      # identified; see the docstring
FIT_FLOOR = 0.99                          # the identification must still hold on every run


def canon(s):
    m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
    if m is None:
        return None
    fr = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
    if len(fr) > 1:
        m = max(fr, key=lambda f: f.GetNumHeavyAtoms())
    try:
        return Chem.MolToSmiles(m)
    except Exception:
        return None


def build(model_dir, random_seed=None):
    from chemprop.models.model import MPNN
    from chemprop.nn import (BondMessagePassing, MeanAggregation, NormAggregation,
                             RegressionFFN, SumAggregation)
    from chemprop.nn.transforms import UnscaleTransform

    sd = torch.load(str(pathlib.Path(model_dir) / "model.pth"), map_location="cpu",
                    weights_only=True)
    agg = {"mean": MeanAggregation, "norm": NormAggregation, "sum": SumAggregation}[AGG]
    mp = BondMessagePassing(d_v=72, d_e=14, d_h=2048, depth=DEPTH, bias=False,
                            dropout=0.0, activation=ACT)
    ffn = RegressionFFN(n_tasks=4, input_dim=2048, hidden_dim=1024, n_layers=3, dropout=0.1,
                        activation=ACT,
                        output_transform=UnscaleTransform(
                            sd["predictor.output_transform.mean"].numpy().ravel(),
                            sd["predictor.output_transform.scale"].numpy().ravel()))
    m = MPNN(mp, agg(), ffn, batch_norm=False)
    if random_seed is not None:
        # Item 343's control, on THIS architecture: the same shape with untrained weights. Only
        # the output transform is kept, so the arm differs from the real one in the weights alone.
        torch.manual_seed(random_seed)
        for mod in m.modules():
            if hasattr(mod, "reset_parameters"):
                mod.reset_parameters()
        m.eval()
        return m
    clean = {k: v for k, v in sd.items() if not k.startswith("metrics")}
    missing, unexpected = m.load_state_dict(clean, strict=False)
    missing = [k for k in missing if not k.startswith("metrics")]
    if missing:
        print("КОНТРОЛЬ: веса не загрузились: %r" % missing[:6], file=sys.stderr)
        raise SystemExit(2)
    m.eval()
    return m


def run(model, smiles, batch=128):
    """pred, emb and ffn in the order given. drop_last is passed explicitly for the reason
    src/embed_medium.py records: chemprop drops a size-1 remainder silently."""
    from chemprop import data, featurizers

    feat = featurizers.SimpleMoleculeMolGraphFeaturizer()
    pts = [data.MoleculeDatapoint.from_smi(s, [0.0] * 4) for s in smiles]
    dl = data.build_dataloader(data.MoleculeDataset(pts, feat), num_workers=0,
                               batch_size=batch, shuffle=False, drop_last=False)
    P, E, F = [], [], []
    with torch.inference_mode():
        for b in dl:
            bmg, V_d, X_d, *_ = b
            h = model.agg(model.message_passing(bmg, V_d), bmg.batch)
            E.append(h.numpy())
            z = h
            for blk in list(model.predictor.ffn)[:-1]:
                z = blk(z)
            F.append(z.numpy())
            P.append(model(bmg, V_d, X_d).numpy())
    P, E, F = (np.concatenate(x, 0) for x in (P, E, F))
    if not (len(P) == len(E) == len(F) == len(smiles)):
        print("КОНТРОЛЬ: %d/%d/%d против %d молекул --- строки разъехались"
              % (len(P), len(E), len(F), len(smiles)), file=sys.stderr)
        raise SystemExit(2)
    return P.astype(np.float32), E.astype(np.float32), F.astype(np.float32)


def main():
    from scipy.stats import spearmanr

    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True, help="the anvil_training directory")
    ap.add_argument("--rows", required=True)
    ap.add_argument("--test", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--batch", type=int, default=128)
    ap.add_argument("--random", type=int, default=None,
                    help="item 343's control: same architecture, untrained weights at this seed")
    a = ap.parse_args()

    md = pathlib.Path(a.model)
    for f in ("model.pth", "data/X_train.csv", "data/y_train.csv"):
        if not (md / f).exists():
            print("нет %s" % (md / f), file=sys.stderr)
            return 3
    model = build(md, a.random)
    print("модель: depth=%d act=%s agg=%s, выход %d задач%s"
          % (DEPTH, ACT, AGG, model.predictor.ffn[-1][-1].out_features,
             "" if a.random is None else "  [СЛУЧАЙНЫЕ ВЕСА, сид %d --- контроль пункта 343]"
             % a.random), flush=True)

    # ---- control 1: the identification must still hold on this invocation ----------------
    X = pd.read_csv(md / "data/X_train.csv")
    Y = pd.read_csv(md / "data/y_train.csv")
    k = np.where(Y.notna().any(axis=1).to_numpy())[0][:384]
    Po, _, _ = run(model, [X.iloc[i, 0] for i in k], a.batch)
    yt = Y.iloc[k].to_numpy(float)
    rh = []
    for j in range(4):
        m_ = ~np.isnan(yt[:, j])
        if m_.sum() >= 30:
            rh.append(float(spearmanr(yt[m_, j], Po[m_, j]).statistic))
    fit = float(np.mean(rh))
    print("контроль 1: модель воспроизводит свои обучающие метки, ро %.4f (%s), порог %.2f"
          % (fit, " ".join("%.3f" % x for x in rh), FIT_FLOOR), flush=True)
    if a.random is None and fit < FIT_FLOOR:
        print("КОНТРОЛЬ 1 ПРОВАЛЕН: ро %.4f --- конфигурация не та, что опознана" % fit,
              file=sys.stderr)
        return 2
    if a.random is not None and fit >= FIT_FLOOR:
        print("КОНТРОЛЬ 1 ПРОВАЛЕН: случайные веса воспроизводят метки при ро %.4f --- "
              "значит они не случайны" % fit, file=sys.stderr)
        return 2

    # ---- control 2: the overlap with the blinded test set, measured before anything else ---
    rows = pd.read_csv(a.rows)
    te = pd.read_csv(a.test)
    ext = set(filter(None, (canon(s) for s in X.iloc[:, 0])))
    tek = [canon(s) for s in te.SMILES]
    if None in tek or not tek:
        print("КОНТРОЛЬ 2: тест не разбирается RDKit", file=sys.stderr)
        return 2
    if sum(1 for s in tek if s in set(tek)) != len(tek):
        print("КОНТРОЛЬ 2 ПРОВАЛЕН: сопоставитель не находит тест в самом себе", file=sys.stderr)
        return 2
    ours = set(filter(None, (canon(s) for s in rows.SMILES)))
    hit_te, hit_tr = len(ext & set(tek)), len(ext & ours)
    print("контроль 2: обучение модели %d структур; пересечение со СЛЕПЫМ тестом %d из 750, "
          "с нашим обучением %d (пункт 290: 64)" % (len(ext), hit_te, hit_tr), flush=True)
    if hit_te:
        print("КОНТРОЛЬ 2 ПРОВАЛЕН: модель видела %d тестовых структур --- это надо решать "
              "правилами соревнования, а не здесь" % hit_te, file=sys.stderr)
        return 2

    for name, df in (("train", rows), ("test", te)):
        smi = list(df.SMILES)
        bad = [i for i, s in enumerate(smi) if Chem.MolFromSmiles(s) is None]
        if bad:
            print("КОНТРОЛЬ: %d неразбираемых строк в %s" % (len(bad), name), file=sys.stderr)
            return 2
        P, E, F = run(model, smi, a.batch)
        # ---- control 3: determinism, and order preserved under a permutation -------------
        P2, E2, _ = run(model, smi[:256], a.batch)
        if not (np.allclose(P2, P[:256], atol=1e-5) and np.allclose(E2, E[:256], atol=1e-5)):
            print("КОНТРОЛЬ 3 ПРОВАЛЕН: повторный проход даёт другое", file=sys.stderr)
            return 2
        pm = np.random.default_rng(0).permutation(len(smi))[:256]
        Pp, _, _ = run(model, [smi[i] for i in pm], a.batch)
        if not np.allclose(Pp, P[pm], atol=1e-5):
            print("КОНТРОЛЬ 3 ПРОВАЛЕН: перестановка не сохраняет соответствие строк",
                  file=sys.stderr)
            return 2
        out = "%s_%s.npz" % (a.out, name)
        np.savez_compressed(out, pred=P, emb=E, ffn=F, smiles=np.array(smi, dtype=object),
                            tasks=np.array(list(Y.columns), dtype=object))
        print("%s: pred %s, emb %s, ffn %s -> %s  (контроль 3: детерминизм и порядок держатся)"
              % (name, P.shape, E.shape, F.shape, out), flush=True)
        print("    предсказания по колонкам модели: " +
              ", ".join("%s %.2f..%.2f (ско %.2f)"
                        % (c.replace("OPENADMET_LOGAC50_", ""), P[:, j].min(), P[:, j].max(),
                           P[:, j].std()) for j, c in enumerate(Y.columns)), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
