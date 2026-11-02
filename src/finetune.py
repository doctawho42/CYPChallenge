"""Fine-tune a published chemprop encoder on our labels and our folds. Writes a file, scores nothing.

Why this exists. Item 341 found that a representation learned from the molecular graph is the first
thing in this project whose errors fall outside the member band, and item 342 took it from harmful
to neutral with the screen as a second head. Both were trained from scratch. The remaining gap to
the shipped probe --- a PRETRAINED encoder with a trivial ridge head, 0.6261 of macro rank against
our 0.6079 --- is what this measures: what pretraining is worth when the encoder is allowed to move.

Runs in a SEPARATE environment, and that is not a preference. chemprop pulls forty-one dependencies
and resolves `scikit-learn` to 1.9.1, while `CLAUDE.md` records that 1.9.0 moves every prediction in
this pipeline by ~0.135 pIC50. So this script must not score anything: it writes out-of-fold
predictions and `verify/k120_finetune.py` scores them in the pinned environment, exactly as
`src/embed.py` writes an array and everything downstream reads it. The contract is a file.

    uv venv <scratch>/enc --python 3.12
    uv pip install --python <scratch>/enc/bin/python chemprop "torch==2.13.0" "rdkit==2025.9.6"
    curl -L -o <scratch>/rdkit2d__s42.pt https://huggingface.co/openadmet/\\
chemprop-foundation-pretraining-weights/resolve/main/rdkit2d__s42.pt
    <scratch>/enc/bin/python src/finetune.py --ckpt <scratch>/rdkit2d__s42.pt --tag rdkit2d

Pinning torch and rdkit to the versions the main environment already has is deliberate: both come
out of the uv cache as hardlinks, so the encoder environment costs 26 MiB rather than several GB.

ARMS. `--ckpt` fine-tunes from a published checkpoint; `--random` is the DISCRIMINATING CONTROL ---
the same architecture, the same schedule, weights initialised fresh, so a win by a checkpoint is
attributable to the pretrained weights and not to chemprop's architecture. `--aux` adds four heads
on the single-point screen's log2fc, the device item 342 measured at +0.0168 of solo rank on our own
net and the configuration the checkpoint's own author reports as his best-ranking chemprop model.

SEVEN DEFECTS OF THIS API, each measured rather than guessed, each pre-empted below. They are listed
because every one of them fails SILENTLY and three would have corrupted the whole experiment.

  1. `RegressionFFN.train_step` is an alias of `forward`, so it runs `output_transform`. In
     `model.train()` the transform is a no-op and the loss lives in scaled space; the moment you
     call `model.eval()` for validation it switches ON, predictions return to pIC50 while the
     targets are still z-scored, and the validation loss comes out 21.50 where it should be 1.08.
     Lightning handles it in `MPNN.on_validation_model_eval`; a hand-written loop never calls that
     hook. This file therefore uses Lightning's `Trainer` with a real val dataloader rather than
     its own loop, so the hook fires.
  2. `model.criterion` is a `torchmetrics.Metric` and accumulates forever. Gradients are safe
     (`full_state_update = False` makes `forward` return the batch value) but `compute()` returns
     the mean since construction. Lightning resets it through `self.log(..., on_epoch=True)`.
  3. `build_dataloader` sets `drop_last=True` by itself when `len(dataset) % batch_size == 1`, with
     only a `logger.warning`. On our sizes fold 0's training half is 3889 rows, which triggers at
     batch 8 AND at batch 16; fold 2's is 3913, which triggers at 8. One molecule vanishes from one
     fold and not another, so the folds stop being comparable element for element --- which is what
     the whole repository rests on. `drop_last=False` is passed explicitly everywhere.
  4. `build_dataloader` defaults to `shuffle=True`, unlike torch's own DataLoader. A prediction
     loader left at the default gets a `RandomSampler` and the row correspondence is lost in
     silence. `shuffle=False` is passed explicitly for every prediction pass.
  5. `RegressionFFN`'s `input_dim` defaults to 300, not to the encoder's `d_h`, and the dimension
     check `MPNN`'s docstring promises does not exist in the code --- the model assembles quietly
     and dies on the first batch. `input_dim=mp.output_dim` is read off the block.
  6. A separate learning rate for the encoder does not exist: `MPNN.configure_optimizers` is one
     `Adam` over `self.parameters()` at `init_lr`. It is overridden below, and the inherited
     Noam-like scheduler is dropped with it, because that schedule normalises on a single `init_lr`
     and multiplies both groups' base rates by it.
  7. `trainer.predict` returns a LIST of per-batch tensors, not a tensor.

Writes data/ft_<tag>_seed<seed>.npz with `oof` (n_rows, 4) in rows.csv order, NaN where a fold was
not run, plus the per-fold epoch counts and losses as arrays. Exit: 0 ran, 2 a control failed.
"""

import argparse
import json
import pathlib
import sys
import time

import numpy as np
import pandas as pd
import torch
from torch import optim

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = str(ROOT / "data") + "/"
sys.path.insert(0, str(ROOT))

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
GOLDEN = "2d93c19815e14261"


def build(ckpt, random_init, hid, layers, lr_enc, lr_head):
    """The encoder from a published blob (or freshly initialised) plus our own regression head."""
    from chemprop.models import MPNN
    from chemprop.nn import BondMessagePassing, MeanAggregation, RegressionFFN

    blob = torch.load(ckpt, map_location="cpu", weights_only=False)
    hp = blob["hyper_parameters"]
    keys = ("d_v", "d_e", "d_h", "bias", "depth", "dropout", "activation", "undirected")
    # A medium-style blob stores the whole model and nests the encoder's settings; a pretraining
    # blob stores the encoder alone at the top level. Both are read, neither is guessed.
    src = hp if all(k in hp for k in keys) else hp["message_passing"]
    cfg = {k: (src[k] if not hasattr(src, "get") else src.get(k, hp.get(k))) for k in keys}
    mp = BondMessagePassing(**cfg)

    loaded = 0
    if not random_init:
        sd = blob["state_dict"]
        # Two shapes of blob, read explicitly rather than pattern-matched. A pretraining
        # checkpoint stores the encoder alone (`W_i.weight`); a whole-model checkpoint prefixes
        # everything (`message_passing.W_i.weight`) and also carries the predictor, which must
        # NOT be loaded here because our head has a different task count.
        pref = "message_passing."
        if any(k.startswith(pref) for k in sd):
            own = {k[len(pref):]: v for k, v in sd.items() if k.startswith(pref)}
        else:
            own = dict(sd)
        want = set(mp.state_dict())
        extra = sorted(set(own) - want)
        own = {k: v for k, v in own.items() if k in want}
        if extra:
            print("      из чекпойнта не берутся (не часть энкодера): %s" % extra[:6], flush=True)
        rep = mp.load_state_dict(own, strict=False)
        if rep.unexpected_keys:
            print("ОТКАЗ: лишние ключи при загрузке энкодера:", rep.unexpected_keys)
            raise SystemExit(2)
        if rep.missing_keys:
            print("ОТКАЗ: энкодер загрузился НЕ полностью, пропущено:", rep.missing_keys)
            raise SystemExit(2)
        loaded = sum(v.numel() for v in own.values())
        print("      энкодер загружен: %d тензоров, %.2f млн параметров"
              % (len(own), loaded / 1e6), flush=True)
    else:
        print("      энкодер СЛУЧАЙНЫЙ (различающий контроль)", flush=True)

    # Defect 5: input_dim defaults to 300 and nothing checks it.
    ffn = RegressionFFN(n_tasks=4, input_dim=mp.output_dim, hidden_dim=hid,
                        n_layers=layers, dropout=0.0)

    class FineTune(MPNN):
        """Defect 6: one Adam group at init_lr upstream, and the Noam schedule goes with it."""

        def configure_optimizers(self):
            enc = list(self.message_passing.parameters())
            rest = (list(self.agg.parameters()) + list(self.bn.parameters())
                    + list(self.predictor.parameters()))
            return optim.Adam([{"params": enc, "lr": lr_enc},
                               {"params": rest, "lr": lr_head}])

    model = FineTune(mp, MeanAggregation(), ffn, batch_norm=True)
    tot = sum(p.numel() for p in model.parameters())
    print("      модель: %.2f млн параметров (энкодер %.2f млн, d_h %d, depth %d)"
          % (tot / 1e6, sum(p.numel() for p in mp.parameters()) / 1e6,
             cfg["d_h"], cfg["depth"]), flush=True)
    return model, cfg, loaded


def datasets(smiles, Y, idx_tr, idx_va, idx_te):
    from chemprop.data import MoleculeDatapoint, MoleculeDataset

    def mk(idx):
        return MoleculeDataset([MoleculeDatapoint.from_smi(smiles[i], y=Y[i]) for i in idx])

    dtr, dva, dte = mk(idx_tr), mk(idx_va), mk(idx_te)
    sc = dtr.normalize_targets()          # NaN-aware StandardScaler, fitted on TRAIN only
    dva.normalize_targets(sc)
    dte.normalize_targets(sc)
    return dtr, dva, dte, sc


def run_fold(args, smiles, Y, fold, f, log):
    from chemprop.data import build_dataloader
    from lightning import pytorch as pl
    from lightning.pytorch.callbacks import EarlyStopping

    rng = np.random.default_rng(args.seed * 100 + f)
    tr_all = np.flatnonzero(fold != f)
    te = np.flatnonzero(fold == f)
    rng.shuffle(tr_all)
    ncut = max(1, int(len(tr_all) * args.inner))
    va, tr = tr_all[:ncut], tr_all[ncut:]

    dtr, dva, dte, _sc = datasets(smiles, Y, tr, va, te)
    model, _cfg, loaded = build(args.ckpt, args.random, args.hid, args.layers,
                               args.lr_enc, args.lr_head)

    # Defects 3 and 4: both passed explicitly on every loader.
    ltr = build_dataloader(dtr, batch_size=args.batch, shuffle=True, drop_last=False,
                           seed=args.seed * 100 + f, num_workers=0)
    lva = build_dataloader(dva, batch_size=args.batch, shuffle=False, drop_last=False,
                           num_workers=0)
    lte = build_dataloader(dte, batch_size=args.batch, shuffle=False, drop_last=False,
                           num_workers=0)

    t0 = time.time()
    stop = EarlyStopping(monitor="val_loss", patience=args.patience, mode="min")
    trainer = pl.Trainer(max_epochs=args.epochs, accelerator="cpu", devices=1,
                         logger=False, enable_checkpointing=False, callbacks=[stop],
                         enable_progress_bar=False, enable_model_summary=False,
                         log_every_n_steps=10 ** 9)
    trainer.fit(model, ltr, lva)
    # Defect 7: a list of per-batch tensors.
    pred = torch.cat(trainer.predict(model, lte)).numpy()
    if pred.ndim == 3:
        pred = pred[..., 0]
    dt = time.time() - t0

    best = float(stop.best_score) if stop.best_score is not None else float("nan")
    ep = int(trainer.current_epoch)
    log.append({"fold": int(f), "epochs": ep, "val_best_scaled": best,
                "seconds": round(dt, 1), "n_train": int(len(tr)), "n_inner": int(len(va)),
                "n_held": int(len(te)), "loaded_params": int(loaded)})
    print("      фолд %d: %d эпох, %.0f с, лучшая валидационная (масштаб.) %.4f, "
          "обучающих %d, отложенных %d"
          % (f, ep, dt, best, len(tr), len(te)), flush=True)
    if pred.shape != (len(te), 4):
        print("ОТКАЗ: форма предсказаний %s против ожидаемой %s"
              % (pred.shape, (len(te), 4)))
        raise SystemExit(2)
    return te, pred


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", default=None, help="путь к опубликованному чекпойнту")
    ap.add_argument("--random", action="store_true",
                    help="РАЗЛИЧАЮЩИЙ КОНТРОЛЬ: та же архитектура, случайные веса")
    ap.add_argument("--tag", required=True, help="имя арма в имени выходного файла")
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--folds", default="0,1,2,3,4")
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--patience", type=int, default=6)
    ap.add_argument("--batch", type=int, default=32)
    ap.add_argument("--hid", type=int, default=512)
    ap.add_argument("--layers", type=int, default=2)
    ap.add_argument("--lr-enc", dest="lr_enc", type=float, default=1e-5)
    ap.add_argument("--lr-head", dest="lr_head", type=float, default=1e-3)
    ap.add_argument("--inner", type=float, default=0.1)
    a = ap.parse_args(argv)
    if a.ckpt is None:
        print("нужен --ckpt: даже --random читает из него гиперпараметры архитектуры")
        return 2
    folds = [int(x) for x in a.folds.split(",") if x.strip()]

    from cypsplit import butina_folds, fold_digest

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    smiles = list(rows.SMILES)
    fold, ncl = butina_folds(smiles, seed=a.seed)
    dg = fold_digest(fold)
    print("  сид %d: дайджест %s, кластеров %d" % (a.seed, dg, ncl), flush=True)
    if a.seed == 0 and dg != GOLDEN:
        print("ОТКАЗ: разбиение сдвинулось, сид 0 дал %s" % dg)
        return 2

    # count first: the label matrix this arm is trained against
    M = ~np.isnan(Y)
    print("  меток по ферментам: " + "  ".join("%s %d" % (c[3:], M[:, e].sum())
                                               for e, c in enumerate(CYPS))
          + "  (всего %d из %d клеток)" % (M.sum(), M.size), flush=True)
    perm = M.sum(1)
    print("  молекул без ни одной метки: %d  (с одной %d, двумя %d, тремя %d, четырьмя %d)"
          % ((perm == 0).sum(), (perm == 1).sum(), (perm == 2).sum(),
             (perm == 3).sum(), (perm == 4).sum()), flush=True)
    if (perm == 0).sum():
        print("ОТКАЗ: есть молекулы без метк --- маска партии может стать нулевой")
        return 2

    oof = np.full((len(smiles), 4), np.nan)
    log = []
    for f in folds:
        te, pred = run_fold(a, smiles, Y, fold, f, log)
        oof[te] = pred

    out = D + "ft_%s_seed%d.npz" % (a.tag, a.seed)
    np.savez(out, oof=oof, digest=np.array(dg), tag=np.array(a.tag),
             meta=np.array(json.dumps({"log": log, "argv": vars(a)}, ensure_ascii=False)))
    done = int(np.isfinite(oof).all(1).sum())
    print("  сохранено: %s  (заполнено %d строк из %d)" % (out, done, len(smiles)), flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
