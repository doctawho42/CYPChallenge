"""Multi-head fine-tune of a published encoder on the union: our labels as the primary task, the
qHTS panel as auxiliary supervision. Writes files, scores nothing.

THE ONE THING THIS CHANGES FROM ITEM 343. That item fine-tuned the SAME checkpoint
(`chemprop_medium`) from the same warm start with ONE head and measured 0.5613 of solo macro rank
against the frozen probe's 0.6261 -- fine-tuning destroyed 0.0648 of what freezing had. The reading
there was that 6525 labels overfit the encoder. This file adds the auxiliary heads and nothing else:
same checkpoint, same encoder-loading code (imported, not copied), same folds, same two-group
optimiser. The label ratio is 1 to 9.5, and if the reading was right that is the difference.

Item 341 sizes the prize off `briford`'s published shape at >= +0.048 of rho. Items 363 and 364
measured that neither ChEMBL nor the panel transfers ordering ALONE, which is the argument for the
union rather than for either source, and item 356 measured the panel as a MEMBER at -0.0092 over the
ensemble -- a different device from this one, and item 368 records why the distinction matters: a
fine-tune makes an encoder, not a member, and an encoder is the only thing that has moved the board.

WHAT IS TRAINED AND WHAT IS SCORED. Three heads over one shared encoder:

    head_main   4 выхода   наш pIC50, маскированная MSE      6525 меток
    head_bin    5 выходов  активен/неактивен из панели, BCE  41530 меток
    head_reg    5 выходов  pAC50 панели, маскированная MSE   20760 меток

The panel's pAC50 is left-truncated at 4.10 against our 1.91 and its means run 5.18--5.36 against
our 4.10--4.96 (item 368), so it stays on its OWN heads: the encoder takes the chemistry and the
primary head keeps our scale. If that separation leaks the arm will look like items 363/364.

Panel-only molecules carry fold -1 and are in every training split. The 182 molecules that are in
both corpora carry their CHALLENGE fold, asserted by `src/union.py`, because a panel row for a
held-out molecule would make the out-of-fold estimate worthless.

TARGETS ARE Z-SCORED HERE, EXPLICITLY, and that is deliberate. Item 343's run shipped predictions in
each fold's scaler space because `RegressionFFN` was built without an `output_transform`, and the
pair came out 4.70 against our usual 0.65. Writing the standardisation by hand with the statistics
saved beside the predictions removes the whole class: there is no hidden transform to forget.

TWO ARMS COME OUT, because the question is about the ENCODER and not about this model's own head:

    `oof`  -- the main head's own out-of-fold predictions, directly comparable to item 343's 0.5613
    `emb`  -- the fold's encoder applied to ALL 4905 challenge molecules, so the shipped RidgeCV
              probe can be run on it out of fold, which is the device that measured 0.6261

Scoring happens in the PINNED environment, never here: chemprop resolves scikit-learn to 1.9.1 and
CLAUDE.md records that 1.9.0 moves every prediction by ~0.135 pIC50. The contract between the two
halves is a file, as `src/embed.py` established.

    <env>/bin/python src/ftunion.py --ckpt <...>/chemprop_medium.pt --union data/union_seed0.npz \\
        --tag unionA --epochs 40 --accelerator gpu --devices 0

Writes <out>/ft_union_<tag>.npz. Exit: 0 ran, 2 a control failed, 3 an input is missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))

import argparse
import json
import time

import numpy as np
import torch
import torch.nn as nn
from torch import optim

N_MAIN, N_BIN, N_REG = 4, 5, 5
NT = N_MAIN + N_BIN + N_REG          # 14 target columns packed into one y vector


def head(din, dout, hid, layers):
    """A plain MLP head. Written out rather than reusing chemprop's FFN classes because those
    carry their own criterion and output transform, which is defect 1 of item 343's file."""
    mods, d = [], din
    for _ in range(max(0, layers)):
        mods += [nn.Linear(d, hid), nn.ReLU()]
        d = hid
    mods += [nn.Linear(d, dout)]
    return nn.Sequential(*mods)


def masked_mse(p, t):
    m = torch.isfinite(t)
    if not m.any():
        return p.sum() * 0.0
    return ((p[m] - t[m]) ** 2).mean()


def masked_bce(p, t):
    m = torch.isfinite(t)
    if not m.any():
        return p.sum() * 0.0
    return nn.functional.binary_cross_entropy_with_logits(p[m], t[m])


def make_net(mp, hid, layers, lam_bin, lam_reg, lr_enc, lr_head):
    import lightning.pytorch as pl
    from chemprop.nn import MeanAggregation

    class UnionNet(pl.LightningModule):
        def __init__(self):
            super().__init__()
            self.mp = mp
            self.agg = MeanAggregation()
            self.bn = nn.BatchNorm1d(mp.output_dim)
            self.h_main = head(mp.output_dim, N_MAIN, hid, layers)
            self.h_bin = head(mp.output_dim, N_BIN, hid, layers)
            self.h_reg = head(mp.output_dim, N_REG, hid, layers)

        def embed(self, bmg):
            return self.bn(self.agg(self.mp(bmg), bmg.batch))

        def forward(self, bmg):
            z = self.embed(bmg)
            return self.h_main(z), self.h_bin(z), self.h_reg(z)

        def _loss(self, batch, tag):
            bmg, Y = batch.bmg, batch.Y
            pm, pb, pr = self(bmg)
            lm = masked_mse(pm, Y[:, :N_MAIN])
            lb = masked_bce(pb, Y[:, N_MAIN:N_MAIN + N_BIN])
            lr = masked_mse(pr, Y[:, N_MAIN + N_BIN:])
            tot = lm + lam_bin * lb + lam_reg * lr
            self.log_dict({tag + "_main": lm, tag + "_bin": lb, tag + "_reg": lr,
                           tag + "_loss": tot}, batch_size=len(Y), prog_bar=False)
            return tot

        def training_step(self, batch, _i):
            return self._loss(batch, "train")

        def validation_step(self, batch, _i):
            return self._loss(batch, "val")

        def configure_optimizers(self):
            enc = list(self.mp.parameters())
            rest = (list(self.agg.parameters()) + list(self.bn.parameters())
                    + list(self.h_main.parameters()) + list(self.h_bin.parameters())
                    + list(self.h_reg.parameters()))
            return optim.Adam([{"params": enc, "lr": lr_enc},
                               {"params": rest, "lr": lr_head}])

    return UnionNet()


def pack(Y, AB, AR, mu, sd, mur, sdr):
    """One (n, 14) target block: our z-scored pIC50, the panel's 0/1, the panel's z-scored pAC50.

    Standardisation uses the TRAINING rows' statistics, passed in, never the fold's own -- fitting a
    scaler on the union of a fold's halves is item 202's defect 2.
    """
    T = np.full((len(Y), NT), np.nan, np.float32)
    T[:, :N_MAIN] = (Y - mu) / sd
    T[:, N_MAIN:N_MAIN + N_BIN] = AB
    T[:, N_MAIN + N_BIN:] = (AR - mur) / sdr
    return T


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", required=True)
    ap.add_argument("--union", required=True)
    ap.add_argument("--tag", required=True)
    ap.add_argument("--out", default="/tmp/d42/out")
    ap.add_argument("--random", action="store_true",
                    help="случайный энкодер: различающий контроль пункта 343")
    ap.add_argument("--no-aux", dest="aux", action="store_false",
                    help="выключить вспомогательные головы: это и есть арм пункта 343")
    ap.add_argument("--folds", default="0,1,2,3,4")
    ap.add_argument("--epochs", type=int, default=40)
    ap.add_argument("--patience", type=int, default=8)
    ap.add_argument("--batch", type=int, default=64)
    ap.add_argument("--hid", type=int, default=512)
    ap.add_argument("--layers", type=int, default=2)
    ap.add_argument("--lr-enc", dest="lr_enc", type=float, default=1e-5)
    ap.add_argument("--lr-head", dest="lr_head", type=float, default=1e-3)
    ap.add_argument("--lam-bin", dest="lam_bin", type=float, default=1.0)
    ap.add_argument("--lam-reg", dest="lam_reg", type=float, default=1.0)
    ap.add_argument("--inner", type=float, default=0.1, help="доля обучения в валидацию")
    ap.add_argument("--accelerator", default="gpu")
    ap.add_argument("--devices", default="0")
    ap.add_argument("--monitor", default="val_loss",
                    help="что слушает ранняя остановка. val_loss это СУММА, где при lam=1 "
                         "вспомогательные члены дают ~90%% массы (62290 меток против 6525), "
                         "так что остановка может срабатывать по чужой задаче. val_main "
                         "слушает только нашу.")
    ap.add_argument("--workers", type=int, default=0,
                    help="воркеры загрузчика. 0 держит featurisation в главном процессе и даёт "
                         "~19%% занятости карты на этой машине; измеряйте, а не угадывайте")
    ap.add_argument("--limit", type=int, default=0, help="только N молекул, для прогона вхолостую")
    a = ap.parse_args(argv)
    if not a.aux:
        a.lam_bin = a.lam_reg = 0.0

    import lightning.pytorch as pl
    from chemprop.data import MoleculeDatapoint, MoleculeDataset, build_dataloader
    import finetune as FT

    for p in (a.ckpt, a.union):
        if not _pl.Path(p).exists():
            print("нет %s" % p)
            return 3
    _pl.Path(a.out).mkdir(parents=True, exist_ok=True)

    z = np.load(a.union, allow_pickle=True)
    smiles = [str(s) for s in z["smiles"]]
    fold, source = z["fold"], z["source"]
    Y, AB, AR = z["Y"].astype(np.float64), z["AB"].astype(np.float64), z["AR"].astype(np.float64)
    dg = str(z["digest"])
    n = len(smiles)
    ch = np.where(source == 0)[0]
    print("  объединение: %d молекул (%d челлендж), дайджест %s" % (n, len(ch), dg), flush=True)
    print("  метки: основных %d, бинарных %d, pAC50 %d; вспомогательные %s"
          % (int(np.isfinite(Y).sum()), int(np.isfinite(AB).sum()), int(np.isfinite(AR).sum()),
             "ВКЛЮЧЕНЫ" if a.aux else "выключены"), flush=True)
    if a.limit:
        keep = np.concatenate([ch[:a.limit], np.where(source == 1)[0][:a.limit]])
        print("  ПРОГОН ВХОЛОСТУЮ: оставляю %d молекул" % len(keep), flush=True)
        smiles = [smiles[i] for i in keep]
        fold, source, Y, AB, AR = fold[keep], source[keep], Y[keep], AB[keep], AR[keep]
        n = len(smiles)
        ch = np.where(source == 0)[0]

    folds = [int(x) for x in a.folds.split(",") if x.strip()]
    oof = np.full((n, N_MAIN), np.nan, np.float64)
    emb = {}
    log = []
    rng = np.random.default_rng(0)

    for f in folds:
        t0 = time.time()
        te = np.where(fold == f)[0]
        trall = np.where((fold != f) & (fold >= 0) | (fold == -1))[0]
        # CONTROL: the held-out fold must not appear in training, and panel-only rows must.
        if set(te) & set(trall):
            print("ОТКАЗ: фолд %d пересекается с обучением — утечка" % f)
            return 2
        if not (fold[trall] == -1).any():
            print("ОТКАЗ: в обучении нет панельных строк — объединение не подано")
            return 2
        va_n = max(1, int(a.inner * len(trall)))
        perm = rng.permutation(len(trall))
        va, tr = trall[perm[:va_n]], trall[perm[va_n:]]

        mu = np.nanmean(Y[tr], 0); sd = np.nanstd(Y[tr], 0)
        mur = np.nanmean(AR[tr], 0); sdr = np.nanstd(AR[tr], 0)
        sd = np.where(np.isfinite(sd) & (sd > 1e-6), sd, 1.0)
        sdr = np.where(np.isfinite(sdr) & (sdr > 1e-6), sdr, 1.0)
        mu = np.where(np.isfinite(mu), mu, 0.0); mur = np.where(np.isfinite(mur), mur, 0.0)
        T = pack(Y, AB, AR, mu, sd, mur, sdr)

        mk = lambda idx: MoleculeDataset([MoleculeDatapoint.from_smi(smiles[i], y=T[i])
                                          for i in idx])
        dtr, dva, dte = mk(tr), mk(va), mk(ch)
        ltr = build_dataloader(dtr, batch_size=a.batch, shuffle=True, drop_last=False,
                               num_workers=a.workers)
        lva = build_dataloader(dva, batch_size=a.batch, shuffle=False, drop_last=False,
                               num_workers=a.workers)
        lte = build_dataloader(dte, batch_size=a.batch, shuffle=False, drop_last=False,
                               num_workers=a.workers)

        model, cfg, loaded = FT.build(a.ckpt, a.random, a.hid, a.layers, a.lr_enc, a.lr_head)
        net = make_net(model.message_passing, a.hid, a.layers, a.lam_bin, a.lam_reg,
                       a.lr_enc, a.lr_head)
        cb = pl.callbacks.EarlyStopping(monitor=a.monitor, patience=a.patience, mode="min")
        dev = [int(x) for x in str(a.devices).split(",")] if a.accelerator == "gpu" else 1
        trainer = pl.Trainer(max_epochs=a.epochs, accelerator=a.accelerator, devices=dev,
                             logger=False, enable_checkpointing=False, callbacks=[cb],
                             enable_progress_bar=False, enable_model_summary=False)
        trainer.fit(net, ltr, lva)
        ep = int(trainer.current_epoch)

        net.eval()
        P, E = [], []
        with torch.no_grad():
            for b in lte:
                # Defect 8, in this file's own numbering: `BatchMolGraph.to(device)` moves the
                # tensors IN PLACE and returns None, unlike every torch module's `.to`. Assigning
                # its result gives None and the next line dies on `.V`. Call it for the effect.
                b.bmg.to(net.device)
                zz = net.embed(b.bmg)
                E.append(zz.cpu().numpy())
                P.append(net.h_main(zz).cpu().numpy())
        P = np.vstack(P) * sd + mu                      # back to pIC50, explicitly
        E = np.vstack(E).astype(np.float32)
        pos = {int(g): k for k, g in enumerate(ch)}
        for i in te:
            oof[i] = P[pos[int(i)]]
        emb["fold%d" % f] = E
        log.append({"fold": f, "n_train": int(len(tr)), "n_val": int(len(va)),
                    "n_panel_in_train": int((fold[tr] == -1).sum()),
                    "epochs": ep, "seconds": round(time.time() - t0, 1),
                    "loaded_params": int(loaded)})
        print("  фолд %d: обучение %d (панели %d), валидация %d, эпох %d, %.0f с"
              % (f, len(tr), int((fold[tr] == -1).sum()), len(va), ep, time.time() - t0),
              flush=True)

    out = "%s/ft_union_%s.npz" % (a.out, a.tag)
    np.savez_compressed(out, oof=oof[ch], emb_folds=np.array(sorted(emb), dtype=object),
                        digest=np.array(dg), tag=np.array(a.tag),
                        meta=np.array(json.dumps({"log": log, "aux": bool(a.aux),
                                                  "lam_bin": a.lam_bin, "lam_reg": a.lam_reg,
                                                  "epochs": a.epochs, "ckpt": a.ckpt})),
                        **{k: v for k, v in emb.items()})
    done = int(np.isfinite(oof[ch]).any(1).sum())
    print("\n  заполнено вне фолда: %d из %d челлендж-строк" % (done, len(ch)), flush=True)
    print("  записано: %s" % out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
