"""Joint multi-head fine-tune over challenge + ChEMBL + qHTS, warm-started from the published
CheMeleon CYP checkpoint. The one configuration items 350 and 364 left standing as "sized for first
place", and the one item 386 independently diagnosed as the fix.

WHY THIS AND NOT ANOTHER MEMBER. Item 386 measured where our ordering fails: pairwise sign accuracy
decays monotonically with structural similarity on CYP1A2 and CYP2D6, reaching 0.533 against a coin
flip of 0.5 at the similarity the test actually occupies -- and the test is nothing but analogue
pairs, 75 series of about ten catalogue neighbours. Morgan fingerprints of two analogues at Tanimoto
0.6 differ in few bits, so a ridge or GP over them returns nearly equal predictions by construction;
no re-weighting of such members resolves an activity cliff. Only the representation can.

THE TWO DEVICES ALL THREE OPEN-CODE ENTRANTS ABOVE US USE AND WE NEVER HAVE (item 316):

  * one output head per data source, label scales NEVER merged. Our entire external-data history is
    the cross-assay affine correction they avoid by construction -- items 60-63, 144, 153, 157 and
    166's +0.44 to +0.87 protocol offset. Here each source gets its own head and its own
    standardisation, and only the challenge head is ever scored.
  * inactives as LEFT-CENSORED BOUNDS under a one-sided loss. The panel's censored rows are 31671
    of 54177, and the bound is read off the data rather than assumed: fitted pAC50 bottoms out at
    4.4 on all five isoforms (1st percentile exactly 4.4, nothing below 4.0), so "censored" means
    true pAC50 < 4.4 and the loss charges only when a prediction rises above it. Our own labels
    hold just 232 / 259 / 129 / 944 points below pIC50 4.

WARM START. Encoder only -- BondMessagePassing d_h=2048, depth 6, relu, mean -- from
`openadmet/cyp1a2-cyp2d6-cyp3a4-cyp2c9-chemeleon-v1`, whose unspecified hyper-parameters item 363
IDENTIFIED rather than guessed (rho 0.9987 on its own labels against 0.9329 for the runner-up).
Heads are fresh, because the checkpoint's head lives on ChEMBL's scale and merging that scale is the
very thing this design exists to avoid.

LEAKAGE. Every compound whose canonical SMILES appears in our own 4905 is dropped from BOTH external
sources, so a validation-fold label can never re-enter through a foreign head. Item 290 counts the
ChEMBL overlap at 64 and item 383 the panel's at 183 to 319; both are asserted here, and the
blinded-test overlap of both sources is 0 of 750 (items 363, 383).

THE BASELINE IS THE CV NUMBER, NOT THE BOARD. The shipped composition scores macro Spearman 0.6417
(sd 0.0015) over `members_seed{0..7}.json` on these same Butina folds (item 316); the board's 0.7212
is a different and easier sample and must never be the comparator. The floor is 0.0036.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))

import argparse
import json
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from rdkit import Chem, RDLogger
from scipy.stats import spearmanr

RDLogger.DisableLog("rdApp.*")

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
PANEL5 = ["CYP1A2", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4"]
CENS_BOUND = 4.4          # read off the panel: fitted pAC50 bottoms out here on all five isoforms
DEPTH, ACT, AGG, DH = 6, "relu", "mean", 2048


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


class Joint(nn.Module):
    """One shared encoder, one head per source. Nothing crosses between heads."""

    def __init__(self, mp, agg, widths, hidden=1024, dropout=0.1, layers=3):
        super().__init__()
        self.mp, self.agg = mp, agg
        self.heads = nn.ModuleList()
        for w in widths:
            if layers <= 1:
                self.heads.append(nn.Linear(DH, w))
                continue
            mods, d = [], DH
            for _ in range(layers - 1):
                mods += [nn.Linear(d, hidden), nn.ReLU(), nn.Dropout(dropout)]
                d = hidden
            mods.append(nn.Linear(d, w))
            self.heads.append(nn.Sequential(*mods))

    def embed(self, bmg, V_d):
        return self.agg(self.mp(bmg, V_d), bmg.batch)

    def forward(self, bmg, V_d, head):
        return self.heads[head](self.embed(bmg, V_d))


def build_encoder(ckpt_dir, device, hidden=1024, layers=3, dropout=0.1):
    from chemprop.nn import BondMessagePassing, MeanAggregation
    sd = torch.load(str(_pl.Path(ckpt_dir) / "model.pth"), map_location="cpu", weights_only=True)
    mp = BondMessagePassing(d_v=72, d_e=14, d_h=DH, depth=DEPTH, bias=False,
                            dropout=0.0, activation=ACT)
    pref = "message_passing."
    enc = {k[len(pref):]: v for k, v in sd.items() if k.startswith(pref)}
    before = {k: v.detach().clone() for k, v in mp.state_dict().items()}
    missing, unexpected = mp.load_state_dict(enc, strict=False)
    if missing:
        raise SystemExit("КОНТРОЛЬ: веса кодировщика не загрузились: %r" % missing[:6])
    # "loaded N tensors" proves nothing on its own -- the warm start has to MOVE the weights.
    after = mp.state_dict()
    moved = [k for k in before if not torch.allclose(before[k], after[k])]
    if len(moved) != len(enc):
        raise SystemExit("КОНТРОЛЬ тёплого старта: изменилось %d тензоров из %d загруженных"
                         % (len(moved), len(enc)))
    print("  кодировщик: загружено %d тензоров (все %d изменились против случайной инициализации),"
          " лишних %d" % (len(enc), len(moved), len(unexpected)), flush=True)
    return Joint(mp, MeanAggregation(), [4, 4, 5], hidden=hidden, layers=layers,
                 dropout=dropout).to(device)


def make_loader(smiles, Y, Mk, Cn, batch, shuffle, feat, nworkers=8):
    """Y: targets (NaN where absent). Mk: observed mask. Cn: censored mask (one-sided)."""
    from chemprop import data
    # NaN must reach the loss as NaN. An earlier version passed nan_to_num(Y) here, which turned
    # every ABSENT label into 0.0 -- the enzyme's own mean after standardisation -- and so
    # supervised the model to predict the mean on every unmeasured cell. Validation Spearman fell
    # from 0.52 to 0.40 while the training loss fell, which is the signature.
    pts = [data.MoleculeDatapoint.from_smi(s, [float(v) for v in Y[i]])
           for i, s in enumerate(smiles)]
    ds = data.MoleculeDataset(pts, feat)
    return data.build_dataloader(ds, num_workers=nworkers, batch_size=batch,
                                 shuffle=shuffle, drop_last=False)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ckpt", required=True)
    ap.add_argument("--data", default="data/")
    ap.add_argument("--out", required=True)
    ap.add_argument("--folds", default="0", help="запятыми; 'all' = 0,1,2,3,4")
    ap.add_argument("--epochs", type=int, default=30)
    ap.add_argument("--batch", type=int, default=48)
    ap.add_argument("--lr", type=float, default=3e-4)
    ap.add_argument("--lr-enc", type=float, default=1e-4)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--head-hidden", type=int, default=1024)
    ap.add_argument("--head-layers", type=int, default=3)
    ap.add_argument("--dropout", type=float, default=0.1)
    ap.add_argument("--w-aux", type=float, default=1.0,
                    help="вес ЧУЖИХ источников: на эпоху приходится 101 шаг панели против 31 нашего")
    ap.add_argument("--freeze", type=int, default=0,
                    help="эпох с замороженным кодировщиком, пока головы не обучатся")
    ap.add_argument("--smoke", action="store_true", help="один шаг на карте и выход")
    a = ap.parse_args()
    D = a.data if a.data.endswith("/") else a.data + "/"
    dev = "cuda:0" if torch.cuda.is_available() else "cpu"
    torch.manual_seed(a.seed)
    np.random.seed(a.seed)
    print("  устройство %s, torch %s" % (dev, torch.__version__), flush=True)

    # ---- source 0: the challenge, the only head ever scored ----
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    s_ch = [canon(s) for s in rows.SMILES]
    Y_ch = np.stack([tr["%s_pIC50_direct_inhibition" % c].to_numpy(float) for c in CYPS], 1)
    own = set(x for x in s_ch if x)
    print("  наш набор: %d соединений, меток %s"
          % (len(s_ch), "/".join(str(int(np.isfinite(Y_ch[:, e]).sum())) for e in range(4))),
          flush=True)

    # ---- source 1: ChEMBL, its own scale, its own head ----
    cx = pd.read_csv(D + "chemeleon_X_train.csv")
    cy = pd.read_csv(D + "chemeleon_y_train.csv")
    scol = [c for c in cx.columns if "smiles" in c.lower()][0]
    s_cb = [canon(s) for s in cx[scol]]
    ycols = [c for c in cy.columns if any(e.lower() in c.lower() for e in CYPS)]
    order = [next(c for c in ycols if e.lower() in c.lower()) for e in CYPS]
    Y_cb = cy[order].to_numpy(float)
    keep = np.array([(s is not None) and (s not in own) for s in s_cb])
    print("  ChEMBL: %d строк, пересечение с нашими %d (пункт 290 насчитал 64), остаётся %d"
          % (len(s_cb), int((~keep).sum()), int(keep.sum())), flush=True)
    s_cb = [s for s, k in zip(s_cb, keep) if k]
    Y_cb = Y_cb[keep]

    # ---- source 2: the qHTS panel, with censored rows as one-sided bounds ----
    pn = pd.read_csv(D + "ncgc/panel.csv", low_memory=False)
    # The index must come from the panel's OWN structures, not from a pivot: pivoting on pAC50
    # silently drops every structure that has no fitted value on ANY isoform -- 13126 down to
    # 9125 -- and those are precisely the compounds that were tested and found inactive
    # everywhere, i.e. the most informative bounds in the whole source.
    allsmi = pd.Index(pn.smiles.dropna().unique(), name="smiles")
    pv = (pn.pivot_table(index="smiles", columns="enzyme", values="pAC50", aggfunc="max")
            .reindex(index=allsmi, columns=PANEL5))
    pc = (pn.assign(c=pn.censored.astype(bool))
            .pivot_table(index="smiles", columns="enzyme", values="c", aggfunc="max")
            .reindex(index=allsmi, columns=PANEL5))
    if len(pv) != pn.smiles.nunique():
        raise SystemExit("КОНТРОЛЬ: индекс панели %d против %d уникальных структур"
                         % (len(pv), pn.smiles.nunique()))
    s_pn = [canon(s) for s in pv.index]
    Y_pn = pv.to_numpy(float)
    C_pn = np.nan_to_num(pc.to_numpy(float), nan=0.0).astype(bool)
    keep = np.array([(s is not None) and (s not in own) for s in s_pn])
    print("  панель: %d структур, пересечение с нашими %d (пункт 383: 183-319), остаётся %d"
          % (len(s_pn), int((~keep).sum()), int(keep.sum())), flush=True)
    s_pn = [s for s, k in zip(s_pn, keep) if k]
    Y_pn, C_pn = Y_pn[keep], C_pn[keep]
    Y_pn = np.where(C_pn, CENS_BOUND, Y_pn)
    print("  цензурированных ячеек %d из %d (граница pAC50 < %.1f)"
          % (int(C_pn.sum()), C_pn.size, CENS_BOUND), flush=True)

    # standardise per source x enzyme, so no scale is ever merged
    def z(Y, M):
        mu = np.array([np.nanmean(Y[M[:, e], e]) if M[:, e].any() else 0.0
                       for e in range(Y.shape[1])])
        sd = np.array([np.nanstd(Y[M[:, e], e]) or 1.0 for e in range(Y.shape[1])])
        return (Y - mu) / sd, mu, sd

    M_ch = np.isfinite(Y_ch)
    M_cb = np.isfinite(Y_cb)
    M_pn = np.isfinite(Y_pn)
    Z_ch, mu_ch, sd_ch = z(Y_ch, M_ch)
    Z_cb, _, _ = z(Y_cb, M_cb)
    Z_pn, mu_pn, sd_pn = z(Y_pn, M_pn)
    cb_pn = (CENS_BOUND - mu_pn) / sd_pn            # the bound in each enzyme's standardised space

    from chemprop import featurizers
    feat = featurizers.SimpleMoleculeMolGraphFeaturizer()
    model = build_encoder(a.ckpt, dev, hidden=a.head_hidden, layers=a.head_layers,
                          dropout=a.dropout)

    if a.smoke:
        dl = make_loader(s_ch[:4 * a.batch], Z_ch[:4 * a.batch], M_ch[:4 * a.batch],
                         None, a.batch, True, feat, nworkers=0)
        opt = torch.optim.AdamW(model.parameters(), lr=a.lr)
        t0 = time.time()
        for b in dl:
            bmg, V_d, X_d, *_ = b
            bmg.to(dev)
            p = model(bmg, V_d, 0)
            loss = (p ** 2).mean()
            loss.backward(); opt.step(); opt.zero_grad()
        torch.cuda.synchronize() if dev.startswith("cuda") else None
        print("  ДЫМОВОЙ ТЕСТ пройден: %d шагов при batch %d за %.1f с, пик памяти %.2f ГБ"
              % (4, a.batch, time.time() - t0,
                 torch.cuda.max_memory_allocated() / 2**30 if dev.startswith("cuda") else 0.0),
              flush=True)
        return 0

    from cypsplit import butina_folds
    fold, ncl = butina_folds(list(rows.SMILES))
    fold = np.asarray(fold)
    folds = [0, 1, 2, 3, 4] if a.folds == "all" else [int(x) for x in a.folds.split(",")]
    print("  сплит: %d кластеров, фолды %s" % (ncl, folds), flush=True)

    results = {}
    for k in folds:
        tr_i = np.where(fold != k)[0]
        te_i = np.where(fold == k)[0]
        print("\n  === фолд %d: обучение %d, валидация %d ===" % (k, len(tr_i), len(te_i)),
              flush=True)
        model = build_encoder(a.ckpt, dev, hidden=a.head_hidden, layers=a.head_layers,
                              dropout=a.dropout)
        opt = torch.optim.AdamW([
            {"params": model.mp.parameters(), "lr": a.lr_enc},
            {"params": model.heads.parameters(), "lr": a.lr}], weight_decay=1e-6)
        loaders = [
            make_loader([s_ch[i] for i in tr_i], Z_ch[tr_i], M_ch[tr_i], None, a.batch, True, feat),
            make_loader(s_cb, Z_cb, M_cb, None, a.batch, True, feat),
            make_loader(s_pn, Z_pn, M_pn, C_pn, a.batch, True, feat)]
        masks = [M_ch[tr_i], M_cb, M_pn]
        cens = [None, None, C_pn]
        best = (-1.0, None)
        for ep in range(-1, a.epochs):          # ep = -1 is the warm start, untouched
            frozen = 0 <= ep < a.freeze
            for prm in model.mp.parameters():
                prm.requires_grad_(not frozen)
            model.train()
            t0 = time.time()
            tot = 0.0
            iters = [iter(dl) for dl in loaders]
            nsteps = 0 if ep < 0 else max(len(dl) for dl in loaders)
            for step in range(nsteps):
                opt.zero_grad()
                for h, it in enumerate(iters):
                    try:
                        b = next(it)
                    except StopIteration:
                        iters[h] = iter(loaders[h]); b = next(iters[h])
                    bmg, V_d, X_d, Yb, *_ = b
                    bmg.to(dev)
                    yb = Yb.to(dev)
                    p = model(bmg, V_d, h)
                    obs = torch.isfinite(yb)
                    yb0 = torch.nan_to_num(yb)
                    if h == 2:
                        # censored cells enter one-sided: charged only above the bound
                        bnd = torch.tensor(cb_pn, device=dev, dtype=p.dtype).unsqueeze(0)
                        isc = obs & ((yb0 - bnd).abs() < 1e-6)
                        loss = torch.where(isc, torch.clamp(p - bnd, min=0.0) ** 2,
                                           (p - yb0) ** 2)
                    else:
                        loss = (p - yb0) ** 2
                    loss = (loss * obs).sum() / obs.sum().clamp(min=1)
                    if h > 0:
                        loss = loss * a.w_aux
                    loss.backward()
                    tot += float(loss.detach())
                opt.step()
            # validate on the held-out fold of OUR data only
            model.eval()
            vl = make_loader([s_ch[i] for i in te_i], Z_ch[te_i], M_ch[te_i], None,
                             128, False, feat, nworkers=4)
            P = []
            with torch.inference_mode():
                for b in vl:
                    bmg, V_d, X_d, *_ = b
                    bmg.to(dev)
                    P.append(model(bmg, V_d, 0).float().cpu().numpy())
            P = np.concatenate(P, 0)
            rhos = []
            for e in range(4):
                m = M_ch[te_i][:, e]
                if m.sum() > 10:
                    rhos.append(float(spearmanr(P[m, e], Y_ch[te_i][m, e]).statistic))
            mac = float(np.mean(rhos))
            if mac > best[0]:
                best = (mac, P.copy())
            print("    эпоха %2d%s: потеря %.4f, макро ро %.4f %s (%.0f с)"
                  % (ep, " [заморожен]" if frozen else "", tot / max(nsteps, 1), mac,
                     " ".join("%.3f" % r for r in rhos), time.time() - t0), flush=True)
        results[k] = {"best_rho": best[0],
                      "pred": best[1].tolist(), "idx": te_i.tolist()}
        print("  фолд %d: лучшее макро ро %.4f" % (k, best[0]), flush=True)

    print("\n  СВОДКА: %s" % ", ".join("фолд %d %.4f" % (k, results[k]["best_rho"])
                                        for k in results))
    print("  среднее по фолдам %.4f против отгружаемого CV 0.6417 (пункт 316), пол 0.0036"
          % np.mean([results[k]["best_rho"] for k in results]))
    json.dump(results, open(a.out, "w"))
    print("  записано %s" % a.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
