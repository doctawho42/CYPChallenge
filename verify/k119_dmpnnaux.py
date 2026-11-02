"""The graph net is overfitting, not under-trained --- so give it the screen as a second head.

WHAT ITEM 341 LEFT. The hand-written D-MPNN reaches macro rank 0.5778 against the shipped
composition's 0.6568, and its errors are the first in this project to fall outside the member band
(0.894 against 0.966--0.982). Right axis, insufficient vehicle. The question is which of the three
obvious strengthenings is the one worth paying for, and the training curves answer it before any of
them is built:

    фолд   эпох из 60   обучающая   валидационная
    0              25        0.395           0.715
    1              30        0.284           0.697
    2              29        0.406           0.642
    3              25        0.392           0.651
    4              25        0.401           0.738

Early stopping fires at 25--30 of 60 with the training loss at half the validation loss. **More
epochs cannot help** --- validation is already turning up. The binding constraint is supervision, not
training time, which is exactly the channel `src/trunk.py` was built around: the single-point screen
carries a reading for 2570--3090 compounds per enzyme whose curve was never run, so the supervision
available to a shared representation goes up 2.76-fold. The trunk measured that channel as real on a
dense net over 2295 columns; a graph net trained on 4905 molecules is the kind of model that should
want it far more.

THE CONTROL IS THE TRUNK'S OWN, and it is better than omitting the head. `lam = 0` keeps the
screening head in the network and gives it no gradient, so the parameter count and the
initialisation order are identical between arms and the loss term is the only difference. Omitting
the head instead would change the random draw and confound the comparison with a different
initialisation --- which is also why `DMPNNAux` subclasses item 341's `DMPNN` and adds its head
AFTER `super().__init__()`: the parent's parameters draw in the same order, so item 341's cached
numbers stay reproducible.

    uv run python verify/k119_dmpnnaux.py --lam 0 --nets 1      # контроль
    uv run python verify/k119_dmpnnaux.py --lam 1 --nets 1      # арм
    uv run python verify/k119_dmpnnaux.py --lam 1 --nets 4      # и усреднение, если арм выиграл

Writes results/logs/k119_dmpnnaux.json. Exit: 0 ran, 2 a control failed, 3 inputs missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn

from cypsplit import butina_folds, fold_digest
from submit import CYPS, load_screen

import k106_probe as K
from k106_probe import score, errcorr, members, probe_oof, combine_add, load_emb
from k118_dmpnn import DMPNN, batch, graph, BATCH, EPOCHS, HID, INNER, LR, PATIENCE

LOG = RES + "logs/k119_dmpnnaux.json"
CACHE = RES + "preds/k119_aux_seed%d_lam%s_n%d.json"
FLOOR = 0.0036


class DMPNNAux(DMPNN):
    """Item 341's network with a second head on the single-point screen.

    The head is created AFTER `super().__init__()` so the parent's parameters draw from the
    generator in the same order they did in item 341; the baseline therefore still reproduces.
    """

    def __init__(self, depth=3, hid=HID, nout=4):
        super().__init__(depth=depth, hid=hid, nout=nout)
        self.head_scr = nn.Sequential(nn.ReLU(), nn.Dropout(0.1), nn.Linear(hid, hid),
                                      nn.ReLU(), nn.Linear(hid, nout))

    def both(self, X, src, dst, ef, rev, mol, n):
        hid = self.Wo.out_features
        if self.depth > 0 and len(src):
            h0 = self.act(self.Wi(torch.cat([X[src], ef], 1)))
            h = h0
            for _ in range(self.depth):
                agg = torch.zeros(len(X), hid).index_add_(0, dst, h)
                h = self.act(h0 + self.Wh(agg[src] - h[rev]))
            m = torch.zeros(len(X), hid).index_add_(0, dst, h)
        else:
            m = torch.zeros(len(X), hid)
        hv = self.act(self.Wo(torch.cat([X, m], 1)))
        pooled = torch.zeros(n, hid).index_add_(0, mol, hv)
        cnt = torch.zeros(n).index_add_(0, mol, torch.ones(len(mol))).clamp(min=1)
        z = pooled / cnt.unsqueeze(1)
        return self.head(z), self.head_scr(z)


def train_fold(gs, Y, M, S, SM, trn, te, lam, seed, log):
    """Train one net on `trn`, predict `te`. Both targets standardised on TRAINING rows only."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    idx = np.where(trn)[0]
    rng.shuffle(idx)
    ncut = max(1, int(len(idx) * INNER))
    inner, outer = idx[:ncut], idx[ncut:]

    mu = np.array([Y[outer][M[outer][:, e], e].mean() for e in range(len(CYPS))])
    sd = np.array([Y[outer][M[outer][:, e], e].std() for e in range(len(CYPS))])
    Z = (Y - mu) / sd
    smu = np.array([S[outer][SM[outer][:, e], e].mean() if SM[outer][:, e].any() else 0.0
                    for e in range(len(CYPS))])
    ssd = np.array([S[outer][SM[outer][:, e], e].std() if SM[outer][:, e].any() else 1.0
                    for e in range(len(CYPS))])
    ZS = (S - smu) / np.where(ssd > 1e-9, ssd, 1.0)

    net = DMPNNAux(depth=3)
    opt = torch.optim.Adam(net.parameters(), lr=LR)
    best, bad, bstate = np.inf, 0, None

    def run(rows, train):
        net.train(train)
        tot = nseen = 0.0
        order = rng.permutation(rows) if train else rows
        for b in range(0, len(order), BATCH):
            sel = order[b:b + BATCH]
            bt = batch([gs[i] for i in sel])
            tz = torch.from_numpy(Z[sel].astype(np.float32))
            tm = torch.from_numpy(M[sel].astype(np.float32))
            ts = torch.from_numpy(ZS[sel].astype(np.float32))
            tsm = torch.from_numpy(SM[sel].astype(np.float32))
            with torch.set_grad_enabled(train):
                p, q = net.both(*bt)
                main = (((p - tz) ** 2) * tm).sum() / tm.sum().clamp(min=1)
                loss = main
                if lam > 0:
                    aux = (((q - ts) ** 2) * tsm).sum() / tsm.sum().clamp(min=1)
                    loss = main + lam * aux
            if train:
                opt.zero_grad(); loss.backward(); opt.step()
            # The selection and the reported number are the MAIN loss in both arms: an arm with an
            # extra term would otherwise early-stop on a different quantity and the comparison
            # would be between two different stopping rules rather than between two models.
            tot += float(main.detach()) * float(tm.sum()); nseen += float(tm.sum())
        return tot / max(nseen, 1)

    curve = []
    for _ in range(EPOCHS):
        tr_loss = run(outer, True)
        va_loss = run(inner, False)
        curve.append((round(tr_loss, 4), round(va_loss, 4)))
        if va_loss < best - 1e-4:
            best, bad = va_loss, 0
            bstate = {k: v.detach().clone() for k, v in net.state_dict().items()}
        else:
            bad += 1
            if bad >= PATIENCE:
                break
    if bstate is not None:
        net.load_state_dict(bstate)
    log.append({"lam": lam, "epochs": len(curve), "loss_first": curve[0],
                "loss_best": round(best, 4), "curve_tail": curve[-3:]})

    net.eval()
    out = np.zeros((int(te.sum()), len(CYPS)))
    rows = np.where(te)[0]
    with torch.no_grad():
        for b in range(0, len(rows), BATCH):
            sel = rows[b:b + BATCH]
            out[b:b + len(sel)] = net.both(*batch([gs[i] for i in sel]))[0].numpy()
    return out * sd + mu, curve[0][0], best


def oof(gs, Y, M, S, SM, fold, lam, seed, nets, log):
    """Out of fold, averaged over `nets` initialisations --- briford's axis, four nets averaged."""
    acc = np.zeros((len(gs), len(CYPS)))
    first, bests = [], []
    for k in range(nets):
        P = np.full((len(gs), len(CYPS)), np.nan)
        for f in range(5):
            trn, te = fold != f, fold == f
            if te.sum() == 0:
                continue
            t0 = time.time()
            pred, l0, lb = train_fold(gs, Y, M, S, SM, trn, te, lam,
                                      seed * 1000 + k * 10 + f, log)
            P[te] = pred
            first.append(l0); bests.append(lb)
            print("      сеть %d фолд %d: %.0f с, потеря %.4f -> %.4f"
                  % (k, f, time.time() - t0, l0, lb), flush=True)
        acc += P
    return acc / nets, float(np.mean(first)), float(np.mean(bests))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--lam", default="1")
    ap.add_argument("--nets", type=int, default=1)
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]
    lam = float(a.lam)

    rows = pd.read_csv(D + "rows.csv")
    z = np.load(D + "feats.npz")
    XFULL = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)
    Yf = np.nan_to_num(Y)
    K.S.LO, K.S.HI = LO, HI

    L2, has = load_screen(rows)
    S, SM = np.nan_to_num(L2), has & ~M          # only where no curve exists for that enzyme
    print("  строю графы", flush=True)
    gs = [graph(s) for s in rows.SMILES]
    if any(g is None for g in gs):
        print("ОТКАЗ: RDKit не разобрал SMILES из rows.csv")
        return 2

    # count first: the supervision this arm actually adds
    add = SM.sum(0)
    cur = M.sum(0)
    print("  добавочная супервизия: " + "  ".join(
        "%s +%d к %d (x%.2f)" % (c[3:], add[e], cur[e], (cur[e] + add[e]) / cur[e])
        for e, c in enumerate(CYPS)), flush=True)
    print("  итого x%.2f  (src/trunk.py: x2.76)"
          % ((cur.sum() + add.sum()) / cur.sum()), flush=True)
    if add.sum() == 0:
        print("ОТКАЗ: скрининг не добавляет ни одной метки")
        return 2

    out = json.load(open(LOG)) if _pl.Path(LOG).exists() else {}
    out.setdefault("meta", {}).update(hid=HID, epochs=EPOCHS, patience=PATIENCE, lr=LR,
                                      batch=BATCH, inner=INNER,
                                      added_labels={c: int(add[e]) for e, c in enumerate(CYPS)})
    out.setdefault("arms", {})

    for s in seeds:
        fold, ncl = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        print("\n  сид %d: дайджест %s, кластеров %d, lam %s, сетей %d"
              % (s, dg, ncl, a.lam, a.nets), flush=True)
        if s == 0 and dg != "2d93c19815e14261":
            print("ОТКАЗ: разбиение сдвинулось, сид 0 дал %s" % dg)
            return 2

        path = CACHE % (s, a.lam, a.nets)
        store = json.load(open(path)) if _pl.Path(path).exists() else {}
        curves = []
        if store.get("digest") == dg:
            P = np.array(store["P"], dtype=float)
            l0, lb, curves = store["loss_first"], store["loss_best"], store.get("curves", [])
            print("      из кэша", flush=True)
        else:
            P, l0, lb = oof(gs, Yf, M, S, SM, fold, lam, s, a.nets, curves)
            json.dump({"digest": dg, "loss_first": l0, "loss_best": lb, "curves": curves,
                       "P": [[float(x) for x in r] for r in P]}, open(path, "w"))

        drop = (l0 - lb) / max(l0, 1e-9)
        print("      основная потеря %.4f -> %.4f (падение %.0f %%)" % (l0, lb, 100 * drop),
              flush=True)
        if drop < 0.05:
            print("ОТКАЗ: потеря почти не упала --- сеть не обучилась")
            return 2

        cols = [P[M[:, e], e] for e in range(len(CYPS))]
        sc = score(cols, Y, LO, HI, M, fold)
        pdz = members(s, XFULL, Y, M, fold, LO, HI)
        probe = [probe_oof(load_emb(rows), Y, M, fold, e)[0] for e in range(len(CYPS))]
        six = combine_add(pdz, probe)
        s_six = score(six, Y, LO, HI, M, fold)
        pm, vm, vs, cal = errcorr(cols, pdz, Y, M)
        seven = [np.mean([Q[e] for k, Q in pdz if K._keep(e, k)] + [probe[e], cols[e]], 0)
                 for e in range(len(CYPS))]
        s7 = score(seven, Y, LO, HI, M, fold)

        key = "lam%s_n%d_seed%d" % (a.lam, a.nets, s)
        out["arms"][key] = {"solo": sc, "six": s_six, "seven": s7,
                            "gain_rank": s7[0] - s_six[0], "gain_pair": s7[1] - s_six[1],
                            "errcorr_vs_five": vm, "loss_first": l0, "loss_best": lb,
                            "epochs": [c["epochs"] for c in curves] or None}
        print("      одиночный макро ранг %.4f, пара %.4f | по ферментам %s"
              % (sc[0], sc[1], " ".join("%.4f" % x for x in sc[2])), flush=True)
        print("      седьмым членом: %.4f -> %.4f (%+.4f ранга), пара %+.4f  пол %.4f"
              % (s_six[0], s7[0], s7[0] - s_six[0], s7[1] - s_six[1], FLOOR), flush=True)
        print("      corr ошибок со средним пяти: %s   (полоса членов %.3f--%.3f)"
              % (" ".join("%.3f" % v for v in vm),
                 min(np.mean(v) for v in cal.values()),
                 max(np.mean(v) for v in cal.values())), flush=True)

        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
        print("      записано: %s" % LOG, flush=True)

    print("\n" + "=" * 88)
    print("СКРИНИНГ ВТОРОЙ ГОЛОВОЙ НА ГРАФОВОЙ СЕТИ (пункт 341 дал 0.5778 без неё)")
    print("=" * 88)
    print("  %-22s %10s %9s %10s %9s" % ("арм", "одиночн.", "пара", "Δ 7-м", "corr"))
    for k, v in sorted(out["arms"].items()):
        print("  %-22s %10.4f %9.4f %+10.4f %9.3f"
              % (k, v["solo"][0], v["solo"][1], v["gain_rank"],
                 float(np.mean(v["errcorr_vs_five"]))))
    print("  %-22s %10.4f %9s %10s %9s" % ("пункт 341, без головы", 0.5778, "0.7139", "-0.0048",
                                           "0.894"))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
