"""A representation LEARNED from the molecular graph --- the one axis 318 items never touched.

WHY THIS AND NOT ANOTHER MEMBER. Read off the committed board snapshot
(`results/leaderboard_2026-10-11T1926Z.json`, 273 entrants), sorting the field three ways puts us
64th by MA-ST-RAE, 71st by MA-R2 and **91st by MA-Spearman**. The board measures ordering ---
Spearman between a row's ST-RAE and its rho is -0.901 across the field, exactly one of the 63
entrants ahead of us has a worse rho than ours, and none has rho below 0.3. So item 300's headline,
that the visible numbers are calibration probes rather than models, was true on 13 September at 221
entrants and is **false now**: the probes have been pushed down and the top is models.

Which means our post-processing is not robbing us, it is flattering us. Fit-free check: 19 entrants
sit within 0.005 of our rho, their median ST-RAE is 0.6944 against our 0.5357, and the best of them
reaches 0.5069 --- so at our ordering ability we already extract more metric than all but five of
them, and the headroom left there is about 0.029 of ST-RAE, nine places. Near the top the same
comparison gives +-0.02: the leaders are not out-processing us, they are out-ordering us by 0.07 of
rho, which is twenty times the rank floor.

THE STRUCTURAL REASON, and it is visible in our own code rather than inferred. Every member sits on
ONE representation: the five are gradient boosting, pooled boosting, a GP, a ridge and `src/trunk.py`
--- and the trunk, for all its torch, is a dense net over the SAME 2295 columns. `chemprop` appears
only frozen: `src/embed*.py` run the organisers' pretrained weights and everything downstream reads
an array. Five estimator families, one representation. The 0.91--0.99 error correlation between
members (item 280) then follows mechanically, and item 288's "decorrelation comes from the ESTIMATOR
FAMILY, not from features" was measured INSIDE those 2295 columns, so it cannot speak about a
different representation. Fine-tuning was never measured here either --- the journal records it as
"low-prior per items 166/189/269/289", a prior, with the only hard evidence a competitor's single
published failure. Meanwhile `briford`, 16th at rho 0.7516 and the best readable entrant, publishes
four Chemprop D-MPNNs trained end to end on all of the data.

A grep of this file for `message passing`, `GNN`, `graph neural` and `graph network` returns zero
against controls `torch` 11, `trunk` 225 and `chemprop` 15. The axis is untouched.

WHY A HAND-WRITTEN D-MPNN AND NOT CHEMPROP. `chemprop`, `torch_geometric`, `dgl` and `deepchem` are
all absent from this environment, and `src/embed.py` says why: chemprop pulls forty-one dependencies
and the lock does not carry it. CLAUDE.md forbids raising the pins as a side effect. torch 2.13.0 and
RDKit are here, and a directed message-passing network is a hundred lines, so it is written out ---
which also keeps the whole thing reproducible under the existing pins.

THE DISCRIMINATING CONTROL. `--depth 0` runs the identical network with message passing REMOVED:
atom features are read, pooled and fed to the same head, so the parameter path is the same and only
the graph propagation is gone. If depth 0 matches depth 3, what we built is another dense net on a
different featurisation and the graph bought nothing --- the same shape of control that killed arm D
in item 339.

    uv run python verify/k118_dmpnn.py --seeds 0 --folds 0        # pilot: time one fold first
    uv run python verify/k118_dmpnn.py --seeds 0

Writes results/logs/k118_dmpnn.json and caches out-of-fold predictions under the fold digest.
Exit: 0 ran, 2 a control failed, 3 inputs missing.
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
from rdkit import Chem
from cypsplit import butina_folds, fold_digest
from submit import CYPS

import k106_probe as K
from k106_probe import score, errcorr, members, probe_oof, combine_add, load_emb

LOG = RES + "logs/k118_dmpnn.json"
CACHE = RES + "preds/k118_dmpnn_seed%d_d%d.json"
FLOOR = 0.0036          # macro-rank noise floor
HID = 300               # chemprop's default width
EPOCHS = 60
PATIENCE = 10
LR = 1e-3
BATCH = 64
INNER = 0.1             # inner validation share, drawn from the training folds only

ELEMS = [6, 7, 8, 9, 15, 16, 17, 35, 53, 5, 14, 34]     # C N O F P S Cl Br I B Si Se
HYB = [Chem.rdchem.HybridizationType.SP, Chem.rdchem.HybridizationType.SP2,
       Chem.rdchem.HybridizationType.SP3, Chem.rdchem.HybridizationType.SP3D,
       Chem.rdchem.HybridizationType.SP3D2]
BT = [Chem.rdchem.BondType.SINGLE, Chem.rdchem.BondType.DOUBLE,
      Chem.rdchem.BondType.TRIPLE, Chem.rdchem.BondType.AROMATIC]


def _onehot(v, choices):
    out = [0.0] * (len(choices) + 1)
    out[choices.index(v) if v in choices else len(choices)] = 1.0
    return out


def atom_feats(a):
    return (_onehot(a.GetAtomicNum(), ELEMS)
            + _onehot(a.GetTotalDegree(), [0, 1, 2, 3, 4, 5])
            + _onehot(a.GetFormalCharge(), [-2, -1, 0, 1, 2])
            + _onehot(a.GetHybridization(), HYB)
            + _onehot(a.GetTotalNumHs(), [0, 1, 2, 3, 4])
            + [1.0 if a.GetIsAromatic() else 0.0,
               1.0 if a.IsInRing() else 0.0,
               a.GetMass() * 0.01])


def bond_feats(b):
    return (_onehot(b.GetBondType(), BT)
            + [1.0 if b.GetIsConjugated() else 0.0,
               1.0 if b.IsInRing() else 0.0]
            + _onehot(int(b.GetStereo()), [0, 1, 2, 3]))


ADIM = len(atom_feats(Chem.MolFromSmiles("CC").GetAtomWithIdx(0)))
BDIM = len(bond_feats(Chem.MolFromSmiles("CC").GetBondWithIdx(0)))


def graph(smiles):
    """One molecule as directed bonds, in the form the message passing wants.

    Returns atom features, per-directed-bond (source atom, bond) features, the reverse-bond index
    and the bond's destination atom. A molecule with no bonds is legal and gets an empty bond set;
    the readout still works because it pools atoms.
    """
    m = Chem.MolFromSmiles(smiles)
    if m is None:
        return None
    X = np.array([atom_feats(a) for a in m.GetAtoms()], dtype=np.float32)
    src, dst, ef, rev = [], [], [], []
    for b in m.GetBonds():
        i, j, f = b.GetBeginAtomIdx(), b.GetEndAtomIdx(), bond_feats(b)
        src += [i, j]; dst += [j, i]; ef += [f, f]
        n = len(src)
        rev += [n - 1, n - 2]
    return (X, np.array(src, dtype=np.int64), np.array(dst, dtype=np.int64),
            np.array(ef, dtype=np.float32).reshape(-1, BDIM), np.array(rev, dtype=np.int64))


def batch(graphs):
    """Several molecules as one disconnected graph, with index offsets."""
    X, src, dst, ef, rev, mol = [], [], [], [], [], []
    ao = bo = 0
    for k, g in enumerate(graphs):
        x, s, d, e, r = g
        X.append(x); ef.append(e)
        src.append(s + ao); dst.append(d + ao); rev.append(r + bo)
        mol.append(np.full(len(x), k, dtype=np.int64))
        ao += len(x); bo += len(s)
    def cat(v, width=None):
        # Always concatenate. Testing `len(v[0])` would look at the FIRST molecule only, so a
        # batch whose first molecule has no bonds would silently lose every other molecule's.
        if not v:
            return np.zeros((0, width) if width else (0,),
                            dtype=np.float32 if width else np.int64)
        return np.concatenate(v)
    return (torch.from_numpy(np.concatenate(X)),
            torch.from_numpy(cat(src)), torch.from_numpy(cat(dst)),
            torch.from_numpy(cat(ef, BDIM)), torch.from_numpy(cat(rev)),
            torch.from_numpy(np.concatenate(mol)), len(graphs))


class DMPNN(nn.Module):
    """Directed message passing over bonds, then a shared head over four enzymes.

    `depth = 0` is the discriminating control: the bond network is never called, so the molecule
    representation is a pooled linear map of the ATOM features alone. Every other tensor, the head
    included, is identical, which is what makes the comparison about the graph and not the capacity.
    """

    def __init__(self, depth=3, hid=HID, nout=4):
        super().__init__()
        self.depth = depth
        self.Wi = nn.Linear(ADIM + BDIM, hid, bias=False)
        self.Wh = nn.Linear(hid, hid, bias=False)
        self.Wo = nn.Linear(ADIM + hid, hid)
        self.head = nn.Sequential(nn.ReLU(), nn.Dropout(0.1), nn.Linear(hid, hid),
                                  nn.ReLU(), nn.Linear(hid, nout))
        self.act = nn.ReLU()

    def forward(self, X, src, dst, ef, rev, mol, n):
        hid = self.Wo.out_features
        if self.depth > 0 and len(src):
            h0 = self.act(self.Wi(torch.cat([X[src], ef], 1)))
            h = h0
            for _ in range(self.depth):
                # sum of incoming messages at each atom, minus the reverse bond: the "directed"
                # part of a D-MPNN, which is what stops a message walking straight back.
                agg = torch.zeros(len(X), hid).index_add_(0, dst, h)
                h = self.act(h0 + self.Wh(agg[src] - h[rev]))
            m = torch.zeros(len(X), hid).index_add_(0, dst, h)
        else:
            m = torch.zeros(len(X), hid)
        hv = self.act(self.Wo(torch.cat([X, m], 1)))
        pooled = torch.zeros(n, hid).index_add_(0, mol, hv)
        cnt = torch.zeros(n).index_add_(0, mol, torch.ones(len(mol))).clamp(min=1)
        return self.head(pooled / cnt.unsqueeze(1))


def train_fold(gs, Y, M, trn, te, depth, seed, log):
    """Train on `trn`, predict `te`. Targets standardised on the TRAINING rows only."""
    torch.manual_seed(seed)
    rng = np.random.default_rng(seed)
    idx = np.where(trn)[0]
    rng.shuffle(idx)
    ncut = max(1, int(len(idx) * INNER))
    inner, outer = idx[:ncut], idx[ncut:]

    mu = np.array([Y[outer][M[outer][:, e], e].mean() for e in range(len(CYPS))])
    sd = np.array([Y[outer][M[outer][:, e], e].std() for e in range(len(CYPS))])
    Z = (Y - mu) / sd

    net = DMPNN(depth=depth)
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
            with torch.set_grad_enabled(train):
                p = net(*bt)
                loss = (((p - tz) ** 2) * tm).sum() / tm.sum().clamp(min=1)
            if train:
                opt.zero_grad(); loss.backward(); opt.step()
            tot += float(loss.detach()) * float(tm.sum()); nseen += float(tm.sum())
        return tot / max(nseen, 1)

    curve = []
    for ep in range(EPOCHS):
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
    log.append({"depth": depth, "epochs": len(curve), "loss_first": curve[0],
                "loss_best": round(best, 4), "curve_tail": curve[-3:]})

    net.eval()
    out = np.zeros((int(te.sum()), len(CYPS)))
    rows = np.where(te)[0]
    with torch.no_grad():
        for b in range(0, len(rows), BATCH):
            sel = rows[b:b + BATCH]
            out[b:b + len(sel)] = net(*batch([gs[i] for i in sel])).numpy()
    return out * sd + mu, curve[0][0], best


def oof(gs, Y, M, fold, depth, seed, folds, log):
    P = np.full((len(gs), len(CYPS)), np.nan)
    first, bests = [], []
    for f in folds:
        trn, te = fold != f, fold == f
        if te.sum() == 0:
            continue
        t0 = time.time()
        pred, l0, lb = train_fold(gs, Y, M, trn, te, depth, seed * 100 + f, log)
        P[te] = pred
        first.append(l0); bests.append(lb)
        print("      глубина %d фолд %d: %.0f с, потеря %.4f -> %.4f"
              % (depth, f, time.time() - t0, l0, lb), flush=True)
    return P, float(np.mean(first)), float(np.mean(bests))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--folds", default="0,1,2,3,4")
    ap.add_argument("--depths", default="3,0")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]
    folds_want = [int(x) for x in a.folds.split(",") if x.strip()]
    depths = [int(x) for x in a.depths.split(",") if x.strip()]

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

    print("  строю графы из SMILES (признаков: атом %d, связь %d)" % (ADIM, BDIM), flush=True)
    gs = [graph(s) for s in rows.SMILES]
    bad = [i for i, g in enumerate(gs) if g is None]
    if bad:
        print("ОТКАЗ: RDKit не разобрал %d SMILES из rows.csv, чего быть не может" % len(bad))
        return 2
    nat = np.array([len(g[0]) for g in gs])
    nbd = np.array([len(g[1]) for g in gs])
    print("      молекул %d, атомов медиана %d (макс %d), направленных связей медиана %d"
          % (len(gs), int(np.median(nat)), nat.max(), int(np.median(nbd))), flush=True)

    out = json.load(open(LOG)) if _pl.Path(LOG).exists() else {}
    out.setdefault("meta", {}).update(hid=HID, epochs=EPOCHS, patience=PATIENCE, lr=LR,
                                      batch=BATCH, inner=INNER, adim=ADIM, bdim=BDIM)
    out.setdefault("seeds", {})

    for s in seeds:
        fold, ncl = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        print("\n  сид %d: дайджест %s, кластеров %d" % (s, dg, ncl), flush=True)
        if s == 0 and dg != "2d93c19815e14261":
            print("ОТКАЗ: разбиение сдвинулось, сид 0 дал %s" % dg)
            return 2
        full = sorted(folds_want) == [0, 1, 2, 3, 4]
        res = out["seeds"].get(str(s), {})

        for depth in depths:
            path = CACHE % (s, depth)
            store = json.load(open(path)) if _pl.Path(path).exists() else {}
            if store.get("digest") == dg and store.get("full") and full:
                P = np.array(store["P"], dtype=float)
                l0, lb, curves = store["loss_first"], store["loss_best"], store.get("curves", [])
                print("      глубина %d из кэша" % depth, flush=True)
            else:
                curves = []
                P, l0, lb = oof(gs, Yf, M, fold, depth, s, sorted(folds_want), curves)
                if full:
                    json.dump({"digest": dg, "full": True, "loss_first": l0, "loss_best": lb,
                               "curves": curves,
                               "P": [[None if not np.isfinite(x) else float(x) for x in r]
                                     for r in P]}, open(path, "w"))

            # CONTROL that must succeed: training moved. A net that never learned would make the
            # comparison with depth 0 meaningless in the same way a dead grep does.
            drop = (l0 - lb) / max(l0, 1e-9)
            print("      глубина %d: средняя потеря %.4f -> %.4f (падение %.0f %%)"
                  % (depth, l0, lb, 100 * drop), flush=True)
            if drop < 0.05:
                print("ОТКАЗ: потеря почти не упала --- сеть не обучилась, сравнивать нечего")
                return 2

            cols = [P[M[:, e], e] for e in range(len(CYPS))]
            if not all(np.isfinite(c).all() for c in cols):
                print("      глубина %d: не все фолды посчитаны, пропускаю оценку" % depth,
                      flush=True)
                res[f"depth{depth}"] = {"loss_first": l0, "loss_best": lb, "partial": True}
                continue
            sc = score(cols, Y, LO, HI, M, fold)
            res[f"depth{depth}"] = {"solo": sc, "loss_first": l0, "loss_best": lb,
                                    "per_enzyme_rho": sc[2]}
            print("      глубина %d: макро ранг %.4f, пара %.4f | по ферментам %s"
                  % (depth, sc[0], sc[1], " ".join("%.4f" % x for x in sc[2])), flush=True)

        if full and all(f"depth{d}" in res and not res[f"depth{d}"].get("partial")
                        for d in depths) and 3 in depths:
            pdz = members(s, XFULL, Y, M, fold, LO, HI)
            probe = [probe_oof(load_emb(rows), Y, M, fold, e)[0] for e in range(len(CYPS))]
            six = combine_add(pdz, probe)
            s_six = score(six, Y, LO, HI, M, fold)
            P3 = np.array(json.load(open(CACHE % (s, 3)))["P"], dtype=float)
            cols3 = [P3[M[:, e], e] for e in range(len(CYPS))]
            pm, vm, vs, cal = errcorr(cols3, pdz, Y, M)
            seven = [np.mean([Q[e] for k, Q in pdz if K._keep(e, k)] + [probe[e], cols3[e]], 0)
                     for e in range(len(CYPS))]
            s7 = score(seven, Y, LO, HI, M, fold)
            res["ensemble"] = {"six": s_six, "seven": s7,
                               "gain_rank": s7[0] - s_six[0], "gain_pair": s7[1] - s_six[1],
                               "errcorr_vs_five": vm, "errcorr_per_member": pm,
                               "member_band": {k: float(np.mean(v)) for k, v in cal.items()}}
            print("      состав: шесть %.4f -> семь %.4f (%+.4f ранга), пара %+.4f  пол %.4f"
                  % (s_six[0], s7[0], s7[0] - s_six[0], s7[1] - s_six[1], FLOOR), flush=True)
            print("      corr ошибок со средним пяти: %s   (полоса членов %.3f--%.3f)"
                  % (" ".join("%.3f" % v for v in vm),
                     min(np.mean(v) for v in cal.values()),
                     max(np.mean(v) for v in cal.values())), flush=True)

        out["seeds"][str(s)] = res
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
        print("      записано: %s" % LOG, flush=True)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
