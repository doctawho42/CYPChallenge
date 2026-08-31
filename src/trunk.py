"""Shared trunk, two heads, and lambda_scr as the switch.

The question this answers: does bringing the primary screen in as a second likelihood
term - not as an input feature - help the pIC50 head?

Why not as an input feature. The blind test set carries Molecule_Name and SMILES and
nothing else; there is no screening reading for those 750 compounds, and a model that
conditions on one at inference cannot be submitted. src/mech3.py already went down the
neighbouring road, feeding a *predicted* screen in as a column, and measured exactly zero.
A shared latent layer is a different object: the gradient from the log2fc head shapes the
representation, and at inference the screening head is simply not called.

What the channel actually buys, measured before writing any of this: per enzyme the
screen carries a reading for 2570-3090 compounds whose curve was never run, so the
supervision available to the trunk goes up by a factor of 2.1 to 3.4 depending on the
enzyme, 2.76 overall. Note it is not extra *rows* - every one of the 4905 molecules has
at least one pIC50 label - it is extra *columns* for rows already present.

The control is the same architecture with lambda_scr = 0. Both heads exist in both arms,
so the parameter count is identical; with lambda_scr = 0 the screening head simply
receives no gradient and cannot influence the trunk. Weight initialisation is seeded on
(split seed, fold) and not on lambda, so the two arms start from identical weights and the
loss term is the only difference between them.

Nothing here touches cypsplit.py, results/preds/oof.json or the golden digest.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from scipy.stats import spearmanr

from cypsplit import butina_folds

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]

# Instrument calibration fitted per enzyme by verify/g1_calib.py. The document's claim is
# that the screen constrains pi *through* this map, not through a free second head: it is
# a fixed monotone saturating function, and knowing where the curve saturates is what
# turns rank into scale. Zero new parameters.
PC0 = 4.305
CAL_E = np.array([0.728, 0.621, 0.867, 0.931], dtype=np.float32)
CAL_H = np.array([1.261, 1.112, 1.243, 1.968], dtype=np.float32)

# Chosen on the lambda_scr = 0 arm alone, never looking at an arm with the screening
# term on, so the choice cannot favour the experiment. It moves the absolute level; it
# cannot move the comparison, which is between two arms sharing these settings exactly.
#
# The criterion was stated before searching: get the trunk into the same league as the
# gradient boosting the rest of the repository uses, because a null result on a trunk
# that is far weaker would say nothing about the channel. It reaches macro ST-RAE 0.7672
# against the boosting's 0.7673.
#
# The one non-obvious choice is BLOCKS. On the full FP+DESC+MECH matrix the best this
# trunk managed was 0.7949 - 2048 fingerprint bits against roughly 1400 labelled rows per
# enzyme is more width than an MLP can use, and regularising harder made it worse, not
# better (weight decay 1e-3 gave 0.8264). Dropping the fingerprint and keeping the 247
# descriptor and mechanistic columns recovered the whole gap at once.
BLOCKS = "DESC+MECH"
HIDDEN = 512
DEPTH = 2
DROPOUT = 0.5
LR = 1e-3
WEIGHT_DECAY = 1e-4
EPOCHS = 200
BATCH = 256


def load(blocks="FP+DESC+MECH"):
    """Feature matrix, both target blocks, and their masks, all in rows.csv order.

    `blocks` selects which feature blocks enter the trunk. 2048 fingerprint bits against
    roughly 1400 labelled rows per enzyme is a lot of width for an MLP; the narrower
    DESC+MECH matrix is better conditioned and may be the more sensitive instrument for
    the question actually being asked. Whatever is chosen applies to both arms.
    """
    z = np.load(D + "feats.npz")
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())

    # Morgan counts are heavily zero-inflated with a long tail; log1p first so that
    # standardising them afterwards does not hand a few busy bits enormous weight.
    parts = {"FP": np.log1p(z["FP"]), "DESC": z["DESC"], "MECH": z["MECH"]}
    X = np.hstack([parts[b] for b in blocks.split("+")]).astype(np.float32)
    X = np.nan_to_num(X, nan=0.0, posinf=0.0, neginf=0.0)

    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    piv = sc.pivot_table(index="Molecule_Name", columns="enzyme", values="log2fc_estimate")
    scr = rows.set_index("Molecule_Name").join(piv)[CYPS].to_numpy(np.float32)

    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(np.float32) for c in CYPS], 1)
    lo = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(np.float32) for c in CYPS], 1)
    hi = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(np.float32) for c in CYPS], 1)
    return X, y, lo, hi, scr, list(rows.SMILES)


class Net(nn.Module):
    """Shared trunk, one head per target block. Both heads always exist."""

    def __init__(self, d_in, hidden=None, depth=None, dropout=None, n_out=4):
        hidden = HIDDEN if hidden is None else hidden
        depth = DEPTH if depth is None else depth
        dropout = DROPOUT if dropout is None else dropout
        super().__init__()
        layers, d = [], d_in
        for _ in range(depth):
            layers += [nn.Linear(d, hidden), nn.ReLU(), nn.Dropout(dropout)]
            d = hidden
        self.trunk = nn.Sequential(*layers)
        self.head_pic = nn.Linear(d, n_out)
        self.head_scr = nn.Linear(d, n_out)

    def forward(self, x):
        h = self.trunk(x)
        return self.head_pic(h), self.head_scr(h)


def g_of_pi(pi, e_, h_):
    """log2fc predicted from pIC50 through the fitted instrument calibration.

    g(pi) = log2(1 - E / (1 + 10^{h (pC0 - pi)})). Differentiable in pi, so the screening
    residual sends gradient back into the same latent the pIC50 head reads - which is what
    "shared latent curve" in the document actually means, as opposed to two free heads
    that are under no obligation to agree about anything.
    """
    # The exponent has to be bounded. Unbounded, 10^(h(pC0 - pi)) overflows to inf as soon
    # as the predicted potency wanders far below pC0, and although the forward value stays
    # finite (e/(1+inf) = 0), the gradient through the overflow is nan and the run dies.
    # This is not hypothetical: the noise ladder lost the calibrated arm on split seed 0 at
    # eta = 0.5 and eta = 1 to exactly this. Clamping at +-30 changes nothing that was
    # already finite - 10^30 and inf give the same inhibited fraction to float precision -
    # and only replaces nan gradients with finite ones.
    inh = e_ / (1.0 + torch.pow(10.0, torch.clamp(h_ * (PC0 - pi), min=-30.0, max=30.0)))
    return torch.log2(torch.clamp(1.0 - inh, min=1e-3))


def masked_mse(pred, target, mask):
    """Mean squared error over observed cells only; zero if nothing is observed.

    Normalised by the number of observed cells rather than summed, so that lambda_scr
    means the same thing regardless of how dense each block happens to be.
    """
    if mask.sum() == 0:
        return pred.sum() * 0.0
    d = (pred - target)[mask]
    return (d * d).mean()


def run_fold(X, y, scr, fold, f, lam, seed, device, mode="twohead", zn=None, eta=0.0):
    """Train one fold and return held-out pIC50 predictions in original units."""
    te = fold == f
    trn = ~te

    my = ~np.isnan(y)
    ms = ~np.isnan(scr)

    # Every statistic below is fitted on the training folds only.
    xm = X[trn].mean(0)
    xs = X[trn].std(0)
    xs[xs < 1e-6] = 1.0
    Xn = (X - xm) / xs

    ym = np.zeros(4, np.float32); ys = np.ones(4, np.float32)
    sm = np.zeros(4, np.float32); ss = np.ones(4, np.float32)
    for e in range(4):
        a = my[:, e] & trn
        if a.sum() > 1:
            ym[e], ys[e] = y[a, e].mean(), max(y[a, e].std(), 1e-6)
        b = ms[:, e] & trn
        if b.sum() > 1:
            sm[e], ss[e] = scr[b, e].mean(), max(scr[b, e].std(), 1e-6)
    yn = (np.nan_to_num(y) - ym) / ys
    sn = (np.nan_to_num(scr) - sm) / ss

    # Degrading the screening channel by a known amount. The division by sqrt(1 + eta^2)
    # is the whole point: without it the noise inflates the target's sd and lambda would
    # silently change meaning, so "less information" would be confounded with "more
    # weight". With it the standardised target keeps sd 1 and only its correlation with
    # pIC50 falls, by exactly 1 / sqrt(1 + eta^2). The noise draw is shared across lambda
    # and across modes, so arms compared at one eta see the same corrupted channel.
    if eta > 0 and zn is not None:
        sn = (sn + eta * zn) / np.sqrt(1.0 + eta * eta)

    t = lambda a: torch.as_tensor(a, device=device)
    Xt, yt, st = t(Xn), t(yn), t(sn)
    myt, mst = t(my), t(ms)
    ymt, yst = t(ym), t(ys)
    smt, sst = t(sm), t(ss)
    cal_e, cal_h = t(CAL_E), t(CAL_H)

    # Seeded on (seed, fold) and NOT on lam, so both arms start from the same weights.
    torch.manual_seed(hash((seed, f)) % (2 ** 31))
    net = Net(X.shape[1]).to(device)
    opt = torch.optim.Adam(net.parameters(), lr=LR, weight_decay=WEIGHT_DECAY)

    idx = np.where(trn)[0]
    g = np.random.default_rng(1000 + seed * 10 + f)
    net.train()
    for _ in range(EPOCHS):
        for b in np.array_split(g.permutation(idx), max(1, len(idx) // BATCH)):
            bt = t(b)
            p, s = net(Xt[bt])
            loss = masked_mse(p, yt[bt], myt[bt])
            if lam > 0:
                if mode == "twohead":
                    loss = loss + lam * masked_mse(s, st[bt], mst[bt])
                else:
                    # One latent: the screening prediction is g(pi_hat), not a free head.
                    # Standardised on the same statistics as the observed screen so that
                    # lambda keeps the same meaning across modes.
                    pi = p * yst + ymt
                    loss = loss + lam * masked_mse(
                        (g_of_pi(pi, cal_e, cal_h) - smt) / sst, st[bt], mst[bt])
            opt.zero_grad(); loss.backward(); opt.step()

    net.eval()
    with torch.no_grad():
        p, _ = net(Xt[t(np.where(te)[0])])
    return p.cpu().numpy() * ys + ym


def fit_predict_test(FP_te, DESC_te, MECH_te, lam=3.0, seed=0, mode="twohead",
                     blocks=None, device="cpu"):
    """Train on the whole training set and predict the blinded test rows.

    Item 120 measured the trunk as a fifth ensemble member at -0.0061 of pair and +0.0054 of
    rank, four seeds out of four on every enzyme, which is more than twice what the ridge
    contributes and the only member that helps all four. Acting on that needs the trunk applied
    to the 750 test structures, and this is that path.

    Nothing about the fitting is new. The test rows are appended with their own fold index and
    with NaN in both target blocks, so `masked_mse` gives them no gradient and every statistic
    `run_fold` computes -- the feature standardisation included -- is still taken over training
    rows alone. `run_fold` itself is called unchanged, which is what makes the regression check
    below possible.

    The one thing a caller must not get wrong: `load()` puts the fingerprint block through
    log1p before standardising, and `src/submit.py` does not. The transform is applied here so
    that a caller passing raw counts, as submit.py builds them, gets the right matrix.
    """
    blocks = blocks or BLOCKS
    X, y, lo, hi, scr, smiles = load(blocks)
    parts_te = {"FP": np.log1p(np.asarray(FP_te, np.float32)),
                "DESC": np.asarray(DESC_te, np.float32),
                "MECH": np.asarray(MECH_te, np.float32)}
    Xte = np.hstack([parts_te[b] for b in blocks.split("+")]).astype(np.float32)
    Xte = np.nan_to_num(Xte, nan=0.0, posinf=0.0, neginf=0.0)
    if Xte.shape[1] != X.shape[1]:
        raise SystemExit(f"ширина не совпала: обучение {X.shape[1]}, тест {Xte.shape[1]}")

    n_tr, n_te = len(X), len(Xte)
    Xa = np.vstack([X, Xte])
    ya = np.vstack([y, np.full((n_te, 4), np.nan, np.float32)])
    sa = np.vstack([scr, np.full((n_te, 4), np.nan, np.float32)])
    fold = np.concatenate([np.zeros(n_tr, int), np.ones(n_te, int)])
    return run_fold(Xa, ya, sa, fold, 1, lam, seed, device, mode=mode)


def check_test_path(seed=0, lam=3.0, mode="twohead", blocks=None, device="cpu"):
    """The test path must be `run_fold` and nothing else. Prove it rather than assert it.

    Holding out fold 1 through the ordinary route and through the appended-rows route has to
    give identical numbers: same weight seed (both use f = 1), same batching seed, same
    standardisation over the complement of fold 1. Any difference means the appended rows are
    reaching a statistic they should not.
    """
    blocks = blocks or BLOCKS
    X, y, lo, hi, scr, smiles = load(blocks)
    fold, _ = butina_folds(smiles, seed=seed)
    direct = run_fold(X, y, scr, fold, 1, lam, seed, device, mode=mode)

    keep = fold != 1
    Xa = np.vstack([X[keep], X[~keep]])
    ya = np.vstack([y[keep], np.full((int((~keep).sum()), 4), np.nan, np.float32)])
    sa = np.vstack([scr[keep], np.full((int((~keep).sum()), 4), np.nan, np.float32)])
    fa = np.concatenate([np.zeros(int(keep.sum()), int), np.ones(int((~keep).sum()), int)])
    viaappend = run_fold(Xa, ya, sa, fa, 1, lam, seed, device, mode=mode)

    d = float(np.max(np.abs(direct - viaappend)))
    print(f"фолд 1, сид {seed}, lam {lam}: максимум |разности| между обычным путём и "
          f"путём с дописанными строками = {d:.3e}")
    print("совпало побитово" if d == 0.0 else
          "НЕ совпало --- дописанные строки куда-то дотягиваются, путь на тест использовать нельзя")
    return d


def evaluate(y, lo, hi, pred):
    """ST-RAE and Spearman per enzyme plus macro, over observed cells."""
    out = {}
    for e, c in enumerate(CYPS):
        m = ~np.isnan(y[:, e])
        out[c] = round(float(strae(y[m, e], pred[m, e],
                                   y_true_upper=hi[m, e], y_true_lower=lo[m, e])), 4)
        out["rho_" + c] = round(float(spearmanr(pred[m, e], y[m, e]).statistic), 3)
    out["MACRO"] = round(float(np.mean([out[c] for c in CYPS])), 4)
    out["MACRO_rho"] = round(float(np.mean([out["rho_" + c] for c in CYPS])), 3)
    return out


def main():
    global HIDDEN, DEPTH, DROPOUT, EPOCHS, WEIGHT_DECAY
    ap = argparse.ArgumentParser()
    ap.add_argument("--lams", default="0,1.0", help="lambda_scr values, comma separated")
    ap.add_argument("--seeds", default="0", help="split seeds, comma separated")
    ap.add_argument("--device", default="mps" if torch.backends.mps.is_available() else "cpu")
    ap.add_argument("--out", default=None, help="по умолчанию preds/trunk_<mode>.json")
    ap.add_argument("--hidden", type=int, default=HIDDEN)
    ap.add_argument("--depth", type=int, default=DEPTH)
    ap.add_argument("--dropout", type=float, default=DROPOUT)
    ap.add_argument("--epochs", type=int, default=EPOCHS)
    ap.add_argument("--wd", type=float, default=WEIGHT_DECAY)
    ap.add_argument("--blocks", default=BLOCKS)
    ap.add_argument("--swap-cal", default="",
                    help="обменять (E, h) между двумя ферментами, например CYP2D6,CYP3A4 - "
                         "единственное вмешательство, которое двигает карту прибора, "
                         "оставляя данные, метки, полосы, разбиение и информативность теми же")
    ap.add_argument("--noise", type=float, default=0.0,
                    help="eta: порча скринингового канала в его же ско, ранг падает в "
                         "sqrt(1+eta^2) раз, масштаб сохраняется")
    ap.add_argument("--mode", default="twohead", choices=["twohead", "calibrated"],
                    help="twohead: free second head. calibrated: screen via g(pi), one latent")
    ap.add_argument("--check-test-path", action="store_true",
                    help="проверить, что путь на тест --- это тот же run_fold, и выйти")
    a = ap.parse_args()

    if a.check_test_path:
        HIDDEN, DEPTH, DROPOUT = a.hidden, a.depth, a.dropout
        EPOCHS, WEIGHT_DECAY = a.epochs, a.wd
        for seed in [int(x) for x in a.seeds.split(",")]:
            check_test_path(seed=seed, lam=float(a.lams.split(",")[-1]),
                            mode=a.mode, blocks=a.blocks, device=a.device)
        return

    HIDDEN, DEPTH, DROPOUT = a.hidden, a.depth, a.dropout
    EPOCHS, WEIGHT_DECAY = a.epochs, a.wd
    # Swapping the instrument map between two enzymes. lambda = 0 does not touch g_of_pi at
    # all, so that arm is unaffected by construction and serves as the leak check.
    global CAL_E, CAL_H
    swap = ""
    if a.swap_cal:
        u, v = [CYPS.index(x.strip()) for x in a.swap_cal.split(",")]
        CAL_E = CAL_E.copy(); CAL_H = CAL_H.copy()
        CAL_E[[u, v]] = CAL_E[[v, u]]
        CAL_H[[u, v]] = CAL_H[[v, u]]
        swap = f"_swap{CYPS[u][3:]}-{CYPS[v][3:]}"
        print(f"карта прибора обменяна: {CYPS[u]} <-> {CYPS[v]}; "
              f"E {CAL_E.tolist()}, h {CAL_H.tolist()}", flush=True)

    if a.out is None:
        tag = ("" if a.noise == 0 else f"_noise{a.noise:g}") + swap
        a.out = RES + f"preds/trunk_{a.mode}{tag}.json"
    lams = [float(v) for v in a.lams.split(",")]
    seeds = [int(v) for v in a.seeds.split(",")]
    X, y, lo, hi, scr, smiles = load(a.blocks)
    print(f"X {X.shape} ({a.blocks}) | устройство {a.device}", flush=True)

    saved, table = {}, []
    for seed in seeds:
        fold, n_cl = butina_folds(smiles, seed=seed)
        # One draw per split seed, reused at every eta and in both modes: the ladder is
        # nested rather than independent, which takes the noise draw out of the comparison.
        zn = np.random.default_rng(90000 + seed).standard_normal(scr.shape).astype(np.float32)
        for lam in lams:
            t0 = time.time()
            pred = np.full_like(y, np.nan)
            for f in range(5):
                te = fold == f
                if te.sum() == 0:
                    continue
                pred[te] = run_fold(X, y, scr, fold, f, lam, seed, a.device, a.mode,
                                    zn=zn, eta=a.noise)
            r = evaluate(y, lo, hi, pred)
            r["seed"], r["lambda"], r["mode"], r["noise"] = seed, lam, a.mode, a.noise
            table.append(r)
            saved[f"{a.mode}|{seed}|{lam}"] = np.where(np.isnan(y), np.nan, pred).tolist()
            print(f"  [{a.mode}] сид {seed} lambda {lam:<4} макро {r['MACRO']:.4f} "
                  f"rho {r['MACRO_rho']:.3f}  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)[["seed", "lambda", *CYPS, "MACRO",
                              *[f"rho_{c}" for c in CYPS], "MACRO_rho"]]
    print()
    print(df.to_string(index=False))
    meta = {"mode": a.mode, "blocks": a.blocks, "hidden": a.hidden, "depth": a.depth,
            "dropout": a.dropout, "epochs": a.epochs, "wd": a.wd, "lr": LR, "batch": BATCH,
            "seeds": seeds, "lams": lams, "noise": a.noise, "swap_cal": a.swap_cal,
            "device": a.device,
            "pc0": PC0, "cal_e": CAL_E.tolist(), "cal_h": CAL_H.tolist(),
            "torch": torch.__version__, "numpy": np.__version__}
    json.dump({"table": table, "preds": saved, "meta": meta}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")


if __name__ == "__main__":
    main()
