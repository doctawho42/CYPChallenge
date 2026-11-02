"""The metric's loss is the epsilon-insensitive loss, and this project has never run the learner
that has it natively. `grep -rn "SVR\\|SVC" src/` returns nothing.

THE ARGUMENT, and it is short. ST-RAE pays exactly zero for a prediction inside the compound's
published band and the distance to the nearer edge outside it (the organisers'
`rae_soft_threshold_absolute_error`). That is the epsilon-insensitive loss with a per-compound
epsilon. This project derived that and built it into boosting as the DEAD ZONE, which is its largest
single gain on record --- +0.033 of rank, item 148, "структура МЕТРИКИ в потере" in item 189's table.
Support vector regression minimises that loss by construction and appears nowhere in the tree.

It was not dismissed, it was never considered. What WAS dismissed, with a measurement, is a
different thing: `src/gp.py` rejects the standard Tanimoto-kernel GP over Morgan fingerprints
because item 90 measured Morgan last of four spaces for carrying activity between neighbours and
standardised RDKit descriptors first by a factor of 2.5, so the kernel went on descriptors. That is
a statement about the SIMILARITY, not about the loss, and an SVR on the same descriptors inherits
it rather than contradicting it.

And the board says a plain learner is enough to beat us on ordering. The entrant at rank 7 is called
`SVM` (rho 0.7710 against our 0.7062) and the one at 22 is called `random-forest`, whose per-enzyme
Spearman is 0.8463 on CYP2C9 and 0.8433 on CYP3A4 against our 0.8231 and 0.8146. Neither needs a
foundation model or a GPU.

WHAT IS MEASURED. An epsilon-SVR member on OUR feature matrix, in item 324's device exactly: per
enzyme, our five Butina folds, standardisation on each fold's TRAINING rows only, the prediction
entering as a sixth member beside the five `_keep` admits, scored after the per-fold affine pair.
Two epsilon arms, which is the whole experiment:

    eps = полоса     epsilon = медиана (hi - lo)/2 по ОБУЧАЮЩИМ строкам фолда --- потеря метрики
    eps = 0          то же самое без нечувствительной зоны --- контроль, который отделяет
                     «SVR помогает» от «потеря метрики помогает»

C and gamma come from a 3-fold inner CV on the training rows under mean absolute error, so no
hyper-parameter sees the outcome; the picks are printed, because an optimum at the edge of its grid
would make the arm a statement about the grid (which is what k134 had to rule out for the probe's
alphas).

THE CONTROLS. Seed 0's Butina digest must be the golden `2d93c19815e14261`. The five-member
reference must reproduce across seeds. The two epsilon arms share everything else, so their
difference is the metric's loss and nothing else. And the per-enzyme deltas are printed, since
`_keep` admits five members on CYP2D6 and two or three elsewhere.

    uv run python verify/k136_svr.py --seeds 0,1

Writes results/logs/k136_svr.json.
Exit: 0 ran, 2 a control failed, 3 an input is missing.
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
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.svm import SVR

from cypsplit import butina_folds, fold_digest
from submit import CYPS, _keep

import k106_probe as K
from k106_probe import score

LOG = RES + "logs/k136_svr.json"
GOLDEN = "2d93c19815e14261"
MEMBERS = ("поферментно", "пул", "GP", "гребневая", "ствол")
FLOOR = 0.0036


def svr_oof(X, Y, LO, HI, M, fold, e, eps_mode, grid):
    m = M[:, e]
    Xm, ym, fm = X[m], Y[m, e], fold[m]
    lo, hi = LO[m, e], HI[m, e]
    p = np.full(int(m.sum()), np.nan)
    picks = []
    for f in range(5):
        trn, te = fm != f, fm == f
        if te.sum() == 0:
            continue
        mu, sd = Xm[trn].mean(0), Xm[trn].std(0)
        keep = sd > 1e-6
        Z = (Xm[:, keep] - mu[keep]) / sd[keep]
        # epsilon from the TRAINING rows' own bands, never the held-out ones
        eps = float(np.median((hi[trn] - lo[trn]) / 2.0)) if eps_mode == "полоса" else 0.0
        gs = GridSearchCV(SVR(kernel="rbf", epsilon=eps), grid,
                          cv=KFold(3, shuffle=True, random_state=0),
                          scoring="neg_mean_absolute_error", n_jobs=1)
        gs.fit(Z[trn], ym[trn])
        p[te] = gs.predict(Z[te])
        b = dict(gs.best_params_)
        b["eps"] = round(eps, 4)
        picks.append(b)
    if not np.isfinite(p).all():
        raise SystemExit("%s: SVR оставил незаполненные строки" % CYPS[e])
    return p, picks


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1")
    ap.add_argument("--emb", default=None,
                    help="вместо нашей матрицы взять замороженное вложение из .npz")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    for f in ("rows.csv", "feats.npz", "cyp-challenge-TRAIN_inhibition.csv"):
        if not _pl.Path(D + f).exists():
            print("нет data/%s --- make setup" % f, file=_sys.stderr)
            return 3

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)
    K.S.LO, K.S.HI = LO, HI

    if a.emb:
        z = np.load(a.emb, allow_pickle=True)
        X = z["train"].astype(np.float64)
        nm = np.asarray(z["molecule_name"], object)
        r = rows.Molecule_Name.to_numpy(object)
        if not bool((nm == r).all()) or bool((np.roll(nm, 1) == r).all()):
            print("КОНТРОЛЬ ПРОВАЛЕН: порядок вложения", file=_sys.stderr)
            return 2
        what = _pl.Path(a.emb).name
    else:
        z = np.load(D + "feats.npz")
        X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float64)
        what = "FP+DESC+MECH"
    print("матрица: %s %s" % (what, X.shape), flush=True)
    d = X.shape[1]
    grid = {"C": [1.0, 10.0], "gamma": ["scale", 0.1 / d]}

    out = {"matrix": what, "shape": list(X.shape), "seeds": seeds, "grid": {k: [str(x) for x in v]
           for k, v in grid.items()}, "per_seed": {}}
    for s in seeds:
        fold, _n = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        if s == 0 and dg != GOLDEN:
            print("ОТКАЗ: сид 0 дал %s" % dg, file=_sys.stderr)
            return 2
        mp = _pl.Path(RES + f"preds/members_seed{s}.json")
        if not mp.exists():
            print("  сид %d: нет %s --- пропуск" % (s, mp.name))
            continue
        mem = json.load(open(mp))
        pdz = [(k, [np.asarray(mem[k][e], float) for e in range(len(CYPS))]) for k in MEMBERS]
        kept = [[Q[e] for k, Q in pdz if _keep(e, k)] for e in range(len(CYPS))]
        five = score([np.mean(kept[e], 0) for e in range(len(CYPS))], Y, LO, HI, M, fold)
        cell = {"digest": dg, "пять": five, "arms": {}}
        print("\n  сид %d: пять = ранг %.4f, пара %.4f" % (s, five[0], five[1]), flush=True)
        for mode in ("полоса", "0"):
            t0 = time.time()
            P, picks = [], {}
            for e, c in enumerate(CYPS):
                pe, pk = svr_oof(X, Y, LO, HI, M, fold, e, mode, grid)
                P.append(pe)
                picks[c] = pk
                print("      %s %-7s %s  (%.0f с)"
                      % (c, mode, pk[0] if pk else "", time.time() - t0), flush=True)
            six = score([np.mean(kept[e] + [P[e]], 0) for e in range(len(CYPS))],
                        Y, LO, HI, M, fold)
            solo = score(P, Y, LO, HI, M, fold)
            cell["arms"]["eps=" + mode] = {
                "six": six, "solo": solo, "picks": picks,
                "d_rank": six[0] - five[0], "d_pair": six[1] - five[1],
                "secs": round(time.time() - t0, 1)}
            print("    eps=%-7s шестым %+.4f ранга, %+.4f пары | одиночн. %.4f | %.0f с"
                  % (mode, six[0] - five[0], six[1] - five[1], solo[0], time.time() - t0),
                  flush=True)
        out["per_seed"][str(s)] = cell
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)

    have = [s for s in seeds if str(s) in out["per_seed"]]
    if not have:
        print("ни одного сида", file=_sys.stderr)
        return 3
    print("\n" + "=" * 86)
    print("%-14s %9s %8s %7s %9s %9s" % ("арм", "Δранг", "sd", "знак", "Δпара", "одиночн."))
    print("=" * 86)
    out["summary"] = {}
    for name in ("eps=полоса", "eps=0"):
        v = np.asarray([out["per_seed"][str(s)]["arms"][name]["d_rank"] for s in have], float)
        pr = np.asarray([out["per_seed"][str(s)]["arms"][name]["d_pair"] for s in have], float)
        so = float(np.mean([out["per_seed"][str(s)]["arms"][name]["solo"][0] for s in have]))
        sd = float(np.std(v, ddof=1)) if len(v) > 1 else float("nan")
        out["summary"][name] = {"d_rank": float(v.mean()), "sd": sd, "d_pair": float(pr.mean()),
                                "solo": so, "sign": [int((v > 0).sum()), len(v)]}
        print("%-14s %+9.4f %8.4f %5d/%d %+9.4f %9.4f"
              % (name, v.mean(), sd, int((v > 0).sum()), len(v), pr.mean(), so))
    g = out["summary"]["eps=полоса"]["d_rank"] - out["summary"]["eps=0"]["d_rank"]
    print("\nпотеря метрики против обычной SVR: %+.4f ранга (пол %.4f)" % (g, FLOOR))
    print("по ферментам, Δранг против пяти:")
    for name in ("eps=полоса", "eps=0"):
        dd = np.mean([[out["per_seed"][str(s)]["arms"][name]["six"][2][e]
                       - out["per_seed"][str(s)]["пять"][2][e] for e in range(len(CYPS))]
                      for s in have], 0)
        out["summary"][name]["d_rank_per"] = [float(x) for x in dd]
        print("  %-14s " % name + "  ".join("%s %+.4f" % (c[3:], dd[e])
                                            for e, c in enumerate(CYPS)))
    json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
    print("\nзапись: %s" % LOG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
