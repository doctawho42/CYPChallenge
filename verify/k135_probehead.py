"""The frozen-encoder probe has always had a LINEAR head. Item 324's control varied the encoder and
item 134's varied nothing else; this varies the head.

WHY. Item 324's device is "a RidgeCV head on frozen dimensions of somebody else's pretrained
encoder", and its discriminating control is the same RidgeCV on OUR OWN DESC+MECH block, which gives
-0.0151 at sign 0/4. That control answers "is it the encoder or the member shape" and answers it
well. It does not answer "is it the encoder or the HEAD": every one of the eighteen frozen
checkpoints in items 324, 335, 341, 342, 345, 347, 349 and 350 was read by the same linear map with
the same three alphas. A 2048-dimensional representation read linearly is a strong baseline and a
deliberate choice, but it has never been compared with anything.

And item 134's closure is the reason this is worth a run rather than a shrug. Neighbour-borrowing,
shared-function transfer, sample size and enzyme level were each refuted BY MEASUREMENT on OUR
feature matrix (items 110, 111, 125, 132) --- a learned representation is a different geometry, and
the refutations do not carry over to it by argument.

WHAT IS MEASURED. The same device in every respect but the head: per enzyme, our five Butina folds,
standardisation and the dead-column filter fitted on each fold's TRAINING rows only (item 202's
defect 2 not repeated), the prediction entering as a sixth member beside the five `_keep` admits,
scored by `verify/k106_probe.py:score` after the per-fold affine pair. Four seeds.

    гребневая (отгружается)   RidgeCV, alphas 1e2/1e3/1e4 --- пункт 324, и k134 подтвердил,
                              что 1e4 это внутренний оптимум, а не край
    ядровая RBF               KernelRidge, (alpha, gamma) выбираются ВНУТРЕННЕЙ CV по
                              обучающим строкам фолда --- ни один параметр не видит исход
    kNN косинусный            k выбирается той же внутренней CV; геометрия вложения, а не
                              наша матрица, где пункт 110 это закрыл
    HistGB                    пункт 335 мерил это на 600 столбцах как +0.0023; здесь 2048

NOTHING IS CHOSEN AGAINST THE OUTCOME. Every head's hyper-parameters come from a 3-fold inner CV on
the training rows of the outer fold, so the selection is nested exactly as RidgeCV's own generalised
cross-validation is. The picks are printed per enzyme, because a parameter pinned at the edge of its
own range would make the arm a statement about the range rather than about the head --- which is
what k134 had to rule out for the alphas.

THE CONTROLS. The linear arm must reproduce item 355's +0.0148 for this embedding, or the device is
not that device. Seed 0's Butina digest must be the golden `2d93c19815e14261`. The embedding's row
order is checked by name against `rows.csv` with the rolled-row control item 345 requires. And the
per-enzyme deltas are printed, since `_keep` admits five members on CYP2D6 and two or three
elsewhere, so a sixth does not enter the four equally.

    uv run python verify/k135_probehead.py --seeds 0,1,2,3

Writes results/logs/k135_probehead.json.
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
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.kernel_ridge import KernelRidge
from sklearn.linear_model import RidgeCV
from sklearn.model_selection import GridSearchCV, KFold
from sklearn.neighbors import KNeighborsRegressor

from cypsplit import butina_folds, fold_digest
from submit import CYPS, _keep

import k106_probe as K
from k106_probe import score

LOG = RES + "logs/k135_probehead.json"
GOLDEN = "2d93c19815e14261"
MEMBERS = ("поферментно", "пул", "GP", "гребневая", "ствол")
ALPHAS = [100.0, 1000.0, 10000.0]
ITEM355 = 0.0148
INNER = 3


def head_linear(Ztr, ytr, Zte):
    r = RidgeCV(alphas=np.asarray(ALPHAS)).fit(Ztr, ytr)
    return r.predict(Zte), {"alpha": float(r.alpha_)}


def head_rbf(Ztr, ytr, Zte):
    # gamma scaled by dimension, the usual parameterisation; the grid is wide and the pick is
    # reported so a pinned optimum is visible rather than silent.
    d = Ztr.shape[1]
    # Deliberately two values each: the full 3x4 grid is 36 kernel builds per fold and ran at
    # about an hour per seed, which buys resolution this arm does not need -- it is here to say
    # whether a nonlinear head beats a linear one at all, not to tune one.
    g = [0.3 / d, 1.0 / d]
    gs = GridSearchCV(KernelRidge(kernel="rbf"),
                      {"alpha": [1.0, 10.0], "gamma": g},
                      cv=KFold(INNER, shuffle=True, random_state=0),
                      scoring="neg_mean_absolute_error", n_jobs=1)
    gs.fit(Ztr, ytr)
    b = gs.best_params_
    return gs.predict(Zte), {"alpha": float(b["alpha"]), "gamma_x_d": float(b["gamma"] * d)}


def head_knn(Ztr, ytr, Zte):
    gs = GridSearchCV(KNeighborsRegressor(metric="cosine", weights="distance"),
                      {"n_neighbors": [5, 10, 25, 50, 100]},
                      cv=KFold(INNER, shuffle=True, random_state=0),
                      scoring="neg_mean_absolute_error", n_jobs=1)
    gs.fit(Ztr, ytr)
    return gs.predict(Zte), {"k": int(gs.best_params_["n_neighbors"])}


def head_gbm(Ztr, ytr, Zte):
    m = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.06, max_depth=None,
                                      min_samples_leaf=20, l2_regularization=1.0,
                                      random_state=0).fit(Ztr, ytr)
    return m.predict(Zte), {}


HEADS = {"гребневая (отгр.)": head_linear, "ядровая RBF": head_rbf,
         "kNN косинусный": head_knn, "HistGB": head_gbm}


def oof_head(fn, E, Y, M, fold, e):
    """probe_oof's loop exactly, with the head swapped."""
    m = M[:, e]
    Xm, ym, fm = E[m], Y[m, e], fold[m]
    p = np.full(int(m.sum()), np.nan)
    picks = []
    for f in range(5):
        trn, te = fm != f, fm == f
        if te.sum() == 0:
            continue
        mu, sd = Xm[trn].mean(0), Xm[trn].std(0)
        keep = sd > 1e-6
        Z = (Xm[:, keep] - mu[keep]) / sd[keep]
        pred, info = fn(Z[trn], ym[trn], Z[te])
        p[te] = pred
        picks.append(info)
    if not np.isfinite(p).all():
        raise SystemExit("%s: голова оставила незаполненные строки" % CYPS[e])
    return p, picks


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--emb", default=D + "emb_surrogate_adme_train.npz")
    ap.add_argument("--seeds", default="0,1,2,3")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    for f in ("rows.csv", "cyp-challenge-TRAIN_inhibition.csv"):
        if not _pl.Path(D + f).exists():
            print("нет data/%s --- make setup" % f, file=_sys.stderr)
            return 3
    if not _pl.Path(a.emb).exists():
        print("нет %s" % a.emb, file=_sys.stderr)
        return 3

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)
    K.S.LO, K.S.HI = LO, HI

    z = np.load(a.emb, allow_pickle=True)
    E = z["train"].astype(np.float64)
    nm = np.asarray(z["molecule_name"], object)
    r = rows.Molecule_Name.to_numpy(object)
    ok, rolled = bool((nm == r).all()), bool((np.roll(nm, 1) == r).all())
    print("вложение %s %s, порядок совпал: %s; сдвинутый контроль (ожидается False): %s"
          % (_pl.Path(a.emb).name, E.shape, ok, rolled), flush=True)
    if not ok or rolled:
        print("КОНТРОЛЬ ПРОВАЛЕН: порядок вложения", file=_sys.stderr)
        return 2

    out = {"emb": _pl.Path(a.emb).name, "seeds": seeds, "per_seed": {}}
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
        for name, fn in HEADS.items():
            t0 = time.time()
            P, picks = [], {}
            for e, c in enumerate(CYPS):
                pe, pk = oof_head(fn, E, Y, M, fold, e)
                P.append(pe)
                picks[c] = pk
            six = score([np.mean(kept[e] + [P[e]], 0) for e in range(len(CYPS))],
                        Y, LO, HI, M, fold)
            solo = score(P, Y, LO, HI, M, fold)
            cell["arms"][name] = {"six": six, "solo": solo, "picks": picks,
                                  "d_rank": six[0] - five[0], "d_pair": six[1] - five[1],
                                  "secs": round(time.time() - t0, 1)}
            print("    %-20s шестым %+.4f ранга, %+.4f пары | одиночн. %.4f | %.0f с"
                  % (name, six[0] - five[0], six[1] - five[1], solo[0], time.time() - t0),
                  flush=True)
            print("        выбор: " + "; ".join(
                "%s %s" % (c[3:], picks[c][0] if picks[c] and picks[c][0] else "—")
                for c in CYPS), flush=True)
        out["per_seed"][str(s)] = cell
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)

    have = [s for s in seeds if str(s) in out["per_seed"]]
    if not have:
        print("ни одного сида", file=_sys.stderr)
        return 3
    print("\n" + "=" * 92)
    print("%-20s %9s %8s %7s %9s %9s" % ("голова", "Δранг", "sd", "знак", "Δпара", "одиночн."))
    print("=" * 92)
    out["summary"] = {}
    for name in HEADS:
        v = np.asarray([out["per_seed"][str(s)]["arms"][name]["d_rank"] for s in have], float)
        pr = np.asarray([out["per_seed"][str(s)]["arms"][name]["d_pair"] for s in have], float)
        so = float(np.mean([out["per_seed"][str(s)]["arms"][name]["solo"][0] for s in have]))
        sd = float(np.std(v, ddof=1)) if len(v) > 1 else float("nan")
        out["summary"][name] = {"d_rank": float(v.mean()), "sd": sd, "d_pair": float(pr.mean()),
                                "solo": so, "sign": [int((v > 0).sum()), len(v)]}
        print("%-20s %+9.4f %8.4f %5d/%d %+9.4f %9.4f"
              % (name, v.mean(), sd, int((v > 0).sum()), len(v), pr.mean(), so))
    lin = out["summary"]["гребневая (отгр.)"]
    print("\nконтроль: линейная голова даёт %+.4f против %+.4f пункта 355 (|Δ| %.4f)"
          % (lin["d_rank"], ITEM355, abs(lin["d_rank"] - ITEM355)))
    if abs(lin["d_rank"] - ITEM355) > 0.003:
        print("КОНТРОЛЬ ПРОВАЛЕН: устройство не воспроизводит пункт 355", file=_sys.stderr)
        return 2
    best = max(out["summary"], key=lambda k: out["summary"][k]["d_rank"])
    print("лучшая голова: %s, %+.4f --- против линейной %+.4f (пол 0.0036)"
          % (best, out["summary"][best]["d_rank"],
             out["summary"][best]["d_rank"] - lin["d_rank"]))
    json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
    print("\nзапись: %s" % LOG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
