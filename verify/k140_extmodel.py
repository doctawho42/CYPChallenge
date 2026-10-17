"""Their data, our learner: the arm that separates "the published model is bad" from "this source
does not transfer".

WHY THIS EXISTS. `verify/k139_oadmetcyp.py` measured the organisers' own CYP-trained model as a
member and every arm failed, the best frozen embedding landing at 0.4372 of standalone rank against
a random-weight control of 0.4833. Two explanations survive that measurement and they imply
opposite next moves: either their particular model generalises badly off its training distribution,
or the 8068 ChEMBL rows themselves carry no ordering that transfers to our unselected library.

The arm that separates them is to hold the data fixed and swap the learner. `verify/k103_extfeats.py`
already builds those 8068 rows in OUR feature space, aligned to `trunk.BLOCKS` by column NAME rather
than by count, with our 64 overlapping molecules dropped and the blinded-test overlap asserted at
zero (8004 rows survive, which is item 319's count). This file fits the project's own learners ---
`гребневая` as RidgeCV and `HistGB` --- on those rows alone and scores their predictions on our
molecules in the same device: standalone, as a sixth member, and as a seventh beside
`surrogate_adme`.

No cross-validation is used on our side and none is needed: the external rows are disjoint from
ours, so a prediction on our rows is out of sample however our folds are drawn. The fold enters
only through `score`'s per-fold affine pair.

THE THREE ROUTES, SAME 8068 ROWS. Item 319 put them through the trunk's third head and got a channel
of +0.0043 against a +0.0105 gate. Item 363 put them through the model their own authors trained and
got 0.4372 standalone, below random weights. This file is the third and last cheap route.

THE CONTROL THAT DECIDES WHETHER A POSITIVE RESULT MEANS ANYTHING. The same learner is fitted again
on the same rows with the external labels PERMUTED within enzyme --- same rows, same feature matrix,
same regulariser, no molecule-to-measurement link. If the shuffled arm scores like the real one, the
arm is measuring our own feature block rather than their labels, which is item 349's defect and the
reason `src/ablncgc.py` carries the same control.

    uv run python verify/k140_extmodel.py --seeds 0,1

Needs data/ext_thirdhead.npz (`uv run python verify/k103_extfeats.py`, about two minutes of RDKit).
Reads the member caches and item 349's probe caches. Writes results/logs/k140_extmodel.json.
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
from sklearn.linear_model import RidgeCV

from cypsplit import butina_folds, fold_digest
from submit import CYPS, _keep
import trunk as T

import k106_probe as K
from k106_probe import score

LOG = RES + "logs/k140_extmodel.json"
EXT = D + "ext_thirdhead.npz"
GOLDEN = "2d93c19815e14261"
MEMBERS = ("поферментно", "пул", "GP", "гребневая", "ствол")
FLOOR, FLOOR_PAIR = 0.0036, 0.007
SCREEN = 0.60
ITEM319 = 0.0043          # the same rows through a third head on the trunk
ITEM139 = {"зонд(эмб 2048)": 0.4372, "прямые (0 параметров)": 0.2881}   # item 363's arms
ALPHAS = np.array([1.0, 10.0, 100.0, 1000.0, 10000.0])


def fit_predict(Xe, ye, Xs, kind, seed=0):
    """One model per enzyme, trained ONLY on the external rows, predicting ours.

    No cross-validation is needed on our side and none is used: the external rows are disjoint
    from ours (k103 drops the 64 overlaps and asserts zero against the blinded test), so a
    prediction on our rows is out of sample however our folds are drawn. The fold enters only
    through score()'s affine pair.
    """
    out = []
    for e in range(len(CYPS)):
        m = np.isfinite(ye[:, e])
        Z, yy = Xe[m], ye[m, e]
        if kind == "гребневая":
            mu, sd = Z.mean(0), Z.std(0)
            keep = sd > 1e-6
            f = lambda B: (B[:, keep] - mu[keep]) / sd[keep]
            mdl = RidgeCV(alphas=ALPHAS).fit(f(Z), yy)
            out.append([mdl.predict(f(X)) for X in Xs])
        elif kind == "HistGB":
            mdl = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.06,
                                                min_samples_leaf=20, l2_regularization=1.0,
                                                random_state=seed).fit(Z, yy)
            out.append([mdl.predict(X) for X in Xs])
        else:
            raise ValueError(kind)
    return [[out[e][i] for e in range(len(CYPS))] for i in range(len(Xs))]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    if not _pl.Path(EXT).exists():
        print("нет %s --- сначала uv run python verify/k103_extfeats.py" % EXT, file=_sys.stderr)
        return 3
    z = np.load(EXT)
    Xe, ye = z["Xe"].astype(np.float64), z["ye"].astype(np.float64)
    X, Y, LO, HI, _scr, _smi = T.load(T.BLOCKS)
    X = np.asarray(X, np.float64)
    M = np.isfinite(Y)
    K.S.LO, K.S.HI = LO, HI
    rows = pd.read_csv(D + "rows.csv")
    print("внешних %s, наших %s; внешних меток: %s"
          % (Xe.shape, X.shape,
             ", ".join("%s %d" % (c, int(np.isfinite(ye[:, e]).sum()))
                       for e, c in enumerate(CYPS))))
    if Xe.shape[1] != X.shape[1]:
        print("КОНТРОЛЬ: ширина %d против %d" % (Xe.shape[1], X.shape[1]), file=_sys.stderr)
        return 2

    # ---- the arms: our own learner on their rows, and the same with their labels shuffled ---
    arms = {}
    for kind in ("гребневая", "HistGB"):
        t0 = time.time()
        arms[kind] = fit_predict(Xe, ye, [X], kind)[0]
        print("  %-12s обучена на внешних строках, %4.0f с" % (kind, time.time() - t0),
              flush=True)
    rng = np.random.default_rng(0)
    ysh = ye.copy()
    for e in range(len(CYPS)):
        m = np.isfinite(ysh[:, e])
        v = ysh[m, e].copy()
        rng.shuffle(v)
        ysh[m, e] = v
    arms["гребневая, метки перемешаны"] = fit_predict(Xe, ysh, [X], "гребневая")[0]
    print("  контроль: тот же учитель на ПЕРЕМЕШАННЫХ внешних метках", flush=True)

    per_seed, used = {}, []
    for s in seeds:
        mem_p = _pl.Path(RES + f"preds/members_seed{s}.json")
        if not mem_p.exists():
            print("  сид %d: нет %s --- пропуск" % (s, mem_p.name))
            continue
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        if s == 0 and dg != GOLDEN:
            print("ОТКАЗ: разбиение сдвинулось, сид 0 дал %s" % dg, file=_sys.stderr)
            return 2
        mem = json.load(open(mem_p))
        if set(mem) != set(MEMBERS):
            print("ОТКАЗ: кэш членов сида %d содержит %s" % (s, sorted(mem)), file=_sys.stderr)
            return 2
        pdz = [(k, [np.asarray(mem[k][e], float) for e in range(len(CYPS))]) for k in MEMBERS]
        kept = [[Q[e] for k, Q in pdz if _keep(e, k)] for e in range(len(CYPS))]
        cell = {"digest": dg, "пять": score([np.mean(kept[e], 0) for e in range(len(CYPS))],
                                            Y, LO, HI, M, fold)}
        sixth = None
        p_sa = _pl.Path(RES + f"preds/k122_surrogate_adme_probe_seed{s}.json")
        if p_sa.exists() and json.load(open(p_sa)).get("digest") == dg:
            sixth = [np.asarray(v, float) for v in json.load(open(p_sa))["P"]]
            cell["пять+sa"] = score([np.mean(kept[e] + [sixth[e]], 0) for e in range(len(CYPS))],
                                    Y, LO, HI, M, fold)
        print("\n  сид %d: пять = %.4f%s" % (s, cell["пять"][0],
              ", пять+sa = %.4f" % cell["пять+sa"][0] if sixth is not None else ""), flush=True)
        for tag, Pfull in arms.items():
            P = [Pfull[e][M[:, e]] for e in range(len(CYPS))]
            rec = {"solo": score(P, Y, LO, HI, M, fold),
                   "sixth": score([np.mean(kept[e] + [P[e]], 0) for e in range(len(CYPS))],
                                  Y, LO, HI, M, fold)}
            if sixth is not None:
                rec["seventh"] = score(
                    [np.mean(kept[e] + [sixth[e], P[e]], 0) for e in range(len(CYPS))],
                    Y, LO, HI, M, fold)
            cell[tag] = rec
            print("    %-28s одиночн. %.4f | шестым %+.4f | седьмым %+.4f"
                  % (tag, rec["solo"][0], rec["sixth"][0] - cell["пять"][0],
                     (rec["seventh"][0] - cell["пять+sa"][0]) if sixth is not None
                     else float("nan")), flush=True)
        per_seed[str(s)] = cell
        used.append(s)

    if not used:
        print("нет кэшей членов", file=_sys.stderr)
        return 3
    if 0 in used and abs(per_seed["0"]["пять"][0] - 0.6434) > 1e-4:
        print("КОНТРОЛЬ ПРОВАЛЕН: пять на сиде 0 дали %.4f против 0.6434"
              % per_seed["0"]["пять"][0], file=_sys.stderr)
        return 2
    print("\nконтроль: пять членов на сиде 0 дают %.4f (k132/k134/k135/k139: 0.6434)"
          % per_seed["0"]["пять"][0])

    print("\n%s" % ("=" * 92))
    print("%-28s %9s %9s | %9s %7s | %9s" % ("арм", "одиночн.", "Δпара", "шестым", "знак",
                                             "седьмым"))
    print("=" * 92)
    out = {"seeds": used, "per_seed": per_seed, "summary": {}}
    for tag in arms:
        sol = np.array([per_seed[str(s)][tag]["solo"][0] for s in used])
        solp = np.array([per_seed[str(s)][tag]["solo"][1] for s in used])
        d6 = np.array([per_seed[str(s)][tag]["sixth"][0] - per_seed[str(s)]["пять"][0]
                       for s in used])
        has7 = all("seventh" in per_seed[str(s)][tag] for s in used)
        d7 = (np.array([per_seed[str(s)][tag]["seventh"][0] - per_seed[str(s)]["пять+sa"][0]
                        for s in used]) if has7 else np.array([np.nan]))
        out["summary"][tag] = {"solo_rank": float(sol.mean()), "solo_pair": float(solp.mean()),
                               "d6": float(d6.mean()),
                               "sign6": [int((d6 > 0).sum()), len(d6)],
                               "d7": float(np.nanmean(d7))}
        print("%-28s %9.4f %9.4f | %+9.4f %5d/%d | %+9.4f"
              % (tag, sol.mean(), solp.mean(), d6.mean(), int((d6 > 0).sum()), len(d6),
                 np.nanmean(d7)))
    print("\nэталоны на ТЕХ ЖЕ 8068 строках:")
    print("  как СТРОКИ через третью голову ствола (пункт 319):      %+.4f канала" % ITEM319)
    for k, v in ITEM139.items():
        print("  как ОБУЧЕННАЯ ИМИ модель, %-22s одиночный %.4f" % (k, v))
    print("экран пункта 348: одиночный ранг %.2f; полы: ранг %.4f, ST-RAE %.3f"
          % (SCREEN, FLOOR, FLOOR_PAIR))
    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
    print("\nзапись: %s" % LOG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
