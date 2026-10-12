"""Can a frozen third-party encoder order CYP2D6 TDI, where our own features cannot?

THE QUESTION, and why it is this one. On the classification board we sit 80th of 144, and the
weakness is entirely in one cell: CYP2D6 at 100th against a field whose 90th percentile is 0.3655
MCC. Our CYP2D6 TDI classifier has out-of-fold **AUC 0.5855** -- the value `src/submit.py`'s own
comment states and this file re-measures, against a control that shuffled labels give 0.4971. MCC
is bounded by ordering, so no threshold rescues an AUC of 0.59: the best MCC over the whole
threshold grid is 0.1272 and the board gives 0.1255, which agree. The threshold is NOT the lever,
and the first version of this investigation thought it was.

WHY THE FROZEN ENCODER AND NOT SOMETHING ELSE. Item 324 measured a ridge head on the frozen
chemprop_medium embedding adding +0.0138 of macro rank to the regression ensemble, sign 8/8, with
a DISCRIMINATING control -- the same head with the same alpha grid on our own DESC+MECH block gives
-0.0151, sign 0/4 -- so the frozen pretrained weights buy it and the architecture buys nothing.
Those embeddings are already on disk. And the axis is untouched for TDI: a grep of
`verify/README.md` for TDI together with embeddings/chemprop/TabICL/TabPFN returns nothing, against
a control that `chemprop` alone appears 11 times. jeremy sits 17th on this board at 0.4025 using
TabICL on the same frozen embeddings, which is why the axis is worth a measurement rather than an
argument.

ARMS, with the control that makes a positive result mean something:

    A  наши признаки          HistGradientBoostingClassifier на FP+DESC+MECH --- то, что отгружено
    B  замороженный энкодер   логистическая на 600 замороженных измерениях chemprop
    C  среднее A и B          объединение вероятностей
    D  РАЗЛИЧАЮЩИЙ КОНТРОЛЬ   та же логистическая на НАШЕМ DESC+MECH

D is the arm that decides what a win by B means. If B beats A and D does not, the frozen weights
are doing the work, as item 324 found for regression. If D also beats A, the gain is the linear
head and has nothing to do with the encoder.

Reported per enzyme: AUC and the best MCC over the threshold grid, out of fold on the Butina split,
over several seeds. **AUC is the primary quantity here**, not MCC: MCC mixes ordering with the
operating point, and the operating point is already handled by the shipped Platt-plus-plug-in rule
(items 233, 235). A gain that shows in MCC but not in AUC is a threshold artefact.

    uv run python verify/k113_tdiemb.py --seeds 0,1,2,3

Writes results/logs/k113_tdiemb.json. Exit: 0 ran, 2 inputs missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES

import argparse
import json

import numpy as np
import pandas as pd
from sklearn.ensemble import HistGradientBoostingClassifier
from sklearn.linear_model import LogisticRegressionCV
from sklearn.metrics import matthews_corrcoef, roc_auc_score

from cypsplit import butina_folds, fold_digest

TDI = ["CYP3A4", "CYP2D6"]
GRID = np.arange(0.05, 0.96, 0.01)
CS = [0.003, 0.01, 0.03, 0.1, 0.3, 1.0]


def best_mcc(y, p):
    return max(matthews_corrcoef(y, (p >= t).astype(int)) for t in GRID)


def oof_gbm(X, y, fold):
    p = np.zeros(len(y))
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0 or y[trn].sum() == 0:
            continue
        p[te] = HistGradientBoostingClassifier(
            max_iter=300, learning_rate=0.06, random_state=0).fit(X[trn], y[trn]
            ).predict_proba(X[te])[:, 1]
    return p


def oof_logit(X, y, fold):
    """Logistic head, out of fold, with C chosen by an INNER cross-validation on the training half.

    The first version picked C by AUC on the training rows themselves, which is in-sample and
    always selects the loosest penalty -- a selection that cannot fail is not a selection.
    `LogisticRegressionCV` does the inner split properly. Standardisation and the dead-column
    filter are fitted on the training rows of each fold only, never on the union, which is item
    202's defect 2.

    This head is NOT where the cost of this file sits, and a PCA-100 version of it was written and
    then withdrawn for exactly that reason. One timed fold of CYP3A4: arm A 366 s, arm B 11 s, arm
    D 8 s. The HistGradientBoosting baseline on 2295 columns is 95 per cent of a seed; cheapening
    the logistic head buys a twentieth of nothing. The arm therefore keeps the 600 frozen
    dimensions item 324 validated, rather than a basis chosen to save time that was never spent
    here. The same timing run also reads load average 255 on eight cores with 7.3 GB of 8 swapped,
    so the 366 s is an upper bound on a thrashing box and not a cost of arm A either.
    """
    p = np.zeros(len(y))
    for f in range(5):
        trn, te = fold != f, fold == f
        if te.sum() == 0 or y[trn].sum() == 0:
            continue
        mu, sd = X[trn].mean(0), X[trn].std(0)
        keep = sd > 1e-6
        Z = (X[:, keep] - mu[keep]) / sd[keep]
        m = LogisticRegressionCV(Cs=CS, cv=4, scoring="roc_auc", max_iter=2000,
                                 n_jobs=-1).fit(Z[trn], y[trn])
        p[te] = m.predict_proba(Z[te])[:, 1]
    return p


def _probs_cached(cyp, seed, digest, m, y, fm, XFULL, EMB, XOWN):
    """Arms A, B and D out of fold, cached on disk under the fold digest.

    Predictions are valid only for the fold vector that produced them, so a digest mismatch
    recomputes rather than silently returning vectors from another split.
    """
    path = RES + "preds/k113_probs_seed%d.json" % seed
    store = json.load(open(path)) if _pl.Path(path).exists() else {}
    if store.get("digest") != digest:
        store = {"digest": digest}
    cell = store.get(cyp)
    if cell and all(len(cell[k]) == len(y) for k in ("A", "B", "D")):
        print("    вероятности %s из кэша" % cyp, flush=True)
        return (np.asarray(cell["A"], float), np.asarray(cell["B"], float),
                np.asarray(cell["D"], float))
    pa = oof_gbm(XFULL[m], y, fm)
    pb = oof_logit(EMB[m], y, fm)
    pd_ = oof_logit(XOWN[m], y, fm)
    store[cyp] = {"A": [float(x) for x in pa], "B": [float(x) for x in pb],
                  "D": [float(x) for x in pd_]}
    _pl.Path(RES + "preds").mkdir(parents=True, exist_ok=True)
    json.dump(store, open(path, "w"))
    return pa, pb, pd_


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--out", default=RES + "logs/k113_tdiemb.json")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    rows = pd.read_csv(D + "rows.csv")
    z = np.load(D + "feats.npz")
    XFULL = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    XOWN = np.hstack([z["DESC"], z["MECH"]])
    e = np.load(D + "emb_chemprop_medium.npz", allow_pickle=True)
    if "train" not in e:
        print("нет ключа train в emb_chemprop_medium.npz:", list(e.keys()))
        return 2
    EMB = e["train"]
    if len(EMB) != len(rows):
        print("вложение %d строк против rows.csv %d" % (len(EMB), len(rows)))
        return 2
    if "molecule_name" in e:
        same = list(map(str, e["molecule_name"])) == list(map(str, rows.Molecule_Name))
        print("порядок вложения совпал с rows.csv: %s  (сдвинутый контроль, ожидается False: %s)" % (
            same, list(map(str, e["molecule_name"]))[1:] == list(map(str, rows.Molecule_Name))[1:-1]))
        if not same:
            print("ПОРЯДОК НЕ СОВПАЛ --- дальше считать нельзя")
            return 2
    tdi = (pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv")
             .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())

    # Flush after every seed and merge with what is already logged. The four-seed run was killed
    # by a background time limit with nothing written, because the file was only saved at the end;
    # a measurement a kill can erase is a measurement that has to be repeated from zero.
    out = {}
    if _pl.Path(a.out).exists():
        try:
            out = json.load(open(a.out))
        except Exception:
            out = {}

    def flush():
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(a.out, "w"), ensure_ascii=False, indent=1)

    for c in TDI:
        col = f"{c}_is_TDI"
        yv = tdi[col].to_numpy(float)
        m = ~np.isnan(yv)
        y = yv[m].astype(int)
        print(f"\n=== {c}: n={m.sum()}, положительных {y.mean():.3f} ===")
        print("%-28s %8s %8s   %8s %8s" % ("арм", "AUC", "MCC", "ΔAUC к A", "ΔMCC к A"))
        per = {nm: [tuple(x) for x in v.get("пары", [])]
               for nm, v in out.get(c, {}).get("арм", {}).items()}
        done = list(out.get(c, {}).get("seeds", []))
        for s in seeds:
            if s in done:
                print("  сид %d уже в логе, пропуск" % s)
                continue
            fold, _ = butina_folds(list(rows.SMILES), seed=s)
            fm = fold[m]
            # The probabilities are cached, not just their summaries. The first version stored
            # AUC and the best MCC over a threshold grid and nothing else, which is an ORACLE
            # threshold -- so the arm could not be chosen on the quantity the board actually
            # pays for (MCC after the shipped Platt-plus-plug-in rule, items 233/235) without
            # paying for the whole sweep again. Same bill item 335 paid twice.
            pa, pb, pd_ = _probs_cached(c, s, fold_digest(fold), m, y, fm, XFULL, EMB, XOWN)
            pc = (pa + pb) / 2
            for nm, p in (("A наши признаки", pa), ("B замороженный энкодер", pb),
                          ("C среднее A и B", pc), ("D контроль DESC+MECH", pd_)):
                per.setdefault(nm, []).append((roc_auc_score(y, p), best_mcc(y, p)))
            done.append(s)
            out[c] = {"n": int(m.sum()), "pos": float(y.mean()), "seeds": done,
                      "арм": {nm: {"auc": float(np.mean([x[0] for x in v])),
                                   "mcc": float(np.mean([x[1] for x in v])),
                                   "пары": [[float(x[0]), float(x[1])] for x in v]}
                              for nm, v in per.items()}}
            flush()
            print("  сид %d посчитан и записан" % s)
        base = np.array(per["A наши признаки"])
        for nm in ("A наши признаки", "B замороженный энкодер", "C среднее A и B",
                   "D контроль DESC+MECH"):
            v = np.array(per[nm])
            print("%-28s %8.4f %8.4f   %+8.4f %+8.4f" % (
                nm, v[:, 0].mean(), v[:, 1].mean(),
                v[:, 0].mean() - base[:, 0].mean(), v[:, 1].mean() - base[:, 1].mean()))
        nn = len(per["A наши признаки"])
        sg = {nm: int(sum(1 for i in range(nn) if np.array(per[nm])[i, 0] > base[i, 0]))
              for nm in per}
        print("  знак по AUC (из %d сидов): %s" % (nn, sg))
        out[c]["знак_auc"] = sg
        flush()
    flush()
    print(f"\nсохранено: {a.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
