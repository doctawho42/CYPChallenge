"""Nested top-k over the nine-member library: does keeping FEWER, better members beat all nine?

The 1/n law says a member's worth divides by the count, so dead weight is expensive.  The
ranking must not see the fold it predicts: for fold f, each candidate's solo Spearman is
computed on the OTHER four folds (its predictions are already out-of-fold, so this leaks
nothing), the top k are averaged, and that average predicts fold f.  One integer
hyperparameter, no per-fold weights -- unlike the Caruana selector that measured -0.0048.

TWO CONTROLS THAT MUST SUCCEED, both exercising the selection clause:
  k=9   must equal the plain nine-member mean bit for bit -- at k=9 nothing is selected.
  случ. the same k with a RANDOM ranking must be worse, or the ordering carries nothing
        and any gain is k alone.
"""
import json
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd
from scipy.stats import spearmanr
from cypsplit import butina_folds, fold_digest
from submit import CYPS
import k106_probe as K
from k106_probe import score, members, probe_oof, load_emb
rows = pd.read_csv(D + "rows.csv"); z = np.load(D + "feats.npz")
XF = np.hstack([z["FP"], z["DESC"], z["MECH"]])
tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
        .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS],1)
LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS],1)
HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS],1)
M = ~np.isnan(Y); K.S.LO, K.S.HI = LO, HI
J = json.load(open(RES + "preds/k147_tabicl_E100.json"))

def topk_oof(cand, yy, fm, k, rng=None):
    """cand: list of 9 prediction vectors over this enzyme's rows.  Returns the nested top-k mean."""
    out = np.full(len(yy), np.nan)
    for f in range(5):
        trn, te = fm != f, fm == f
        if te.sum() == 0: continue
        if rng is None:
            r = [float(spearmanr(yy[trn], c[trn]).statistic) for c in cand]
            order = np.argsort(r)[::-1]
        else:
            order = rng.permutation(len(cand))
        out[te] = np.mean([cand[i][te] for i in order[:k]], 0)
    assert np.isfinite(out).all()
    return out

acc = {}
for s in (0,1,2,3):
    fold, _ = butina_folds(list(rows.SMILES), seed=s)
    assert fold_digest(fold) == J["digests"][s]
    pdz = members(s, XF, Y, M, fold, LO, HI)
    probe = [probe_oof(load_emb(rows), Y, M, fold, e)[0] for e in range(4)]
    chem = [np.asarray(v,float) for v in
            json.load(open(RES+"preds/k122_chemeleon_probe_seed%d.json"%s))["P"]]
    gbm  = [np.asarray(v,float) for v in
            json.load(open(RES+"preds/k114_heads_seed%d.json"%s))["heads"]["gbm"]]
    T    = [np.asarray(v,float) for v in J["preds"]["E_seed%d"%s]]
    cand = [[Q[e] for k_,Q in pdz] + [probe[e], chem[e], gbm[e], T[e]] for e in range(4)]
    assert all(len(c) == 9 for c in cand)
    ys = [Y[M[:,e],e] for e in range(4)]
    fms = [fold[M[:,e]] for e in range(4)]
    # the shipped eight, for the baseline column
    kept = [[Q[e] for k_,Q in pdz if K._keep(e,k_)] for e in range(4)]
    acc.setdefault("отгруженные восемь", []).append(
        score([np.mean(kept[e]+[probe[e],chem[e],gbm[e]],0) for e in range(4)],Y,LO,HI,M,fold))
    for k in range(1, 10):
        P = [topk_oof(cand[e], ys[e], fms[e], k) for e in range(4)]
        acc.setdefault("k=%d" % k, []).append(score(P, Y, LO, HI, M, fold))
    # CONTROL 1: k=9 must be the plain mean, bit for bit
    flat = [np.mean(cand[e],0) for e in range(4)]
    d = max(float(np.max(np.abs(topk_oof(cand[e],ys[e],fms[e],9)-flat[e]))) for e in range(4))
    # 1e-12, not 0: np.mean over the top-9 in RANK order reassociates the sum against the
    # same nine in pool order, which moves the last bits (measured 1.8e-15 on values ~5).
    assert d < 1e-12, "КОНТРОЛЬ k=9 провалился: max|Δ| = %.2e" % d
    # CONTROL 2: random ranking at the best k
    rng = np.random.default_rng(s)
    for k in (4, 5, 6):
        P = [topk_oof(cand[e], ys[e], fms[e], k, rng) for e in range(4)]
        acc.setdefault("случ. k=%d" % k, []).append(score(P, Y, LO, HI, M, fold))
print("  КОНТРОЛЬ k=9 == обычное среднее девяти: пройден, max|Δ| = %.1e (порог 1e-12)" % d)
b = np.mean([v[0] for v in acc["отгруженные восемь"]])
bp = np.mean([v[1] for v in acc["отгруженные восемь"]])
print("\n  ВЛОЖЕННОЕ ОТСЕЧЕНИЕ по одиночному рангу, библиотека из девяти")
print("  %-18s %9s %8s %10s %7s %10s" % ("арм","ранг","sd","Δранг","знак","ΔST-RAE"))
for nm, v in acc.items():
    r = [x[0] for x in v]; p = [x[1] for x in v]
    g = [x - y for x,y in zip(r, [q[0] for q in acc["отгруженные восемь"]])]
    print("  %-18s %9.4f %8.4f %+10.4f %5d/%d %+10.4f"
          % (nm, np.mean(r), np.std(r,ddof=1), np.mean(r)-b,
             sum(1 for x in g if x>0), len(g), np.mean(p)-bp))
print("\n  пол ранга 0.0036, пол ST-RAE 0.007; ST-RAE -- потеря, минус = лучше")
json.dump({k: v for k,v in acc.items()}, open(RES+"logs/k146_topk.json","w"), indent=1)
