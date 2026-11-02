"""Do the two sub-floor TabICL members ADD?  Item 346's precedent says sub-floor members can.

Each alone is under the floor on the base of eight: E +0.0014, S +0.0028.  Item 346 added three
such members and the sum cleared both floors, so the sum is the measurement, not either part.

THE BASELINE IS ITEM 346's EIGHT, which is what `--stack346` builds and what items 347 and 349
used, so these deltas are row-comparable with theirs.  The shipped DEFAULT is still six; that
distinction is item 347's and it is kept here.

CONTROL that must succeed and exercises the suspect clause: the two TabICL arms must not be the
same vector.  They come from one runner, one checkpoint and one PCA cap, differing only in the
input key, so a bundle mix-up would give identical predictions and a spurious "sum".
"""
import itertools
import json
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd
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
JE = json.load(open(RES + "preds/k147_tabicl_E100.json"))     # E, pca 100
JS = json.load(open(RES + "preds/k147_tabicl_S100.json"))    # S, pca 100
acc = {}
for s in (0,1,2,3):
    fold, _ = butina_folds(list(rows.SMILES), seed=s); dg = fold_digest(fold)
    assert JE["digests"][s] == dg and JS["digests"][s] == dg
    pdz = members(s, XF, Y, M, fold, LO, HI)
    probe = [probe_oof(load_emb(rows), Y, M, fold, e)[0] for e in range(4)]
    chem = [np.asarray(v,float) for v in
            json.load(open(RES+"preds/k122_chemeleon_probe_seed%d.json"%s))["P"]]
    gbm  = [np.asarray(v,float) for v in
            json.load(open(RES+"preds/k114_heads_seed%d.json"%s))["heads"]["gbm"]]
    TE = [np.asarray(v,float) for v in JE["preds"]["E_seed%d" % s]]
    TS = [np.asarray(v,float) for v in JS["preds"]["S_seed%d" % s]]
    d = max(float(np.max(np.abs(TE[e]-TS[e]))) for e in range(4))
    assert d > 0.1, "КОНТРОЛЬ: армы E и S почти совпали, max|Δ| = %.2e -- перепутан вход" % d
    kept = [[Q[e] for k,Q in pdz if K._keep(e,k)] for e in range(4)]
    base = lambda e: kept[e] + [probe[e], chem[e], gbm[e]]
    arms = {"восемь (база)":      [np.mean(base(e),0)               for e in range(4)],
            "+E":                 [np.mean(base(e)+[TE[e]],0)       for e in range(4)],
            "+S":                 [np.mean(base(e)+[TS[e]],0)       for e in range(4)],
            "+E+S (сумма)":       [np.mean(base(e)+[TE[e],TS[e]],0) for e in range(4)],
            "+E+S, зонд долой":   [np.mean(kept[e]+[chem[e],gbm[e],TE[e],TS[e]],0)
                                   for e in range(4)]}
    for k, P in arms.items():
        acc.setdefault(k, []).append(score(P, Y, LO, HI, M, fold))
    cc = float(np.mean([abs(np.corrcoef(Y[M[:,e],e]-TE[e], Y[M[:,e],e]-TS[e])[0,1])
                        for e in range(4)]))
    acc.setdefault("_corrES", []).append(cc)
print("  КОНТРОЛЬ: армы E и S различны, max|Δ| = %.2f (порог 0.1)" % d)
print("  корреляция ошибок E с S: %.3f" % np.mean(acc.pop("_corrES")))
b  = np.mean([v[0] for v in acc["восемь (база)"]])
bp = np.mean([v[1] for v in acc["восемь (база)"]])
print("\n  %-20s %9s %8s %10s %7s %10s" % ("состав","ранг","sd","Δранг","знак","ΔST-RAE"))
for k, v in acc.items():
    r = [x[0] for x in v]; p = [x[1] for x in v]
    g = [x - y for x,y in zip(r, [q[0] for q in acc["восемь (база)"]])]
    print("  %-20s %9.4f %8.4f %+10.4f %5d/%d %+10.4f"
          % (k, np.mean(r), np.std(r,ddof=1), np.mean(r)-b,
             sum(1 for x in g if x>0), len(g), np.mean(p)-bp))
e_, s_ = (np.mean([x[0] for x in acc["+%s"%t]])-b for t in ("E","S"))
tot = np.mean([x[0] for x in acc["+E+S (сумма)"]])-b
print("\n  по отдельности %+.4f и %+.4f, арифметическая сумма %+.4f, измеренная %+.4f"
      % (e_, s_, e_+s_, tot))
print("  пол ранга 0.0036, пол ST-RAE 0.007")
json.dump(acc, open(RES+"logs/k148_tabiclsum.json","w"), default=float, indent=1)
