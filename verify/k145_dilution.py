"""Does a member's contribution scale as 1/n?  Substitute TabICL for the probe in
compositions of n = 1..8 equally-weighted members and read Delta*n.

The companion pool is FIXED at all five members from members() regardless of _keep, so n is
the same for every enzyme.  That is not the shipped composition -- it is the instrument for
the scaling law.  The shipped arm (_keep, n = 6/7/8) is printed separately below.
"""
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
J = json.load(open(RES + "preds/k147_tabicl_E100.json"))
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
    # fixed companion pool, same length for every enzyme: five members, chemeleon, gbm
    pool = [[Q[e] for k,Q in pdz] for e in range(4)]
    for e in range(4):
        pool[e] = pool[e] + [chem[e], gbm[e]]
    assert all(len(p) == 7 for p in pool), [len(p) for p in pool]
    for n in range(1, 9):
        comp = [p[:n-1] for p in pool]
        a = score([np.mean(comp[e]+[probe[e]],0) for e in range(4)], Y, LO, HI, M, fold)[0]
        b = score([np.mean(comp[e]+[T[e]],0)     for e in range(4)], Y, LO, HI, M, fold)[0]
        acc.setdefault(n, []).append((a, b))
    # the shipped arm, _keep-admitted, for the record
    kept = [[Q[e] for k,Q in pdz if K._keep(e,k)] for e in range(4)]
    for nm, extra in (("шесть", []), ("семь", [chem]), ("восемь", [chem, gbm])):
        ex = lambda e: [x[e] for x in extra]
        a = score([np.mean(kept[e]+ex(e)+[probe[e]],0) for e in range(4)], Y,LO,HI,M,fold)[0]
        b = score([np.mean(kept[e]+ex(e)+[T[e]],0)     for e in range(4)], Y,LO,HI,M,fold)[0]
        acc.setdefault(nm, []).append((a, b))
print("  ЗАКОН РАЗБАВЛЕНИЯ: замена зонда на TabICL в составе из n равновесных членов")
print("  %3s %10s %10s %10s %7s %12s" % ("n","с зондом","с TabICL","Δранг","знак","Δ·n"))
for n in range(1, 9):
    v = acc[n]; d = [b-a for a,b in v]
    print("  %3d %10.4f %10.4f %+10.4f %5d/%d %12.4f"
          % (n, np.mean([a for a,_ in v]), np.mean([b for _,b in v]),
             np.mean(d), sum(1 for x in d if x>0), len(d), np.mean(d)*n))
print("\n  отгруженная ветка (_keep, число членов по ферментам [3,2,5,2]):")
for nm, n in (("шесть",6), ("семь",7), ("восемь",8)):
    v = acc[nm]; d = [b-a for a,b in v]
    print("  %6s %10.4f %10.4f %+10.4f %5d/%d %12.4f"
          % (nm, np.mean([a for a,_ in v]), np.mean([b for _,b in v]),
             np.mean(d), sum(1 for x in d if x>0), len(d), np.mean(d)*n))
g = np.mean([b-a for a,b in acc[1]])
print("\n  одиночный зазор g = %+.4f; предсказание закона Δ ≈ g/n" % g)
print("  чтобы Δ превысил пол 0.0036 при n=8, нужен зазор g = %.4f, то есть член с рангом %.4f"
      % (0.0036*8, np.mean([a for a,_ in acc[1]]) + 0.0036*8))
print("  состав из восьми сам имеет ранг %.4f" % np.mean([a for a,_ in acc["восемь"]]))
json.dump({str(k): v for k,v in acc.items()}, open(RES+"logs/k145_dilution.json","w"), indent=1)
