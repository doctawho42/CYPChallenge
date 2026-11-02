"""All 32 subsets of the five optional slots on top of the kept five.  Seeds 0-3.

This is SELECTION: picking the best of 32 on the same four seeds the floor was measured on is
exactly what item 316's fitted weights and item 280's Caruana selector did before they failed to
transfer.  So the table is the finding and the winner is a HYPOTHESIS, to be confirmed on seeds
4-7, which this script never reads.  Printed sorted, with the shipped arms marked.
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
JE = json.load(open(RES + "preds/k147_tabicl_E100.json")); JS = json.load(open(RES + "preds/k147_tabicl_S100.json"))
NAMES = ["зонд", "chemeleon", "gbm", "TabICL-E", "TabICL-S"]
acc = {}
for s in (0,1,2,3):
    fold, _ = butina_folds(list(rows.SMILES), seed=s); dg = fold_digest(fold)
    assert JE["digests"][s] == dg and JS["digests"][s] == dg
    pdz = members(s, XF, Y, M, fold, LO, HI)
    opt = [[probe_oof(load_emb(rows), Y, M, fold, e)[0] for e in range(4)],
           [np.asarray(v,float) for v in
            json.load(open(RES+"preds/k122_chemeleon_probe_seed%d.json"%s))["P"]],
           [np.asarray(v,float) for v in
            json.load(open(RES+"preds/k114_heads_seed%d.json"%s))["heads"]["gbm"]],
           [np.asarray(v,float) for v in JE["preds"]["E_seed%d"%s]],
           [np.asarray(v,float) for v in JS["preds"]["S_seed%d"%s]]]
    kept = [[Q[e] for k,Q in pdz if K._keep(e,k)] for e in range(4)]
    for r in range(6):
        for combo in itertools.combinations(range(5), r):
            P = [np.mean(kept[e] + [opt[i][e] for i in combo], 0) for e in range(4)]
            acc.setdefault(combo, []).append(score(P, Y, LO, HI, M, fold))
SIX   = (0,)                 # the shipped default: kept five + probe
EIGHT = (0, 1, 2)            # item 346's stack, what --stack346 builds
b  = np.mean([v[0] for v in acc[EIGHT]]); bp = np.mean([v[1] for v in acc[EIGHT]])
s6 = np.mean([v[0] for v in acc[SIX]]);   p6 = np.mean([v[1] for v in acc[SIX]])
print("  отгруженный ПО УМОЛЧАНИЮ состав из шести:  ранг %.4f" % s6)
print("  восемь пункта 346 (--stack346), база ниже:  ранг %.4f\n" % b)
print("  %-42s %4s %8s %8s %9s %7s %9s %9s"
      % ("слоты поверх пяти оставленных","n","ранг","sd","Δк восьми","знак","Δк шести","ΔST-RAE"))
rk = sorted(acc.items(), key=lambda kv: -np.mean([x[0] for x in kv[1]]))
for combo, v in rk[:12]:
    r = [x[0] for x in v]; p = [x[1] for x in v]
    g = [x - y for x,y in zip(r, [q[0] for q in acc[EIGHT]])]
    mark = "  <- по умолчанию" if combo == SIX else ("  <- пункт 346" if combo == EIGHT else "")
    print("  %-42s %4d %8.4f %8.4f %+9.4f %5d/%d %+9.4f %+9.4f%s"
          % ("+".join(NAMES[i] for i in combo) or "(ничего)", 5+len(combo),
             np.mean(r), np.std(r,ddof=1), np.mean(r)-b,
             sum(1 for x in g if x>0), len(g), np.mean(r)-s6, np.mean(p)-bp, mark))
print("  ...")
for combo in (SIX, EIGHT):
    v = acc[combo]; r = [x[0] for x in v]
    print("  %-42s %4d %8.4f %8.4f %+9.4f %5s %+9.4f %+9.4f%s"
          % ("+".join(NAMES[i] for i in combo), 5+len(combo), np.mean(r),
             np.std(r,ddof=1), np.mean(r)-b, "", np.mean(r)-s6,
             np.mean([x[1] for x in v])-bp,
             "  <- по умолчанию" if combo == SIX else "  <- пункт 346"))
print("\n  пол ранга 0.0036, пол ST-RAE 0.007; лучший из 32 -- ГИПОТЕЗА, проверка на сидах 4-7")
json.dump({",".join(map(str,k)): v for k,v in acc.items()},
          open(RES+"logs/k149_sweep.json","w"), default=float, indent=1)
