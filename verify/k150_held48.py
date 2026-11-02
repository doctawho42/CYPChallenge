"""Confirmation on seeds 4-7, which the 32-arm sweep never read.

The sweep picked its winner on seeds 0-3, which is selection; items 280 and 316 both measured
things that won in-sample and did not transfer.  These four seeds are the out-of-sample test.  Only
the arms needing no pending cache are here: the shipped six, and the kept five with one or both
TabICL arms.  The arms that also want the chemeleon probe and gbm600 follow when those caches fill.

CONTROL that must succeed: the fold digest of each held-out seed must match the digest the box
stamped into the prediction file, or the predictions describe a different split.
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
J = json.load(open(RES + "preds/k150_tabicl_held48.json"))
SEEDS = (4,5,6,7)
acc = {}
for i, s in enumerate(SEEDS):
    fold, _ = butina_folds(list(rows.SMILES), seed=s); dg = fold_digest(fold)
    assert J["digests"][i] == dg, "сид %d: %s в файле против %s здесь" % (s, J["digests"][i], dg)
    pdz = members(s, XF, Y, M, fold, LO, HI)
    probe = [probe_oof(load_emb(rows), Y, M, fold, e)[0] for e in range(4)]
    chem = [np.asarray(v,float) for v in
            json.load(open(RES+"preds/k122_chemeleon_probe_seed%d.json"%s))["P"]]
    TE = [np.asarray(v,float) for v in J["preds"]["E_seed%d" % s]]
    TS = [np.asarray(v,float) for v in J["preds"]["S_seed%d" % s]]
    kept = [[Q[e] for k,Q in pdz if K._keep(e,k)] for e in range(4)]
    arms = {"шесть: пять+зонд":       [np.mean(kept[e]+[probe[e]],0)        for e in range(4)],
            "пять+TabICL-E":          [np.mean(kept[e]+[TE[e]],0)           for e in range(4)],
            "пять+TabICL-S":          [np.mean(kept[e]+[TS[e]],0)           for e in range(4)],
            "пять+E+S":               [np.mean(kept[e]+[TE[e],TS[e]],0)     for e in range(4)],
            "шесть+E+S":              [np.mean(kept[e]+[probe[e],TE[e],TS[e]],0) for e in range(4)],
            "пять+chem+E+S":          [np.mean(kept[e]+[chem[e],TE[e],TS[e]],0) for e in range(4)],
            "шесть+chem+E+S":         [np.mean(kept[e]+[probe[e],chem[e],TE[e],TS[e]],0)
                                       for e in range(4)],
            "пять+chem+E":            [np.mean(kept[e]+[chem[e],TE[e]],0) for e in range(4)],
            "пять+chem+S":            [np.mean(kept[e]+[chem[e],TS[e]],0) for e in range(4)],
            "одиночный TabICL-E":     TE,
            "одиночный TabICL-S":     TS,
            "одиночный зонд":         probe}
    for k, P in arms.items():
        acc.setdefault(k, []).append(score(P, Y, LO, HI, M, fold))
print("  КОНТРОЛЬ: дайджесты сидов 4-7 совпали с теми, что проштампованы в файле предсказаний")
b  = np.mean([v[0] for v in acc["шесть: пять+зонд"]])
bp = np.mean([v[1] for v in acc["шесть: пять+зонд"]])
print("\n  ОТЛОЖЕННЫЕ СИДЫ 4-7, база -- отгруженные по умолчанию шесть (ранг %.4f)\n" % b)
print("  %-22s %9s %8s %10s %7s %10s" % ("арм","ранг","sd","Δк шести","знак","ΔST-RAE"))
for k, v in acc.items():
    r = [x[0] for x in v]; p = [x[1] for x in v]
    g = [x - y for x,y in zip(r, [q[0] for q in acc["шесть: пять+зонд"]])]
    print("  %-22s %9.4f %8.4f %+10.4f %5d/%d %+10.4f"
          % (k, np.mean(r), np.std(r,ddof=1), np.mean(r)-b,
             sum(1 for x in g if x>0), len(g), np.mean(p)-bp))
print("\n  для сравнения, те же армы на сидах 0-3 (из k149_sweep.json):")
sw = json.load(open(RES+"logs/k149_sweep.json"))
m = {"шесть: пять+зонд": "0", "пять+TabICL-E": "3", "пять+TabICL-S": "4",
     "пять+E+S": "3,4", "шесть+E+S": "0,3,4", "пять+chem+E+S": "1,3,4",
     "шесть+chem+E+S": "0,1,3,4", "пять+chem+E": "1,3", "пять+chem+S": "1,4"}
b0 = np.mean([x[0] for x in sw["0"]])
for k, key in m.items():
    r = [x[0] for x in sw[key]]
    print("    %-22s %9.4f  Δк шести %+.4f" % (k, np.mean(r), np.mean(r)-b0))
json.dump(acc, open(RES+"logs/k150_held48.json","w"), default=float, indent=1)
