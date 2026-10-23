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
from k106_probe import score, errcorr, members, probe_oof, load_emb
rows = pd.read_csv(D + "rows.csv"); z = np.load(D + "feats.npz")
XF = np.hstack([z["FP"], z["DESC"], z["MECH"]])
tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
        .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS],1)
LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS],1)
HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS],1)
M = ~np.isnan(Y); K.S.LO, K.S.HI = LO, HI
J = json.load(open(RES + "preds/k147_tabicl_E100.json"))
rows_out = {}
for s in (0,1,2,3):
    fold, _ = butina_folds(list(rows.SMILES), seed=s); dg = fold_digest(fold)
    assert dg == J["digests"][s]
    pdz = members(s, XF, Y, M, fold, LO, HI)
    probe = [probe_oof(load_emb(rows), Y, M, fold, e)[0] for e in range(len(CYPS))]
    cj = json.load(open(RES+"preds/k122_chemeleon_probe_seed%d.json"%s))
    hj = json.load(open(RES+"preds/k114_heads_seed%d.json"%s))
    chem = [np.asarray(v,float) for v in cj["P"]]; gbm = [np.asarray(v,float) for v in hj["heads"]["gbm"]]
    T = [np.asarray(v,float) for v in J["preds"]["E_seed%d"%s]]
    kept = [[Q[e] for k,Q in pdz if K._keep(e,k)] for e in range(len(CYPS))]
    arms = {
      "восемь (база)":            [np.mean(kept[e]+[probe[e],chem[e],gbm[e]],0) for e in range(4)],
      "ЗАМЕНА зонда на TabICL":   [np.mean(kept[e]+[T[e],chem[e],gbm[e]],0) for e in range(4)],
      "девятым (добавление)":     [np.mean(kept[e]+[probe[e],chem[e],gbm[e],T[e]],0) for e in range(4)],
      "замена + добавление зонда":[np.mean(kept[e]+[T[e],chem[e],gbm[e],probe[e]],0) for e in range(4)],
      "ЗАМЕНА chemeleon":         [np.mean(kept[e]+[probe[e],T[e],gbm[e]],0) for e in range(4)],
      "ЗАМЕНА gbm":               [np.mean(kept[e]+[probe[e],chem[e],T[e]],0) for e in range(4)],
    }
    for k,P in arms.items():
        sc = score(P, Y, LO, HI, M, fold)
        rows_out.setdefault(k, []).append({"rank": sc[0], "pair": sc[1]})
b = np.mean([x["rank"] for x in rows_out["восемь (база)"]])
bp = np.mean([x["pair"] for x in rows_out["восемь (база)"]])
print("  %-28s %9s %8s %10s %7s %10s" % ("состав","ранг","sd","Δранг","знак","Δпара"))
for k,v in rows_out.items():
    r = [x["rank"] for x in v]; p = [x["pair"] for x in v]
    g = [x - y for x,y in zip(r, [x["rank"] for x in rows_out["восемь (база)"]])]
    print("  %-28s %9.4f %8.4f %+10.4f %5d/%d %+10.4f"
          % (k, np.mean(r), np.std(r,ddof=1), np.mean(r)-b,
             sum(1 for x in g if x>0), len(g), np.mean(p)-bp))
print("\n  пол ранга 0.0036, пол ST-RAE 0.007")
json.dump(rows_out, open(RES+"logs/k144_subst.json","w"), ensure_ascii=False, indent=1)
