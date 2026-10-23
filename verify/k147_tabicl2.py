"""Score the second TabICL run.  The question is the LAW, not the member.

Today's 1/n measurement says a member's worth is its solo advantage divided by the member count,
with no interaction term, measured on a candidate whose errors correlate 0.898 with the five.  The
one thing that could reopen the family is a DECORRELATED member showing a rising Delta*n.  S is
surrogate_adme, the only embedding on disk with no member standing on it inside the composition.
E at other PCA caps asks the separate question of whether solo rank is tunable at all.
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
from k106_probe import score, errcorr, members, probe_oof, load_emb
rows = pd.read_csv(D + "rows.csv"); z = np.load(D + "feats.npz")
XF = np.hstack([z["FP"], z["DESC"], z["MECH"]])
tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
        .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS],1)
LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS],1)
HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS],1)
M = ~np.isnan(Y); K.S.LO, K.S.HI = LO, HI
cands = {}
for tag, fn, key in (("TabICL S pca100","k147_tabicl_S100","S"), ("TabICL E pca50","k147_tabicl_E050","E"),
                     ("TabICL E pca200","k147_tabicl_E200","E"), ("TabICL E pca100","k147_tabicl_E100","E")):
    try: cands[tag] = (json.load(open(RES + "preds/%s.json" % fn)), key)
    except FileNotFoundError: print("  %s: файла нет, пропуск" % tag)
acc = {}
for s in (0,1,2,3):
    fold, _ = butina_folds(list(rows.SMILES), seed=s)
    dg = fold_digest(fold)
    pdz = members(s, XF, Y, M, fold, LO, HI)
    probe = [probe_oof(load_emb(rows), Y, M, fold, e)[0] for e in range(4)]
    chem = [np.asarray(v,float) for v in
            json.load(open(RES+"preds/k122_chemeleon_probe_seed%d.json"%s))["P"]]
    gbm  = [np.asarray(v,float) for v in
            json.load(open(RES+"preds/k114_heads_seed%d.json"%s))["heads"]["gbm"]]
    kept = [[Q[e] for k,Q in pdz if K._keep(e,k)] for e in range(4)]
    eight = [np.mean(kept[e]+[probe[e],chem[e],gbm[e]],0) for e in range(4)]
    s8 = score(eight, Y, LO, HI, M, fold)
    acc.setdefault("восемь (база)", []).append(dict(solo=s8[0], pair=s8[1], d9=0.0, corr5=np.nan,
                                                    corrp=np.nan, dn2=np.nan, dn8=np.nan))
    pool = [[Q[e] for k,Q in pdz] + [chem[e], gbm[e]] for e in range(4)]
    for tag, (J, key) in cands.items():
        assert J["digests"][s] == dg, "%s: разбиение не то" % tag
        T = [np.asarray(v,float) for v in J["preds"]["%s_seed%d" % (key, s)]]
        sT = score(T, Y, LO, HI, M, fold)
        nine = [np.mean(kept[e]+[probe[e],chem[e],gbm[e],T[e]],0) for e in range(4)]
        d9 = score(nine, Y, LO, HI, M, fold)[0] - s8[0]
        _p, vm, _v, _c = errcorr(T, pdz, Y, M)
        cp = float(np.mean([abs(np.corrcoef(Y[M[:,e],e]-T[e], Y[M[:,e],e]-probe[e])[0,1])
                            for e in range(4)]))
        # the law: substitute for the probe at n = 2 and n = 8, read Delta*n
        dn = {}
        for n in (2, 8):
            comp = [p[:n-1] for p in pool]
            a = score([np.mean(comp[e]+[probe[e]],0) for e in range(4)], Y,LO,HI,M,fold)[0]
            b = score([np.mean(comp[e]+[T[e]],0)     for e in range(4)], Y,LO,HI,M,fold)[0]
            dn[n] = (b - a) * n
        acc.setdefault(tag, []).append(dict(solo=sT[0], pair=sT[1], d9=d9,
                                            corr5=float(np.mean(vm)), corrp=cp,
                                            dn2=dn[2], dn8=dn[8]))
print("\n  %-17s %8s %8s %9s %7s %8s %8s %9s %9s"
      % ("арм","одиноч.","sd","Δ девятым","знак","corr 5","corr зонд","Δ·n n=2","Δ·n n=8"))
for tag, v in acc.items():
    f = lambda k: np.mean([x[k] for x in v])
    d = [x["d9"] for x in v]
    print("  %-17s %8.4f %8.4f %+9.4f %5d/%d %8.3f %9.3f %+9.4f %+9.4f"
          % (tag, f("solo"), np.std([x["solo"] for x in v], ddof=1), f("d9"),
             sum(1 for x in d if x>0), len(d), f("corr5"), f("corrp"), f("dn2"), f("dn8")))
print("\n  пол ранга 0.0036.  ЗАКОН: Δ·n постоянно по n и равно одиночному зазору над зондом.")
print("  зонд одиночный 0.6264; порог для девятого члена при n=8 -- одиночный ранг 0.6552")
json.dump(acc, open(RES+"logs/k147_tabicl2.json","w"), default=float, indent=1)
