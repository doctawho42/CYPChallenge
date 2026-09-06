"""The TDI label, and the closed form of it that took 211 journal entries to write down.

The rule is stated piecewise below because that is how it was first found. It FOLDS:

    is_TDI  <=>  (Delta > log10 2)  AND  (pi_TDI > 4 + log10 2)          Delta = pi_TDI - pi_dir

Elementwise identical to the piecewise form on every labelled row with both arms, and both
reproduce the published label exactly: 2334/2334 on CYP3A4, 1493/1493 on CYP2D6. The proof is two
lines. If pi_dir > 4 then Delta > log10 2 already implies pi_TDI > 4.301, so the second conjunct is
slack and the first decides. If pi_dir <= 4 then pi_TDI > 4.301 already implies Delta > 0.301, so
the first is slack and the second decides. Not two cases -- one conjunction seen from two sides.

**Why this docstring exists.** The folded form makes visible what the piecewise form hides: the
shift is NECESSARY for every positive, and the second conjunct is a GATE that can only strike rows
out, never make one positive. Reading the mechanism off the piecewise form instead produced two
wrong accounts in a row (items 229 and 232) of where CYP3A4's classification score comes from, and
item 23 had already stated the conjunction in prose 211 entries earlier. The check that settles it
is one boolean comparison, and it now runs below every time this file does.

Consequences worth carrying, all measured (items 234, 238):

    фермент   ворота закрыты   MCC одних ворот   MCC одного сдвига   восстановлено воротами
    CYP3A4             0.398           +0.5667             +0.6875                   53.7 %
    CYP2D6             0.058           +0.1302             +0.9010                    1.2 %

CYP3A4's gate strikes out 39.8 per cent of rows and finding them is a pure potency question;
CYP2D6's strikes out 5.8 per cent, so its label is very nearly the shift alone. AUC 0.745 against
0.588 follows from that. And thresholding a SHRUNK prediction at the true threshold is the wrong
rule -- on CYP2D6 the out-of-fold optimal cut on the predicted arm is 4.69 to 4.98, not 4.301.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES
import pandas as pd, numpy as np
from sklearn.metrics import matthews_corrcoef
t=pd.read_csv(D+"cyp-challenge-TRAIN_TDI.csv")
for c in ["CYP3A4","CYP2D6"]:
    d=t[[f"{c}_is_TDI",f"{c}_pIC50_direct_inhibition",f"{c}_pIC50_TDI_condition"]].dropna()
    di=d[f"{c}_pIC50_direct_inhibition"]; td=d[f"{c}_pIC50_TDI_condition"]; y=d[f"{c}_is_TDI"].astype(bool)
    rule=np.where(di>4, (td-di)>np.log10(2), td>4.301)
    # The folded form, checked against the piecewise one on every run so the equivalence is a
    # fact of the repository rather than a claim in a docstring.
    conj=((td-di)>np.log10(2)) & (td>4.0+np.log10(2))
    assert (conj==rule).all(), f"{c}: свёрнутая форма разошлась с кусочной"
    gate=(td>4.0+np.log10(2)).to_numpy(); shift=((td-di)>np.log10(2)).to_numpy()
    assert not (y.to_numpy() & ~shift).any(), f"{c}: положительная метка без сдвига"
    assert not (y.to_numpy() & ~gate).any(), f"{c}: положительная метка при закрытых воротах"
    print(f"{c}: rule reproduces label in {np.mean(rule==y)*100:.2f}% of {len(d)} rows; MCC={matthews_corrcoef(y,rule):.3f}")
    mis=d[rule!=y]
    print("   mismatches:",len(mis))
    if len(mis): print(mis.head(5).round(3).to_string())
    # how much would a perfect direct-inhibition model + perfect TDI-arm model give?
    print(f"   base rate={y.mean():.3f};  label vs (direct pIC50>4.6): MCC={matthews_corrcoef(y,di>4.6):.3f}")
    print(f"   label vs (TDI-arm pIC50>4.8): MCC={matthews_corrcoef(y,td>4.8):.3f}")
    best=max(((matthews_corrcoef(y,td>th),th) for th in np.arange(4.0,6.0,0.05)))
    print(f"   best single-threshold on TDI-arm pIC50: MCC={best[0]:.3f} at {best[1]:.2f}")
    best=max(((matthews_corrcoef(y,(td-di)>th),th) for th in np.arange(0.0,1.2,0.02)))
    print(f"   best single-threshold on shift: MCC={best[0]:.3f} at {best[1]:.2f}")
    # Item 234: the conjunction, and the one number that explains the gap between the endpoints.
    print(f"   folded form agrees elementwise; gate closed on {1-gate.mean():.3f} of rows, "
          f"all negative; among open, positives {y[gate].mean():.3f}")
    print(f"   gate alone MCC={matthews_corrcoef(y,gate):.4f}, "
          f"shift alone MCC={matthews_corrcoef(y,shift):.4f}")
