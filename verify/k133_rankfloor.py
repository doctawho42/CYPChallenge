"""The exact ST-RAE floor at FIXED ordering, and it is what corrects item 354's headline.

WHY. Item 354 read 0.0872 of macro ST-RAE off the board as "unclaimed at the ordering we already
have", by matching on Spearman: every entrant within 0.03 of our rho on an enzyme, best ST-RAE among
them. That number is an upper bound and it conflates two things. ST-RAE under the BEST monotone
placement depends only on the rank PERMUTATION, and matching Spearman does not match the
permutation --- two predictors at rho 0.75 can order the compounds differently, and one ordering can
admit a much lower floor than another by putting its errors where the bands are wide. So part of the
0.0872 is placement skill we lack and part is a permutation we do not have.

WHAT THIS COMPUTES. The floor itself, exactly, on our own labels: the weakly monotone placement of
our out-of-fold predictions that minimises the interval hinge. That is isotonic regression under a
convex pointwise loss, so PAVA applies with the block optimum solved in closed form --- for the
hinge the derivative sum w_i(1[c>hi_i] - 1[c<lo_i]) is non-decreasing in c, so the block optimum is
the first zero crossing over the pooled band edges, found by one sort.

    фермент    как есть   пол, w=1   пол, широкая маргиналь
    CYP1A2       0.7594     0.7131            0.7554
    CYP2C9       0.5549     0.5016            0.5175
    CYP2D6       0.8403     0.8157            0.9100
    CYP3A4       0.4135     0.3836            0.3837
    МАКРО        0.6421     0.6035            0.6416

**So placement's ENTIRE headroom at our ordering, on our own label distribution, is 0.0386 of macro
--- in sample.** Item 316 already measured what happens out of fold: a monotone piecewise-linear
placement is 0.0011 to 0.0014 WORSE than the affine pair at every bin count from 8 to 128, sign 1/4,
and the in-sample 0.0223 it quoted was pure selection contamination. So the realistic expected gain
from placement alone is at or below that 0.0386, and previously it has been negative.

WHAT THIS DOES NOT SETTLE, and it is the reason item 354's finding is narrowed rather than killed.
This floor is computed on OUR label distribution, where there is no distribution mismatch by
construction --- so it cannot see the thing item 354 actually found, which is that the test's labels
are 1.02 to 1.48 times wider than ours and the shipped transform is fitted to the narrow one. A
mismatch correction is not in-sample placement. The honest reading of the two together: the scale
correction is worth something between nothing and 0.0872, the in-sample floor says the placement
part of it is at most 0.0386, and only the board can say which.

Also worth reading off the table: our board score of 0.5357 is BELOW our own out-of-fold floor of
0.6035, so the test is an easier sample than our cross-validation and no out-of-fold floor transfers
to the board arithmetically.

    uv run python verify/k133_rankfloor.py

Reads results/preds/oof_submitted.json, results/logs/k131_scaleshift.json, data/rows.csv and the
inhibition table. Prints only; the numbers above are the record.
"""
import sys, pathlib
R=pathlib.Path("/home/user/CYPChallenge")
sys.path.insert(0,str(R)); sys.path.insert(0,str(R/"src")); sys.path.insert(0,str(R/"verify"))
from cyppaths import D, RES
import json, numpy as np, pandas as pd
from submit import CYPS
from reweight import widen
from k106_probe import strae

def block_opt(lo, hi, w):
    """argmin_c sum w_i*(max(c-hi_i,0)+max(lo_i-c,0)).  Derivative is non-decreasing:
    g(c) = sum w_i[c>hi_i] - sum w_i[c<lo_i], from -W to +W. Take the first zero crossing."""
    ev = np.concatenate([np.stack([lo, w], 1), np.stack([hi, w], 1)])
    ev = ev[np.argsort(ev[:, 0], kind="stable")]
    g = -w.sum() + np.cumsum(ev[:, 1])
    j = int(np.searchsorted(g, 0.0, side="left"))
    return float(ev[min(j, len(ev)-1), 0])

def iso_hinge(order, lo, hi, w):
    """Monotone-in-`order` placement minimising the interval hinge. PAVA with the block
    optimum above; blocks keep their member index lists so a merge is exact."""
    idx = np.argsort(order, kind="stable")
    blocks = []            # (indices, c)
    for i in idx:
        cur = [i]
        c = block_opt(lo[cur], hi[cur], w[cur])
        while blocks and blocks[-1][1] >= c - 1e-12:
            prev, _ = blocks.pop()
            cur = prev + cur
            c = block_opt(lo[cur], hi[cur], w[cur])
        blocks.append((cur, c))
    q = np.empty(len(order))
    for members, c in blocks:
        q[members] = c
    return q

rows=pd.read_csv(D+"rows.csv")
tr=(pd.read_csv(D+"cyp-challenge-TRAIN_inhibition.csv")
    .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
Y =np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS],1)
LO=np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS],1)
HI=np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS],1)
M=~np.isnan(Y)
P=[np.asarray(v,float) for v in json.load(open(RES+"preds/oof_submitted.json"))["P"]]
k=json.load(open(RES+"logs/k131_scaleshift.json"))
ratio={c: k["spread"][c]["bound"]/k["spread"][c]["train_sd"] for c in CYPS}

print("Пол ST-RAE при НЕИЗМЕННОМ порядке (лучшая монотонная раскладка), в выборке\n")
print("%-8s %9s %9s %9s | %9s %9s" % ("фермент","как есть","пол, w=1","пол, шир.","доска: наш","лучший при ρ"))
board={"CYP1A2":(0.5656,0.4955),"CYP2C9":(0.4258,0.3620),"CYP2D6":(0.6786,0.5357),"CYP3A4":(0.4729,0.4008)}
acc=[[],[],[]]
for e,c in enumerate(CYPS):
    m=M[:,e]; y,lo,hi,p=Y[m,e],LO[m,e],HI[m,e],P[e]
    u=np.ones(len(y))/len(y)
    cur=strae(y,p,y_true_upper=hi,y_true_lower=lo)
    q1=iso_hinge(p,lo,hi,u); f1=strae(y,q1,y_true_upper=hi,y_true_lower=lo)
    ww=widen(y,ratio[c]); ww=ww/ww.sum()
    q2=iso_hinge(p,lo,hi,ww); f2=strae(y,q2,y_true_upper=hi,y_true_lower=lo)
    # rank must be preserved (weakly: monotone, ties allowed)
    from scipy.stats import spearmanr
    r0,r1=spearmanr(y,p).statistic, spearmanr(y,q1).statistic
    acc[0].append(cur); acc[1].append(f1); acc[2].append(f2)
    print("%-8s %9.4f %9.4f %9.4f | %9.4f %9.4f   (ρ %.4f -> %.4f)"
          % (c,cur,f1,f2,board[c][0],board[c][1],r0,r1))
print("%-8s %9.4f %9.4f %9.4f | %9.4f %9.4f"
      % ("МАКРО",np.mean(acc[0]),np.mean(acc[1]),np.mean(acc[2]),
         np.mean([v[0] for v in board.values()]),np.mean([v[1] for v in board.values()])))
print("\nпервое место на доске: 0.3706")
