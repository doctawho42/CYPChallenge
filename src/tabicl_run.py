"""An in-context tabular learner as a member, on two representations. Writes files, scores nothing.

Item 335 closed the tabular-foundation axis BY PROXY and said so in its own first sentence: three
non-linear heads on the cached 600 columns "ask the same question" as a foundation head, for twelve
CPU-minutes "instead of a dependency and a build". The proxy is the weak step. A gradient-boosting
machine is not a stand-in for a transformer pretrained on millions of synthetic tabular tasks doing
approximate Bayesian inference in context, whose advantage is specifically the small-n regime --
1285 to 2335 rows per enzyme, which is ours. The reason item 335 gave for the proxy was COST, and a
GPU has removed it.

TabICL rather than TabPFN: TabPFN 9.2 gates its weights behind a registered account and a license
acceptance, which is the user's act and not this file's. TabICL is the same family by different
authors and downloads its checkpoint without one.

TWO ARMS, because item 349 measured that two learners on the same target class make the same errors,
and the composition already holds a ridge AND a boosting head on the embedding:

    E  --  на замороженном вложении (600 -> PCA 100), рецепт места 41
    X  --  на НАШЕЙ матрице (2295 -> PCA 100), где живут пять членов и ни один не in-context

PCA-capping is place 41's own control, in its words: "with only one to two thousand labelled rows
per endpoint, an uncapped concatenation of several thousand embedding dimensions lets those models
memorise the fold rather than learn from it." Scaler and PCA are fitted on the FITTING rows of each
fold only; doing either on the union is item 202's defect 2.

A ridge on the identical post-PCA input is fitted in the same loop, so the comparison is the LEARNER
and not the preprocessing -- without it a difference could be the PCA.

THE RECIPE, written out because a log whose script cannot be re-run is item 345's k123 defect and a
dependency that lives only in one shell's history is the same defect one step out. Measured on the
GPU box, 2026-10-23:

    uv pip install --python <env>/bin/python tabicl          # 2.2.0, no account, no licence

    tabicl 2.2.0   torch 2.13.0+cu130   scikit-learn 1.9.1   numpy 2.5.4   scipy 1.18.1   py3.12

The checkpoint downloads itself on first use from Hugging Face `jingang/TabICL`, snapshot
4dcd344ece2c00be9e831fdd35bed57b5ad83e19, file `tabicl-regressor-v2-20260212.ckpt`, **114324594
bytes**. No token and no accepted licence are involved, which is the whole reason this family is
TabICL and not TabPFN: TabPFN 9.2 gates its weights behind a registered account and an accepted
licence, and accepting one is the user's act, not this repository's.

**Note the scikit-learn: 1.9.1, not the pinned 1.3.2.** That is why this file writes a file and
scores nothing. `CLAUDE.md` records that 1.9.0 moves every prediction in this repository by ~0.135
pIC50 and doubles the measured effect of the mechanistic block, so every number that reaches the
journal is computed under the pin, by `verify/k147_tabicl2.py` and the scripts around it.
"""
import argparse, json, pathlib, time
import numpy as np


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--inp", default="/tmp/d42/data/tabpfn_in.npz")
    ap.add_argument("--out", default="/tmp/d42/out")
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--arms", default="E,X")
    ap.add_argument("--pca", type=int, default=100)
    ap.add_argument("--device", default="cuda")
    # Held-out seeds arrive as a separate 40 KB file rather than a second copy of the
    # 46 MB bundle: its fold<N> keys are merged in, so --seeds can name them.
    ap.add_argument("--extra", default="")
    a = ap.parse_args()
    from sklearn.decomposition import PCA
    from sklearn.preprocessing import StandardScaler
    from sklearn.linear_model import RidgeCV
    from tabicl import TabICLRegressor
    import torch

    z = np.load(a.inp, allow_pickle=True)
    Y = z["Y"].astype(np.float64)
    M = np.isfinite(Y)
    CY = [str(c) for c in z["cyps"]]
    # Any 2-D key that is not a bundle fixture is a representation, so a second bundle can
    # carry other embeddings without touching this file again.
    FIX = {"Y", "cyps", "digests"}
    reps = {k: z[k].astype(np.float64) for k in z.files
            if k not in FIX and not k.startswith("fold") and z[k].ndim == 2}
    EX = dict(np.load(a.extra, allow_pickle=True)) if a.extra else {}
    digs = [str(d) for d in z["digests"]]
    if "digests48" in EX:
        digs = [str(d) for d in EX["digests48"]]
    out, ridge_out = {}, {}
    for s in [int(x) for x in a.seeds.split(",") if x.strip()]:
        fold = EX["fold%d" % s] if ("fold%d" % s) in EX else z["fold%d" % s]
        for arm in [x for x in a.arms.split(",") if x.strip()]:
            R = reps[arm]
            P, Q = [], []
            t0 = time.time()
            for e in range(len(CY)):
                m = M[:, e]
                Rm, ym, fm = R[m], Y[m, e], fold[m]
                p = np.full(int(m.sum()), np.nan)
                q = np.full(int(m.sum()), np.nan)
                for f in range(5):
                    trn, te = fm != f, fm == f
                    if te.sum() == 0:
                        continue
                    sc = StandardScaler().fit(Rm[trn])
                    A, B = sc.transform(Rm[trn]), sc.transform(Rm[te])
                    k = min(a.pca, A.shape[1], A.shape[0] - 1)
                    pc = PCA(n_components=k, random_state=0).fit(A)
                    A, B = pc.transform(A), pc.transform(B)
                    mdl = TabICLRegressor(device=a.device, random_state=0)
                    mdl.fit(A, ym[trn])
                    p[te] = mdl.predict(B)
                    q[te] = RidgeCV(alphas=[100., 1000., 10000.]).fit(A, ym[trn]).predict(B)
                    del mdl
                    torch.cuda.empty_cache()
                if not np.isfinite(p).all() or not np.isfinite(q).all():
                    raise SystemExit("%s/%s: незаполненные строки" % (arm, CY[e]))
                P.append(p); Q.append(q)
                print("    %s %s %s за %.0f с" % (arm, CY[e], "сид%d" % s, time.time() - t0),
                      flush=True)
            out["%s_seed%d" % (arm, s)] = [v.tolist() for v in P]
            ridge_out["%s_seed%d" % (arm, s)] = [v.tolist() for v in Q]
            print("  арм %s сид %d готов, %.0f с" % (arm, s, time.time() - t0), flush=True)
    pathlib.Path(a.out).mkdir(parents=True, exist_ok=True)
    fn = "%s/tabicl_oof.json" % a.out
    json.dump({"preds": out, "ridge_same_input": ridge_out, "pca": a.pca, "digests": digs},
              open(fn, "w"))
    print("\n  записано %s" % fn)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
