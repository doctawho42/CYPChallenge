"""Where is rank exchange possible at all? Three oracles and the one number item 133 skipped.

The affine pair is a strictly increasing map, so it preserves Spearman exactly (item 77, and
the reason "measure rank" is the rule). A correction that is monotone WITHIN a group is
likewise order-preserving inside that group, so it can only move rank BETWEEN groups. That
splits every achievable rank gain into two disjoint parts, and this file bounds each.

    A   true order INSIDE groups, model order BETWEEN     ceiling of within-group permutation
    B   true group LEVELS, model order inside             ceiling of group-wise shifts
    C   how much of the discordance lives on pairs the ensemble members disagree about

Item 128 is the precedent and the procedure: bound a class of proposals with one oracle before
building any member of it. It was applied to the metric numerator and never to rank, though
rank is the criterion that survives post-processing.

**And the number item 133 did not compute.** That item measured the model's within-series
spread against tau, what chemistry allows, and found ratios of 0.29 to 0.58 -- the model is two
to three times over-smoothed inside series. It read this as "shrinking would compound the
error" and closed the series layer. The opposite reading is that expansion is licensed, and it
is not equivalent: **a scale ratio cannot tell a compressed signal from noise.** Deviations
that are pure noise would produce exactly the same ratio, and expanding them would amplify
noise into the global ordering.

The quantity that separates the two is the within-group regression slope of the true deviation
on the model's,

    beta = rho_within * sd_true / sd_model,

the classical attenuation factor. Expansion by k is licensed only where beta > 1. Item 133's
ratios put sd_true / sd_model at 1.3 to 2.9, so the threshold on rho_within is roughly 0.35 to
0.77 depending on enzyme -- a demanding bar, and one nobody has checked.

Groups are analog series: Butina clusters of the training molecules at a tighter threshold than
the cross-validation split uses, matching the 0.60-0.85 similarity band k34_series.py works in.
The split's own 0.35 clusters are far too coarse to be analog series.

Reads results/preds/oof.json and data/rows.csv. Writes nothing. Minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator
from rdkit.ML.Cluster import Butina

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]


def series(smiles, cut):
    """Analog series: Butina at a tight cutoff. Same generator settings as cypsplit."""
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in smiles]
    n = len(fps)
    d = []
    for i in range(1, n):
        d.extend(1.0 - np.array(DataStructs.BulkTanimotoSimilarity(fps[i], fps[:i])))
    cl = Butina.ClusterData(d, n, cut, isDistData=True)
    g = np.zeros(n, int)
    for k, c in enumerate(cl):
        for i in c:
            g[i] = k
    return g


def oracle_A(p, y, g):
    """Model order between groups, TRUE order inside. The model's own values are reassigned
    within each group by the true ranking, so the between-group structure and the value
    distribution are untouched and only who-gets-which changes."""
    q = p.copy()
    for k in np.unique(g):
        i = np.where(g == k)[0]
        if len(i) < 2:
            continue
        q[i[np.argsort(np.argsort(y[i]))]] = np.sort(p[i])
    return q


def oracle_B(p, y, g):
    """TRUE group levels, model order inside: the group mean is replaced by the true one and
    the model's deviations are carried over unchanged."""
    q = p.copy()
    for k in np.unique(g):
        i = np.where(g == k)[0]
        q[i] = y[i].mean() + (p[i] - p[i].mean())
    return q


def within_beta(p, y, g):
    """Attenuation slope of true deviation on model deviation, pooled over groups of size>=2."""
    dp, dy = [], []
    for k in np.unique(g):
        i = np.where(g == k)[0]
        if len(i) < 2:
            continue
        dp.append(p[i] - p[i].mean())
        dy.append(y[i] - y[i].mean())
    if not dp:
        return np.nan, np.nan, np.nan, 0
    dp, dy = np.concatenate(dp), np.concatenate(dy)
    if dp.std() < 1e-9:
        return np.nan, np.nan, np.nan, len(dp)
    rho = float(np.corrcoef(dp, dy)[0, 1])
    return float(dy.std() / dp.std()), rho, rho * float(dy.std() / dp.std()), len(dp)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cut", type=float, default=0.40,
                    help="расстояние Бутины: 0.40 = сходство 0.60, нижний край полосы k34")
    ap.add_argument("--key", default="FP+DESC+MECH")
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    oof = json.load(open(RES + "preds/oof.json"))

    print(f"строю серии при расстоянии {a.cut} (сходство {1-a.cut:.2f})", flush=True)
    g_all = series(list(rows.SMILES), a.cut)
    sz = np.bincount(g_all)
    print(f"  {len(sz)} серий на {len(rows)} молекул; "
          f"в сериях >=2 членов: {int(sz[sz >= 2].sum())} молекул в {int((sz >= 2).sum())} сериях\n",
          flush=True)

    T = []
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy()
        p = np.asarray(oof[f"{a.key}|{c}"], float)
        g = g_all[m]
        # переиндексация групп внутри маски + сколько молекул реально в парах
        u, g = np.unique(g, return_inverse=True)
        s = np.bincount(g)
        paired = int(s[s >= 2].sum())

        base = float(spearmanr(y, p).statistic)
        rA = float(spearmanr(y, oracle_A(p, y, g)).statistic)
        rB = float(spearmanr(y, oracle_B(p, y, g)).statistic)
        sc, rho_w, beta, npair = within_beta(p, y, g)
        T.append({"фермент": c, "n": int(m.sum()), "в сериях": paired,
                  "база": base, "оракул A": rA, "оракул B": rB,
                  "A-база": rA - base, "B-база": rB - base,
                  "sd_ист/sd_мод": sc, "rho внутри": rho_w, "beta": beta})

    df = pd.DataFrame(T).set_index("фермент")
    print(df[["n", "в сериях", "база", "оракул A", "оракул B", "A-база", "B-база"]]
          .round(4).to_string())
    print()
    print(df[["sd_ист/sd_мод", "rho внутри", "beta"]].round(3).to_string())
    print(f"\nмакро: A-база {df['A-база'].mean():+.4f}, B-база {df['B-база'].mean():+.4f}")
    print(f"макро beta {df['beta'].mean():.3f}")
    print("""
Как читать. A и B --- ПОТОЛКИ, недостижимые по построению: в них подставлена истина. Опыт
этого файла (пункты 126, 128) в том, что реальное составляет от потолка малую долю, поэтому
маленький потолок закрывает направление, а большой его лишь не закрывает.

beta --- единственное число здесь, которое НЕ потолок, а прямая проверка. Расширение внутри
серии на коэффициент k монотонно внутри группы, поэтому внутрисерийный порядок оно не меняет
вовсе; вся его польза межсерийная и существует ровно тогда, когда внутрисерийные отклонения
модели информативны. beta > 1 --- расширение лицензировано, оптимальный k около beta.
beta < 1 --- отклонения зашумлены сильнее, чем сжаты, и расширение усилит шум: тогда пункт 133
закрыл серийную ось целиком, а не только стягивание, и читать его можно было в одну сторону.""")


if __name__ == "__main__":
    main()
