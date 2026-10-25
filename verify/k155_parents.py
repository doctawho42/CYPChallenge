"""Recover the test set's construction from the data, and check it against what the organisers
published.

The challenge's About tab states the design in one sentence: "Hit expansion was performed on the top
25 hits from CYP1A2, CYP2C9, and CYP3A4 by purchasing the top 10 chemisimilars (ECFP4 Tanimoto
similarity) per hit compound from the Enamine US catalog", giving 25 x 3 x 10 = 750 compounds
"fully assayed in 12-point DRC mode against all four CYPs".

If that is true of the data in front of us, then three things must hold, and each is falsifiable:

  1. every test compound's nearest TRAINING neighbour by ECFP4 is its parent, so the 750 nearest
     neighbours concentrate on about 75 distinct training compounds rather than spreading over
     thousands;
  2. those parents are the top 25 by direct-inhibition pIC50 for CYP1A2, CYP2C9 and CYP3A4, and
     NOT for CYP2D6, which the page excludes from the expansion;
  3. each parent carries about ten children.

Item 376 measured that the test sits at a median nearest-neighbour similarity of 0.5873 against
0.4500 for leave-one-out on our own data, and read it as generic density. If this file's claims
hold, that density is not generic at all: it is 75 tight series, which is a different object and
carries a different lever -- the test's ACTIVITY distribution becomes derivable from the design
rather than guessed, and item 377 showed that guessing it wrong cost 71 places.

THE NULL THAT MAKES CONCENTRATION READABLE. 750 random training compounds, matched to their nearest
OTHER training compound, must NOT concentrate on 75 parents. Without it, "the neighbours concentrate"
could be an artefact of any clustered chemical library.

THE CONTROL THAT MUST SUCCEED. Every training compound's maximum similarity to the full training
set INCLUDING itself is exactly 1.000, which proves the matrix and the argmax before any claim is
read off them.

ECFP4 is Morgan with radius 2, which is what `src/feats.py` and `cypsplit.py` already use, so the
geometry here is commensurable with the split's 0.35 threshold.

Writes results/logs/k155_parents.json and results/preds/k155_parent_of_test.csv.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json
from collections import Counter

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

from submit import CYPS

TOPN = 25
PER_HIT = 10
EXPANDED = ["CYP1A2", "CYP2C9", "CYP3A4"]      # the page's three; CYP2D6 is excluded


def fps(smis, gen):
    out = []
    for s in smis:
        m = Chem.MolFromSmiles(s)
        if m is None:
            raise SystemExit("RDKit не разобрал %r" % s)
        out.append(gen.GetFingerprint(m))
    return out


def main():
    rows = pd.read_csv(D + "rows.csv")
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    FTR, FTE = fps(list(rows.SMILES), gen), fps(list(te.SMILES), gen)
    print("  обучение %d, тест %d, ECFP4 = Morgan r=2 2048" % (len(FTR), len(FTE)))

    S = np.vstack([np.asarray(DataStructs.BulkTanimotoSimilarity(f, FTR), dtype=np.float32)
                   for f in FTR])
    assert np.allclose(S.max(1), 1.0), "КОНТРОЛЬ: максимум с собой не 1"
    print("  КОНТРОЛЬ: максимум с собой = 1.000 на всех %d строках" % len(FTR))

    STE = np.vstack([np.asarray(DataStructs.BulkTanimotoSimilarity(f, FTR), dtype=np.float32)
                     for f in FTE])
    nn_idx = STE.argmax(1)
    nn_sim = STE.max(1)
    cnt = Counter(nn_idx.tolist())
    print("\n  1. КОНЦЕНТРАЦИЯ: 750 тестовых указывают на %d различных обучающих" % len(cnt))
    cov = np.cumsum(sorted(cnt.values(), reverse=True)) / len(te)
    for k in (25, 50, 75, 100, 150):
        if k <= len(cnt):
            print("     верхние %3d родителей покрывают %5.1f%% теста" % (k, 100 * cov[k - 1]))

    # the null: 750 random training compounds against their nearest OTHER training compound
    rng = np.random.default_rng(0)
    pick = rng.choice(len(rows), len(te), replace=False)
    Sx = S.copy()
    np.fill_diagonal(Sx, -1.0)
    nn_null = Sx[pick].argmax(1)
    cnull = Counter(nn_null.tolist())
    covn = np.cumsum(sorted(cnull.values(), reverse=True)) / len(te)
    print("     НУЛЬ (750 случайных обучающих): %d различных, верхние 75 покрывают %.1f%%"
          % (len(cnull), 100 * covn[min(74, len(cnull) - 1)]))

    # 2. are the parents the top-25 hits of the three expanded isoforms?
    top = {}
    for c in CYPS:
        v = tr["%s_pIC50_direct_inhibition" % c].to_numpy(float)
        order = np.argsort(np.where(np.isfinite(v), v, -np.inf))[::-1]
        top[c] = list(order[:TOPN])
    parents = [i for i, _ in cnt.most_common(75)]
    print("\n  2. РОДИТЕЛИ против топ-%d по pIC50, на каждый фермент:" % TOPN)
    hits = {}
    for c in CYPS:
        inter = len(set(parents) & set(top[c]))
        hits[c] = inter
        mark = "  <- расширялся" if c in EXPANDED else "  <- НЕ расширялся"
        print("     %-8s %2d из %d верхних родителей входят в топ-%d%s"
              % (c, inter, TOPN, TOPN, mark))
    union3 = set().union(*[set(top[c]) for c in EXPANDED])
    print("     объединение топ-25 трёх расширявшихся: %d соединений, из них родителей %d"
          % (len(union3), len(set(parents) & union3)))

    # 3. children per parent
    sizes = np.array([cnt[p] for p in parents])
    print("\n  3. ДЕТЕЙ НА РОДИТЕЛЯ среди верхних 75: медиана %.0f, среднее %.1f, диапазон %d..%d"
          % (np.median(sizes), sizes.mean(), sizes.min(), sizes.max()))
    print("     сходство ребёнка с ближайшим обучающим: медиана %.4f, 10%%=%.4f, 90%%=%.4f"
          % (np.median(nn_sim), *np.percentile(nn_sim, [10, 90])))

    out = pd.DataFrame({"Molecule_Name": te.Molecule_Name,
                        "parent_name": rows.Molecule_Name.to_numpy()[nn_idx],
                        "parent_sim": nn_sim})
    for c in CYPS:
        out["parent_%s" % c] = tr["%s_pIC50_direct_inhibition" % c].to_numpy(float)[nn_idx]
    out.to_csv(RES + "preds/k155_parent_of_test.csv", index=False)
    json.dump({"n_distinct_parents": len(cnt), "coverage_top75": float(cov[74]),
               "null_distinct": len(cnull), "null_coverage_top75": float(covn[74]),
               "parents_in_top25": hits, "children_median": float(np.median(sizes)),
               "nn_sim_median": float(np.median(nn_sim))},
              open(RES + "logs/k155_parents.json", "w"), ensure_ascii=False, indent=1)
    print("\n  записано %spreds/k155_parent_of_test.csv" % RES)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
