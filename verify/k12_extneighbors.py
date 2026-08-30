"""How close the external compounds get to ours, and what the zero overlap with the test means.

Two things are checked and they are not the same thing.

The first is leakage. Exact overlap with the blinded test is zero of 750, but exact matching
misses the case that matters more for a competition: a compound that is not identical to a
test compound but is one methyl away from it. So this file measures, for every test compound,
the highest Tanimoto similarity to anything in the external set, and reports the tail rather
than the average - the average is uninformative, the fraction above 0.7 is not.

The second is what the zero itself implies, and this is where the arithmetic misleads. Overlap
with our training set is 64 of 4905, or 1.3%, which on 750 test compounds predicts about ten
matches; the probability of drawing none is about five in a hundred thousand. So the zero is
not chance. Two explanations survive, and they lead to opposite conclusions:

  A. the test was assembled to avoid publicly known compounds. Then the test is depleted of
     public chemistry, and the model this document rests on - the test is anchors plus their
     analogues, drawn from the same pool as the training set - is weakened;
  B. the released external file was deduplicated against the blinded test before publication.
     Then the zero says nothing whatever about the test's composition, and nothing changes.

B is the more ordinary act: an organiser publishing a baseline alongside a blind challenge
removes their own test compounds from the training file they release. A cannot even be
executed cleanly, since "publicly known" has no boundary.

The measurement separates them. Under A the test compounds must be FARTHER from public
chemistry than our training compounds are, since that is the property A selects on. Under B
there is no reason for any difference. So the two similarity distributions - test to external,
training to external - decide it, and the direction of the difference is the whole answer.

Note the one thing this cannot rule out: under B the deduplication was by exact structure, so
near-duplicates of test compounds would have survived it. The tail of the test-to-external
distribution is therefore also the leakage check, and it has to be read before anything else.

Reads data/rows.csv, the blinded test, the external CSV. Prints only. Takes about two minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D
from cypsplit import cluster_ids

import argparse

import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator

DRAWS = 2000


def fps(smiles):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    out = []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        if m is not None:
            out.append(gen.GetFingerprint(m))
    return out


def nearest(query, ref):
    """Highest similarity from each query molecule to anything in ref."""
    return np.array([max(DataStructs.BulkTanimotoSimilarity(q, ref)) for q in query])


def describe(name, v):
    print(f"{name:26s} {np.median(v):8.3f} {np.percentile(v, 75):8.3f} "
          f"{(v > 0.7).mean():9.4f} {(v > 0.9).mean():9.4f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ext-x", required=True)
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    ext = pd.read_csv(a.ext_x)

    print("считаю отпечатки", flush=True)
    f_tr, f_te = fps(list(rows.SMILES)), fps(list(te.SMILES))
    f_ex = fps(list(ext.OPENADMET_CANONICAL_SMILES))
    print(f"  обучение {len(f_tr)}, тест {len(f_te)}, внешние {len(f_ex)}\n", flush=True)

    d_te = nearest(f_te, f_ex)
    d_tr = nearest(f_tr, f_ex)
    # Для масштаба: насколько близко наши соединения лежат друг к другу.
    d_own = np.array([sorted(DataStructs.BulkTanimotoSimilarity(q, f_tr))[-2] for q in f_tr])

    print(f"{'':26s} {'медиана':>8s} {'75%':>8s} {'доля>0.7':>9s} {'доля>0.9':>9s}")
    describe("тест -> внешние", d_te)
    describe("обучение -> внешние", d_tr)
    describe("обучение -> обучение", d_own)

    n_tr, n_te, k = len(rows), len(te), 64
    p = k / n_tr
    print(f"\nдоля пересечения с обучением {p:.4%}, ожидание на тесте {p * n_te:.1f} совпадений,"
          f"\nвероятность ровно нуля при биномиальной модели {(1 - p) ** n_te:.2e}")

    diff = np.median(d_te) - np.median(d_tr)
    print(f"\nтест ближе к внешнему набору, чем обучение, на {diff:+.3f} по медиане")

    # Интервал. Обе стороны пересчитываются, а не одна: держать одну сторону
    # фиксированной --- ровно та ошибка, что была найдена в verify/k8_kernel.py.
    # Обучение ресемплится КЛАСТЕРАМИ, потому что аналоги ходят сериями и
    # независимыми соединения там не являются.
    cid, ncl = cluster_ids(list(rows.SMILES))
    cid = np.asarray(cid)[: len(d_tr)]
    groups = [np.where(cid == g)[0] for g in np.unique(cid)]
    rng = np.random.default_rng(0)
    boot = np.empty(DRAWS)
    for i in range(DRAWS):
        pick = rng.integers(0, len(groups), len(groups))
        idx = np.concatenate([groups[j] for j in pick])
        te_i = rng.integers(0, len(d_te), len(d_te))
        boot[i] = np.median(d_te[te_i]) - np.median(d_tr[idx])
    q = np.percentile(boot, [2.5, 97.5])
    print(f"95% интервал по {DRAWS} выборкам ({ncl} кластеров): [{q[0]:+.3f}, {q[1]:+.3f}], "
          f"доля отрицательных {float((boot < 0).mean()):.3f}")

    print("\nГипотеза A-сильная --- тест собирали в обход публичных данных вообще --- требует "
          "ЗНАКА МИНУС здесь.\nГипотеза A-слабая --- тест почистили только от точных "
          "совпадений --- предсказывает ноль и от B\nнаблюдениями не отличается.")
    if q[0] > 0:
        print("\nВердикт: A-сильная опровергнута, интервал целиком положителен.")
    elif diff > 0:
        print("\nВердикт: знак против A-сильной, но интервал ноль накрывает.")
    else:
        print("\nВердикт: знак согласуется с A-сильной.")
    print("A-слабая не опровергается ничем и не должна: под ней из теста ушло бы не более "
          f"{k / n_tr:.1%}\nсоединений, против двукратного обеднения основаниями, "
          "которым объясняется сдвиг на CYP2D6.")


if __name__ == "__main__":
    main()
