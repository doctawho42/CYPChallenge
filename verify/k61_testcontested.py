"""The contested-pair rate on the blinded test set: the one reliability statistic that needs no labels.

Item 207 measured that the five members' pairwise disagreement locates the ensemble's error --
contested pairs are 35 per cent of pairs and carry 52 per cent of the discordance, a lift of
1.56 -- and that nothing in the ensemble fixes it. What survived was not the cascade but the
statistic: **whether two members order a pair the same way is computable without any label at
all**, so unlike every other diagnostic in this file it can be evaluated on the 750 blinded
test molecules directly.

That matters because of item 206. Our cross-validation is structurally blind to the test's
analog structure -- 10.2 per cent of training molecules sit in a series against 69.1 per cent of
test molecules -- so a statement about the test that does not depend on labels is worth more
here than usual.

The question this answers: **is the submission less internally consistent on the test than it is
out of fold?** If the contested rate is the same, the ensemble meets the test in the regime it
was measured in, and item 207's 57-to-62-per-cent accuracy on contested pairs transfers. If the
test rate is markedly higher, a larger share of the submitted ordering sits in the regime where
the ensemble is barely better than a coin, and the out-of-fold rank estimate is optimistic for a
reason that has nothing to do with covariate shift in the usual sense.

**Four members, not five.** The trunk scores 0.4956 to 0.5146 on contested pairs (item 207) --
a coin on every enzyme -- so including it inflates the contested rate by adding a near-random
vote rather than a disagreeing opinion. Both sides are therefore computed on the same four
sklearn members, which is what makes train and test comparable. The five-member figure would be
larger and would mean less.

Targets for the dead-zone pass are read from the member cache written by k58 rather than
recomputed, so this costs sixteen full fits rather than a second full out-of-fold pass.

Reads the k58 cache, data/feats.npz, data/rows.csv and the blinded test file. Writes nothing.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time

import numpy as np
import pandas as pd

import submit as SB

CYPS = SB.CYPS
FOUR = ["поферментно", "пул", "GP", "гребневая"]


def contested(P):
    """Доля спорных пар среди всех: члены расходятся в упорядочении пары."""
    n = P.shape[1]
    iu = np.triu_indices(n, 1)
    sm = np.stack([np.sign(p[iu[0]] - p[iu[1]]) for p in P])
    return float((~np.all(sm == sm[0], axis=0)).mean()), len(iu[0])


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True)
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)

    desc_names = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
    mech_names = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    te, Xte = SB.test_features(desc_names, mech_names)
    print(f"тест {Xte.shape}, обучение {X.shape}", flush=True)

    C = json.load(open(a.cache))
    plain = {k: [np.asarray(v, float) for v in P] for k, P in C["parts"]}
    dzp = {k: [np.asarray(v, float) for v in P] for k, P in C["dzp"]}

    T = []
    for e, c in enumerate(CYPS):
        m = mask[:, e]
        Pte = []
        for kind in FOUR:
            t0 = time.time()
            tgt = np.clip(plain[kind][e], LO[m, e], HI[m, e])
            Xtr_d, Xte_d = SB._dz_design(kind, X[m], e, Xte)
            Pte.append(SB._dz_fit(kind, Xtr_d, tgt, Xte_d))
            print(f"    {c} {kind:12s} {time.time()-t0:.0f} с", flush=True)
        f_te, n_te = contested(np.stack(Pte))
        f_tr, n_tr = contested(np.stack([dzp[k][e] for k in FOUR]))
        T.append({"фермент": c, "обуч. пар": n_tr, "спорных обуч.": f_tr,
                  "тест пар": n_te, "спорных тест": f_te,
                  "отношение": f_te / max(f_tr, 1e-12)})

    df = pd.DataFrame(T).set_index("фермент")
    print()
    print(df.round(4).to_string())
    print(f"\nмакро: обучение {df['спорных обуч.'].mean():.4f}, "
          f"тест {df['спорных тест'].mean():.4f}, "
          f"отношение {df['отношение'].mean():.3f}")
    print("""
Как читать. Отношение около единицы --- тест встречает ансамбль в том же режиме внутренней
согласованности, в котором его мерили вне фолда, и точность 57-62 % на спорных парах из
пункта 207 переносится как есть. Отношение заметно выше единицы --- на тесте бОльшая доля
подаваемого порядка лежит в области, где ансамбль едва лучше монетки, и оценка ранга вне
фолда оптимистична по причине, не сводящейся к обычному сдвигу ковариат.

Величина сама по себе ничего не исправляет: пункт 207 показал, что чинить спорные пары
внутри имеющегося ансамбля нечем. Это измерение НАДЁЖНОСТИ, а не источник поправки, и его
единственное достоинство --- что оно про тест и не требует ни одной метки.""")


if __name__ == "__main__":
    main()
