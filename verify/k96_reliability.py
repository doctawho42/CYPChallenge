"""A per-compound reliability artefact, built only from quantities that need no labels.

What this is for. Item 207 found that member disagreement locates the ensemble's pairwise error
-- contested pairs are 35 per cent of pairs and carry 52 per cent of the discordance, a lift of
1.56 -- and that nothing in the ensemble repairs it: majority vote loses to the mean on four
enzymes of four, and item 217 confirmed that on four seeds with the mean beaten in one cell of
sixteen. So the cascade is closed. What item 207 said survives is the *diagnostic*, and item 208
showed it transfers: the contested rate is computable on the 750 blinded molecules with no labels
at all, and its macro ratio test-to-train is 1.012.

This file turns that aggregate into a **per-compound** quantity, and adds a second label-free one
beside it: how far each compound sits from the training set it was predicted from. Together they
answer "where should this submission not be trusted", using only the test set itself.

**It is not a correction and does not move the score.** Item 207 measured that the pairs where
members argue are hard rather than mis-ordered -- the mean already scores 0.58 there, well above a
coin -- so repairing them "would have to be right where five diverse models jointly are not, which
is a demand for new information, not a rearrangement of what is present". Anything here that looked
like a score gain would be that mistake in a new wrapper.

**Four members, not five.** The trunk scores 0.4956 to 0.5146 on contested pairs (item 207) -- a
coin on every enzyme -- so it inflates the rate with a random vote rather than a disagreeing
opinion. Item 208 excluded it for the same reason, and both sides must be computed the same way for
the numbers to be comparable.

**The control that matters.** The aggregate train-side rate computed here must reproduce item 208's
column (0.2820 / 0.2683 / 0.3139 / 0.2319). If it does not, this file is computing a different
quantity under the same name, and nothing below it means anything. The check prints first and says
so either way.

Scope. The train side is arithmetic over the members cached by k84 (`members_seed{s}.json`,
dead-zone-passed) and costs nothing. The **test-side contested share is deliberately not here**: it
needs the four members refitted on the test rows against `clip(plain, LO, HI)`, and the cache holds
only the passed predictions, not the plain ones. That half waits for a run of its own. What the
test side gets here is the neighbour statistic, which needs only fingerprints.

Reads committed predictions and features; writes results/logs/k96_reliability.json.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES

import json

import numpy as np
import pandas as pd

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
FOUR = ["поферментно", "пул", "GP", "гребневая"]
# Item 208, column "спорных обуч." -- the control this file must reproduce.
ITEM_208_TRAIN = {"CYP1A2": 0.2820, "CYP2C9": 0.2683, "CYP2D6": 0.3139, "CYP3A4": 0.2319}
TOL = 0.005


def contested_matrix(P):
    """Per-compound contested share and the aggregate, from members x compounds.

    A pair is contested when the members do not agree on the sign of the difference, which is
    exactly k61's definition. Done as a full n x n sign comparison rather than over the upper
    triangle, because the per-compound share needs every pair an individual sits in.
    """
    n = P.shape[1]
    sgn = np.sign(P[:, :, None] - P[:, None, :]).astype(np.int8)
    agree = np.all(sgn == sgn[0], axis=0)
    contested = ~agree
    np.fill_diagonal(contested, False)
    per_compound = contested.sum(1) / (n - 1)
    iu = np.triu_indices(n, 1)
    return per_compound, float(contested[iu].mean()), len(iu[0])


def tanimoto_max(Q, Pool, chunk=256):
    """For each query row, the largest binary Tanimoto against the pool, and its index."""
    a = (Q > 0).astype(np.float32)
    b = (Pool > 0).astype(np.float32)
    na, nb = a.sum(1), b.sum(1)
    best = np.zeros(a.shape[0], np.float32)
    who = np.zeros(a.shape[0], np.int64)
    for i in range(0, a.shape[0], chunk):
        j = min(i + chunk, a.shape[0])
        inter = a[i:j] @ b.T
        sim = inter / np.maximum(na[i:j, None] + nb[None, :] - inter, 1e-9)
        who[i:j] = sim.argmax(1)
        best[i:j] = sim.max(1)
    return best, who


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    FP = np.load(D + "feats.npz")["FP"]
    TFP = np.load(D + "test_feats.npz")["FP"]
    mem = json.load(open(RES + "preds/members_seed0.json"))

    print("Поштучная надёжность: величины, не требующие меток")
    print("Пункт 207 -- диагностика выживает, каскад нет. Пункт 208 -- она переносится на тест.\n")

    print("КОНТРОЛЬ. Агрегатная доля спорных пар на обучении против пункта 208.")
    print("Если строка не сходится, ниже считается не та величина.\n")
    print("    %-8s %8s %10s %10s %9s %s"
          % ("фермент", "n", "здесь", "пункт 208", "разн.", "вердикт"))
    per_enz, ok_all = {}, True
    for e, cyp in enumerate(CYPS):
        P = np.stack([np.asarray(mem[k][e], float) for k in FOUR])
        share, agg, npairs = contested_matrix(P)
        d = agg - ITEM_208_TRAIN[cyp]
        ok = abs(d) <= TOL
        ok_all &= ok
        per_enz[cyp] = {"per_compound": share, "aggregate": agg, "n_pairs": npairs}
        print("    %-8s %8d %10.4f %10.4f %+9.4f %s"
              % (cyp, P.shape[1], agg, ITEM_208_TRAIN[cyp], d, "сходится" if ok else "РАСХОДИТСЯ"))
    print("\n    вердикт контроля: %s (допуск %.3f)"
          % ("ПРОЙДЕН" if ok_all else "НЕ ПРОЙДЕН -- дальше не читать", TOL))

    print("\n\nРаспределение поштучной доли спорных пар на обучении.\n")
    print("    %-8s %8s %8s %8s %8s %8s"
          % ("фермент", "мин", "квант25", "медиана", "квант75", "макс"))
    for cyp in CYPS:
        s = per_enz[cyp]["per_compound"]
        print("    %-8s %8.3f %8.3f %8.3f %8.3f %8.3f"
              % (cyp, s.min(), np.quantile(s, .25), np.median(s), np.quantile(s, .75), s.max()))

    print("\n\nБлизость тестовых молекул к обучающей выборке (по 750, поферментно).")
    print("Пул ограничен строками, размеченными для этого фермента: сосед без метки\n"
          "ничего не говорит о том, на чём модель училась предсказывать этот фермент.\n")
    print("    %-8s %8s %9s %9s %9s %9s"
          % ("фермент", "пул", "медиана", "квант25", "доля<0.4", "доля<0.3"))
    test_out = {}
    for e, cyp in enumerate(CYPS):
        col = cyp + "_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        sim, who = tanimoto_max(TFP, FP[m])
        nb_label = tr.loc[m, col].to_numpy(float)[who]
        test_out[cyp] = {"sim": sim.tolist(), "neighbour_label": nb_label.tolist()}
        print("    %-8s %8d %9.3f %9.3f %9.3f %9.3f"
              % (cyp, m.sum(), float(np.median(sim)), float(np.quantile(sim, .25)),
                 float((sim < 0.4).mean()), float((sim < 0.3).mean())))

    out = {
        "что": "поштучная надёжность: спорные пары (обучение) и близость к обучению (тест)",
        "контроль_208": {"пройден": bool(ok_all), "допуск": TOL,
                         "здесь": {c: per_enz[c]["aggregate"] for c in CYPS},
                         "пункт_208": ITEM_208_TRAIN},
        "обучение_спорные_поштучно": {c: per_enz[c]["per_compound"].tolist() for c in CYPS},
        "тест_близость": test_out,
        "молекулы_теста": te.Molecule_Name.tolist(),
    }
    with open(RES + "logs/k96_reliability.json", "w") as f:
        json.dump(out, f)
    print("\n    записано: %slogs/k96_reliability.json" % RES)
    print("""
Как читать. Обе величины измеряют НЕНАДЁЖНОСТЬ, а не ошибку, и ни одна не двигает счёт.
Высокая доля спорных пар означает, что члены расходятся в упорядочении этого соединения
относительно остальных; низкое сходство с обучением означает, что предсказание опирается
на дальних соседей. Пункт 207 измерил, что чинить первое внутри имеющегося ансамбля нечем,
поэтому это диагностика, а не источник поправки.

Тестовая доля спорных пар здесь НЕ считается: она требует доучивания четырёх членов на
тестовых строках против clip(plain, lo, hi), а кэш держит только предсказания после прохода.""")


if __name__ == "__main__":
    main()
