"""When does band-aware training help, on data where the answer is known by construction.

Why synthetic. The dead-zone objective is the largest reproducible effect in this file --- +0.033 of
rank as a single model over four seeds (item 148), +0.0167 in the submitted five-member ensemble
(item 164) --- and item 148 also found that its mechanism is the **opposite** of the one it was
proposed with: the gain comes from **narrow** bands, not wide ones, because a wide band means the
error is already free and the majorised target equals the current prediction.

That is a claim about a general property of tolerance-band losses, and it cannot be settled on four
enzymes. Here the band structure is generated rather than observed, so the quantity that is supposed
to drive the effect can be swept directly.

**The controlled variable is band width relative to model error**, not band width. That follows from
item 148's mechanism: what matters is whether the prediction lands outside its own band. Bands are
therefore drawn as `w_i = kappa * sigma * u_i`, with `sigma` the achievable residual scale and `u_i`
a positive draw, and `kappa` swept.

    kappa -> 0     полоса вырождается в точку, мёртвая зона совпадает с L1, выигрыш обязан исчезнуть
    kappa -> inf   всё внутри полосы, ошибка у всех бесплатна, выигрыш обязан исчезнуть тоже
    середина       если механизм пункта 148 верен, выигрыш имеет МАКСИМУМ где-то здесь

A single-peaked curve is the pre-registered prediction and a monotone one refutes the mechanism.
Two further conditions are swept because both occur in the real data and are confounded there:
whether band width correlates with the label (in our assay it does --- the band is 3.92 sigma and
sigma tracks the label at R-squared 0.93 to 0.98), and how heavy the width distribution is.

The learner and the majorise-minimise step are the ones the repository uses, so the synthetic result
is about the same procedure and not a stand-in for it: one pass, target `clip(p, lo, hi)`, refit
under L1, predictions taken out of fold.

Writes nothing. ~10 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, tutorial
tutorial()

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.tree import DecisionTreeRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

NTREE, LR, DEPTH, MF = 120, 0.08, 5, 0.5
N, P, FOLDS = 1500, 60, 5
KAPPAS = [0.0, 0.25, 0.5, 1.0, 2.0, 4.0, 8.0]


def boost(Xtr, ytr, Xte, seed, l1=False):
    rng = np.random.default_rng(seed)
    s = np.full(len(ytr), float(np.median(ytr) if l1 else ytr.mean()))
    pred = np.full(len(Xte), s[0])
    for _ in range(NTREE):
        r = np.sign(ytr - s) if l1 else (ytr - s)
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MF,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(Xtr, r)
        s = s + LR * t.predict(Xtr)
        pred = pred + LR * t.predict(Xte)
    return pred


def cv(X, y, fold, seed, l1=False, target=None):
    tgt = y if target is None else target
    p = np.zeros(len(y))
    for f in range(FOLDS):
        trn, te = fold != f, fold == f
        p[te] = boost(X[trn], tgt[trn], X[te], seed * 10 + f, l1)
    return p


def make(seed, corr_with_y, heavy):
    """Данные с известной достижимой ошибкой и управляемой шириной полосы."""
    rng = np.random.default_rng(seed)
    X = rng.normal(size=(N, P))
    w = rng.normal(size=P) * (np.arange(P) < 12)          # только 12 колонок несут сигнал
    f = X @ w + 1.5 * np.tanh(X[:, 0] * X[:, 1]) + 0.8 * X[:, 2] ** 2
    f = (f - f.mean()) / f.std()
    y = f + rng.normal(scale=0.6, size=N)                  # шум = достижимый предел
    u = rng.lognormal(0, 0.9, N) if heavy else np.abs(rng.normal(1.0, 0.25, N))
    if corr_with_y:
        # Ширина растёт со значением метки --- как в анализе, где полоса есть функция метки.
        u = u * (0.5 + 1.0 * (y - y.min()) / (y.max() - y.min() + 1e-9))
    return X, y, u / u.mean()


def main():
    rows = []
    for seed in range(3):
        for corr in (False, True):
            for heavy in (False, True):
                X, y, u = make(seed, corr, heavy)
                fold = np.arange(N) % FOLDS
                base = cv(X, y, fold, seed)                    # L2, один раз на набор
                l1p = cv(X, y, fold, seed, l1=True)
                sig = float(np.std(y - base))
                for k in KAPPAS:
                    half = 0.5 * k * sig * u
                    lo, hi = y - half, y + half
                    dz = cv(X, y, fold, seed, l1=True, target=np.clip(base, lo, hi))
                    row = dict(seed=seed, corr=corr, heavy=heavy, kappa=k,
                               доля_вне=float(np.mean((base < lo) | (base > hi))))
                    for nm, p in (("L2", base), ("L1", l1p), ("МЗ", dz)):
                        row[f"{nm} пара"] = float(strae(y, p, y_true_upper=hi, y_true_lower=lo))
                        row[f"{nm} ранг"] = float(spearmanr(y, p).statistic)
                    rows.append(row)
        print(f"  сид {seed} готов", flush=True)

    df = pd.DataFrame(rows)
    g = df.groupby("kappa").agg({"доля_вне": "mean", "L2 ранг": "mean", "L1 ранг": "mean",
                                 "МЗ ранг": "mean", "L2 пара": "mean", "МЗ пара": "mean"})
    g["выигрыш МЗ, ранг"] = g["МЗ ранг"] - g["L2 ранг"]
    g["выигрыш МЗ, пара"] = g["МЗ пара"] - g["L2 пара"]
    print("\n=== по ширине полосы относительно ошибки модели ===")
    print(g.round(4).to_string())

    print("\n=== разложение по двум условиям (kappa усреднено по 0.25..2.0) ===")
    m = df[df.kappa.between(0.25, 2.0)].copy()
    m["выигрыш"] = m["МЗ ранг"] - m["L2 ранг"]
    print(m.groupby(["corr", "heavy"])["выигрыш"].agg(["mean", "std", "size"]).round(4).to_string())

    print("""
Как читать. Предрегистрировано: выигрыш мёртвой зоны обязан быть ОДНОВЕРШИННЫМ по kappa ---
исчезать и при вырожденной полосе (там она совпадает с L1), и при очень широкой (там ошибка у
всех уже бесплатна, проекция ничего не меняет). Монотонная кривая опровергает механизм пункта
148, по которому выигрыш берётся с УЗКИХ полос.

Столбец «доля_вне» --- сколько предсказаний лежит за пределами своей полосы. По механизму
именно он, а не сама ширина, должен управлять выигрышем, и максимум выигрыша должен стоять
там, где эта доля заметна, но не равна единице.

Разложение внизу отвечает, зависит ли эффект от того, что ширина полосы связана с меткой ---
в нашем анализе связана, и на реальных данных эти два условия разделить нельзя.""")


if __name__ == "__main__":
    main()
