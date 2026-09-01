"""Does CYP3A4 need a two-site dose-response, and does it need it alone.

Why this and not another calibration patch. Five independent measurements point at one enzyme:
item 129 (two populations), 167 (nothing ever clears its floor), 169 (transporting the calibration
across halves costs +1.365 against +0.019 to +0.205 elsewhere, with the slope going 0.503 on the
lower half of the labels to 2.516 on the upper), 170 (the only enzyme whose median slope is not one
under the fitted amplitude), 171 (the only enzyme whose layer offset survives a change of
parameterisation). Items 173 and 175 add two more hints from the trunk's fitted offset.

**A five-fold swing of the Hill slope inside one enzyme is not a mis-set constant; it is the wrong
functional form.** A single-site sigmoid cannot produce it at any (E, h).

And CYP3A4 has the textbook reason for exactly that shape. Its cavity is the largest of the four,
it binds two ligands at once, and its dose-response is classically described by homotropic
cooperativity rather than by one site. CYP2D6's pocket is narrow with a salt bridge -- one site, one
sigmoid. So "the form breaks on CYP3A4 and nowhere else" is what the biochemistry says about these
four enzymes and what five measurements found without being told to look.

The models, both fitted to the same cells by the same least squares:

    односайтовая   I = E / (1 + 10^{h (pC0 - pi)})                       2 параметра
    двухсайтовая   I = E [ f/(1 + 10^{h(pC0-pi)}) + (1-f)/(1 + 10^{h(pC0-pi-D)}) ]   4 параметра

The second is a mixture of two sigmoids with a shared slope, a fraction `f` on the high-affinity
site and an affinity separation `D`. It reduces exactly to the first at f = 1 or D = 0, so the
comparison is nested and the extra parameters cannot make it worse in-sample.

Because it cannot be worse in-sample, the comparison is on **corrected AIC**, which charges for the
two extra parameters, and on a **five-fold cross-validated** residual, which charges for whatever
the AIC does not. Both are reported; a gain that appears on AIC and not on cross-validation is
over-fitting with a formal blessing.

Pre-registered, all three outcomes:

    лучше на 3A4 и НЕ лучше на трёх остальных -> кооперативность, найденная из данных
    лучше везде                               -> просто больше параметров, закрыто
    не лучше нигде                            -> форма ни при чём, и тогда пункт 129
                                                 (две кампании) остаётся единственным
                                                 объяснением адреса

Reads the inhibition, Emax and screening tables. Writes nothing. ~1 minute.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D

import numpy as np
import pandas as pd
from scipy.optimize import least_squares

from cypsplit import butina_folds

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
PC0 = 4.305


def pred1(t, pi):
    E, h = t
    return np.log2(np.clip(1 - E / (1 + 10 ** np.clip(h * (PC0 - pi), -30, 30)), 1e-6, None))


def pred2(t, pi):
    E, h, f, Dd = t
    f = np.clip(f, 0.0, 1.0)
    a = 1 / (1 + 10 ** np.clip(h * (PC0 - pi), -30, 30))
    b = 1 / (1 + 10 ** np.clip(h * (PC0 - pi - Dd), -30, 30))
    return np.log2(np.clip(1 - E * (f * a + (1 - f) * b), 1e-6, None))


def fit(pred, x0, lo, hi, pi, y):
    r = least_squares(lambda t: pred(t, pi) - y, x0, bounds=(lo, hi))
    return r.x, (pred(r.x, pi) - y)


def aicc(res, k):
    n = len(res)
    ll = -n / 2 * np.log(np.maximum(np.mean(res ** 2), 1e-12))
    a = 2 * k - 2 * ll
    return a + (2 * k * (k + 1)) / max(n - k - 1, 1)


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    fold, _ = butina_folds(list(rows.SMILES))

    print(f"{'фермент':8s} {'n':>5s} | {'односайт ско':>13s} {'двухсайт ско':>13s} "
          f"{'AICc разн.':>11s} | {'CV одно':>9s} {'CV два':>9s} {'CV разн.':>9s}")
    out = []
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        sub = sc[sc.enzyme == c].set_index("Molecule_Name")["log2fc_estimate"]
        pi = tr[col].to_numpy(float)
        y = sub.reindex(rows.Molecule_Name).to_numpy(float)
        m = np.isfinite(pi) & np.isfinite(y)
        pi, y, fi = pi[m], y[m], fold[m]

        t1, r1 = fit(pred1, [1.0, 1.0], [0.2, 0.2], [1.2, 4.0], pi, y)
        t2, r2 = fit(pred2, [t1[0], t1[1], 0.5, 0.5],
                     [0.2, 0.2, 0.0, -3.0], [1.2, 4.0, 1.0, 3.0], pi, y)
        d_aic = aicc(r2, 4) - aicc(r1, 2)          # отрицательное = двухсайтовая лучше

        cv1, cv2 = [], []
        for f in range(5):
            trn, te = fi != f, fi == f
            if te.sum() < 10 or trn.sum() < 50:
                continue
            a1, _ = fit(pred1, [1.0, 1.0], [0.2, 0.2], [1.2, 4.0], pi[trn], y[trn])
            a2, _ = fit(pred2, [a1[0], a1[1], 0.5, 0.5],
                        [0.2, 0.2, 0.0, -3.0], [1.2, 4.0, 1.0, 3.0], pi[trn], y[trn])
            cv1.append(np.mean((pred1(a1, pi[te]) - y[te]) ** 2))
            cv2.append(np.mean((pred2(a2, pi[te]) - y[te]) ** 2))
        s1, s2 = np.sqrt(np.mean(cv1)), np.sqrt(np.mean(cv2))
        print(f"{c:8s} {len(pi):5d} | {r1.std():13.4f} {r2.std():13.4f} {d_aic:11.1f} | "
              f"{s1:9.4f} {s2:9.4f} {s2-s1:+9.4f}")
        out.append((c, t2, s2 - s1, d_aic))

    print(f"\n{'фермент':8s} {'E':>7s} {'h':>7s} {'доля 1-го сайта':>16s} {'разделение D':>13s}")
    for c, t2, _, _ in out:
        print(f"{c:8s} {t2[0]:7.3f} {t2[1]:7.3f} {t2[2]:16.3f} {t2[3]:13.3f}")

    print("""
Как читать. Решает столбец «CV разн.», а не AICc: вложенная модель с двумя лишними параметрами
не может проиграть внутри выборки, поэтому только отложенный остаток говорит, есть ли там
что-то. Отрицательное значение --- двухсайтовая лучше.

Три исхода предрегистрированы. Выигрыш ТОЛЬКО на CYP3A4 --- кооперативность, найденная данными,
и тогда её стоит вносить в прибор. Выигрыш везде --- это лишние параметры, и вопрос закрыт.
Отсутствие выигрыша нигде --- форма ни при чём, и адрес CYP3A4 объясняется составом (пункт 129),
а не физикой.

Параметры внизу читаются только если выигрыш есть. Доля первого сайта около единицы или
разделение около нуля означают вырождение в односайтовую, каким бы ни был остаток.""")


if __name__ == "__main__":
    main()
