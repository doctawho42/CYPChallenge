"""Does the disagreement between the two measurement layers predict where the model will be wrong.

Why this is the one gap on the award's weakest axis. The organisers name three: architectural
novelty, creative data usage, novel uncertainty quantification. The first two this file has
material for. On the third it has a **finding about the benchmark's own UQ** -- the shipped credible
interval is 3.92 sigma and sigma is a function of the label at R-squared 0.93 to 0.98 (item 114), so
the interval is a deterministic function of the answer rather than a per-compound uncertainty, which
makes ST-RAE a potency-weighted absolute error (114) with a U-shaped penalty (115) and means any
team calibrating "uncertainty" against these bands is calibrating against the label.

What it does not have is a per-compound uncertainty that is **not** a function of the answer. There
is exactly one candidate: the disagreement between the curve and the single screening point under
the explicit instrument equation. Item 170 measured that its spread exceeds its own propagated error
by 2.9 to 5.8, so the quantity is real.

**Item 179 closed it as a training weight and that does not close this.** "Does it improve the fit"
and "does it predict where the model errs" are different questions and only the first was asked. The
first was answered decisively -- weighting by it loses to the base and to its own shuffle -- and
item 179's reading was that the disagreeing compounds carry signal rather than noise. That reading
cuts both ways here, which is why it needs measuring: if those compounds are informative they may
also be the ones the model gets right, and then the residual will not track the error at all.

**The confounder, named before the run because it has already burned this file once.** Both
residuals contain the label. The Hill residual is `g(y) - log2fc` and the model residual is `y - p`,
so a compound with an extreme label can have both large by construction -- exactly the mechanism
that made the microsomal-binding correlation look real in item 171 until potency was removed
isotonically and the sign reversed. So every correlation here is reported raw **and** after removing
the label, with `increasing="auto"`, since the default direction silently produced a near-constant
fit in item 171 and left the correlation numerically unchanged.

    положительна после контроля  ->  найдена посоединённая неопределённость, укоренённая в
                                     НЕЗАВИСИМОМ измерении, а не в метке
    около нуля                   ->  ноль пункта 179 обобщается с весов на неопределённость,
                                     и линия закрыта окончательно

Reads results/preds/oof.json and the screening table. Writes nothing. ~2 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES

import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.isotonic import IsotonicRegression

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
PC0 = 4.305
FIT_E = {"CYP1A2": 0.728, "CYP2C9": 0.621, "CYP2D6": 0.867, "CYP3A4": 0.931}
FIT_H = {"CYP1A2": 1.261, "CYP2C9": 1.112, "CYP2D6": 1.243, "CYP3A4": 1.968}


def resid_out(x, y):
    """Убрать x из y изотоникой. increasing='auto' обязателен --- умолчание True в пункте
    171 подогнало почти константу к убывающей связи и оставило корреляцию неизменной."""
    iso = IsotonicRegression(increasing="auto", out_of_bounds="clip").fit(x, y)
    return y - iso.predict(x)


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    P = json.load(open(RES + "preds/oof.json"))

    print(f"{'фермент':8s} {'n':>5s} | {'сырое':>8s} {'после снятия метки':>19s} | "
          f"{'sd |ост.Хилла|':>15s} {'sd |ост.модели|':>16s}")
    keep = {}
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy()
        sub = sc[sc.enzyme == c].set_index("Molecule_Name")["log2fc_estimate"]
        l2 = sub.reindex(rows.Molecule_Name).to_numpy(float)[m]
        p = np.asarray(P[f"FP+DESC+MECH|{c}"], float)
        I = FIT_E[c] / (1 + 10 ** np.clip(FIT_H[c] * (PC0 - y), -30, 30))
        rH = np.log2(np.clip(1 - I, 1e-6, None)) - l2
        ok = np.isfinite(rH)
        aH, aM, yy = np.abs(rH[ok]), np.abs(y[ok] - p[ok]), y[ok]
        raw = spearmanr(aH, aM).statistic
        ctl = spearmanr(resid_out(yy, aH), resid_out(yy, aM)).statistic
        print(f"{c:8s} {ok.sum():5d} | {raw:+8.3f} {ctl:+19.3f} | "
              f"{aH.std():15.3f} {aM.std():16.3f}")
        keep[c] = (aH, aM, yy)

    print("\n=== по децилям |остатка Хилла|: растёт ли ошибка модели ===")
    print(f"{'фермент':8s} " + "".join(f"{'д'+str(i+1):>7s}" for i in range(5))
          + "   отношение д5/д1")
    for c in CYPS:
        aH, aM, _ = keep[c]
        q = np.quantile(aH, np.linspace(0, 1, 6))
        med = []
        for i in range(5):
            sel = (aH >= q[i]) & (aH <= q[i + 1])
            med.append(np.median(aM[sel]) if sel.sum() > 5 else np.nan)
        print(f"{c:8s} " + "".join(f"{v:7.3f}" for v in med)
              + f"   {med[4]/max(med[0],1e-9):14.2f}")

    print("""
Как читать. Решает столбец «после снятия метки», а не сырой: оба остатка содержат метку по
построению, и сырая корреляция может быть целиком потенцией --- ровно так пункт 171 выглядел
подтверждением, пока потенцию не убрали и знак не перевернулся.

Положительная контролируемая корреляция означает, что расхождение двух анализов помечает
соединения, на которых модель ошибётся, --- то есть посоединённую неопределённость,
укоренённую в НЕЗАВИСИМОМ измерении, а не в ответе. Именно этого у отгруженной полосы нет:
она есть 3.92*sigma, а sigma есть функция метки при R-квадрат 0.93..0.98.

Децили внизу --- та же величина без предположения о монотонности: если ошибка модели растёт
от первого дециля к пятому, связь есть, какой бы ни была её форма.""")


if __name__ == "__main__":
    main()
