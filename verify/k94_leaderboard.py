"""Item 299: what the published rules and the live boards actually say.

Three things this reproduces, because item 299 quotes all of them:

  1. Per-enzyme ST-RAE against band occupancy on the SHIPPED out-of-fold arm. The point
     is not the score (item 294 already has it) but the coupling: ST-RAE is normalised by
     the spread of y_true and its error is clipped to the credible band, so a wider band
     and a wider spread both lower it. That is what makes an out-of-fold number on the
     training set incomparable to a leaderboard number on a differently-composed test set.

  2. MCC against decision threshold on the shipped TDI bundle. This was written to test a
     specific suspicion -- that the shipped positive rate (0.38 / 0.48) is wrong because
     the board's own numbers imply a test prevalence near 0.16 -- and it refutes it:
     matching prevalence LOWERS macro MCC, and the MCC-optimal rate is itself near 0.47.

  3. The board snapshot as a dated constant, since the boards move and the item's claims
     are anchored to what was visible on 13 September 2026.

Runtime: a few seconds. Reads only committed predictions; writes nothing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, tutorial
tutorial()

import json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.metrics import matthews_corrcoef
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]

# Observed on https://openadmet-cyp-challenge.hf.space, 13 September 2026, ~13:00 MSK.
# Kept as data because the boards move; the item's comparisons are anchored to this.
BOARD_REGRESSION = [
    # username, MA-ST-RAE, MA-MAE, MA-R2, MA-Spearman, MA-Kendall, proprietary, open code
    ("admetchallenger", 0.4075, 0.6101, -0.0109, 0.7747, 0.5973, True, True),
]
BOARD_TDI = [
    # username, MA-MCC, MA-Accuracy, MA-Precision, MA-Recall, MA-F1
    ("PIPL-Bio", 0.4502, 0.8620, 0.5741, 0.5138, 0.5371),
    ("Asidsal11", 0.4410, 0.8622, 0.7116, 0.4133, 0.4940),
    ("harinu123", 0.4407, 0.8578, 0.5318, 0.5380, 0.5321),
    ("dopamine", 0.4212, 0.8587, 0.6705, 0.3838, 0.4801),
    ("deoxys", 0.4183, 0.8360, 0.5146, 0.5514, 0.5258),
    ("jeremy", 0.4025, 0.8433, 0.4933, 0.5198, 0.5029),
    ("JSHCF", 0.3846, 0.8457, 0.5467, 0.4379, 0.4777),
    ("SVK091", 0.3737, 0.8447, 0.5375, 0.4066, 0.4551),
    ("EleanorRigby", 0.3708, 0.8208, 0.4220, 0.5711, 0.4721),
    ("Gashaw", 0.3694, 0.8303, 0.4702, 0.4912, 0.4627),
    ("vyshuG15", 0.3610, 0.8053, 0.4398, 0.5605, 0.4845),
    ("ppsqq", 0.3576, 0.8205, 0.4363, 0.5234, 0.4725),
]


def implied_prevalence(accuracy, precision, recall):
    """Invert accuracy from precision and recall to get the positive rate.

    With prevalence p: TP = recall*p, FP = TP*(1/precision - 1), TN = 1 - p - FP, so
    accuracy = TP + TN = 1 - p*(1 + recall/precision - 2*recall). Solve for p.

    This is an ESTIMATE and its weakness is stated wherever it is used: the board's
    precision and recall are already macro-averaged over two enzymes and rounded to four
    places, and a macro-average of ratios does not invert to a single prevalence exactly.
    """
    k = 1.0 + recall / precision - 2.0 * recall
    return (1.0 - accuracy) / k


def band_table():
    """ST-RAE and rank per enzyme, beside band width and band occupancy."""
    tr = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
    rows = pd.read_csv(D + "rows.csv")
    tr = tr.set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index()
    P = json.load(open("results/preds/oof_submitted.json"))["P"]

    print("Занятость полосы и ST-RAE на ПОДАВАЕМОМ OOF")
    print("(меньше ST-RAE лучше; полоса и разброс оба его снижают)\n")
    print("    %-7s %6s %8s %8s %9s %8s %8s"
          % ("фермент", "n", "|y-ср|", "ширина", "в полосе", "ST-RAE", "ранг"))
    srae, rank = [], []
    for e, cyp in enumerate(CYPS):
        col = cyp + "_pIC50_direct_inhibition"
        lo, hi = col + "_conf_low", col + "_conf_high"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy(float)
        l = tr.loc[m, lo].to_numpy(float)
        h = tr.loc[m, hi].to_numpy(float)
        p = np.asarray(P[e], float)
        assert len(p) == len(y), "%s: %d predictions against %d labels" % (cyp, len(p), len(y))
        s = float(strae(y, p, y_true_upper=h, y_true_lower=l))
        r = float(spearmanr(y, p).statistic)
        srae.append(s)
        rank.append(r)
        print("    %-7s %6d %8.3f %8.3f %9.3f %8.4f %8.4f"
              % (cyp, m.sum(), np.abs(y - y.mean()).mean(), (h - l).mean(),
                 np.mean((p >= l) & (p <= h)), s, r))
    print("    %-7s %6s %8s %8s %9s %8.4f %8.4f"
          % ("макро", "", "", "", "", np.mean(srae), np.mean(rank)))
    print("\n    ранг 0.6434 воспроизводит пункт 294 (0.6433) --- сверка конвейера")
    return float(np.mean(srae)), float(np.mean(rank))


def mcc_table():
    """MCC at the shipped rule, at its own optimum, and at prevalence-matched."""
    d = json.load(open("results/preds/bundle_oof.json"))
    print("\n\nMCC против порога на подаваемом TDI-пучке")
    print("(проверка подозрения, что доля положительных 0.38/0.48 --- дефект)\n")
    print("    %-7s %6s %8s %7s %9s %7s %11s %7s"
          % ("фермент", "prev", "MCC подача", "доля", "MCC опт", "доля", "MCC под prev", "доля"))
    got = {}
    for key, v in d.items():
        cyp = key.split("|")[1]
        y = np.asarray(v["y"], int)
        cal = np.asarray(v["cal"], float)
        dec = np.asarray(v["dec"], int)
        prev = float(y.mean())

        m_dec = float(matthews_corrcoef(y, dec))
        best_s, best_t = -9.0, None
        for t in np.unique(np.round(cal, 4)):
            q = (cal >= t).astype(int)
            if q.sum() == 0 or q.sum() == len(q):
                continue
            s = float(matthews_corrcoef(y, q))
            if s > best_s:
                best_s, best_t = s, float(t)
        r_best = float((cal >= best_t).mean())

        t_prev = float(np.quantile(cal, 1.0 - prev))
        q_prev = (cal >= t_prev).astype(int)
        m_prev = float(matthews_corrcoef(y, q_prev))

        got[cyp] = (m_dec, best_s, m_prev)
        print("    %-7s %6.3f %8.4f %7.3f %9.4f %7.3f %11.4f %7.3f"
              % (cyp, prev, m_dec, dec.mean(), best_s, r_best, m_prev, q_prev.mean()))

    macro = [float(np.mean([got[c][i] for c in got])) for i in range(3)]
    print("    %-7s %6s %8.4f %7s %9.4f %7s %11.4f"
          % ("макро", "", macro[0], "", macro[1], "", macro[2]))
    return macro


def shipped_rates():
    print("\n\nДоля True в подаваемом файле TDI")
    d = pd.read_csv("results/submission/tdi_submission.csv")
    for c in ["CYP2D6_is_TDI", "CYP3A4_is_TDI"]:
        print("    %-16s %.4f" % (c, d[c].mean()))


def boards(macro_srae, macro_rank, macro_mcc):
    print("\n\nБорды на 13 сентября 2026 (снимок в константах этого файла)\n")
    print("    регрессия: %d участник(ов)" % len(BOARD_REGRESSION))
    for u, s, mae, r2, rho, tau, prop, code in BOARD_REGRESSION:
        print("      %-16s ST-RAE %.4f  MAE %.4f  R2 %+.4f  rho %.4f  tau %.4f  пропр.%s  код%s"
              % (u, s, mae, r2, rho, tau, "да" if prop else "нет", "да" if code else "нет"))
    print("      %-16s ST-RAE %.4f                          rho %.4f   <- наш OOF"
          % ("(мы, OOF)", macro_srae, macro_rank))

    print("\n    TDI: %d участник(ов), MA-MCC %.4f ... %.4f"
          % (len(BOARD_TDI), BOARD_TDI[-1][1], BOARD_TDI[0][1]))
    print("      %-16s MCC %.4f   <- наш OOF (пучок, пункт 250)" % ("(мы, OOF)", macro_mcc))

    ps = [implied_prevalence(a, p, r) for _, _, a, p, r, _ in BOARD_TDI]
    print("\n    prevalence теста, выведенная из accuracy/precision/recall каждой строки:")
    print("      медиана %.4f, разброс %.4f-%.4f по %d строкам"
          % (float(np.median(ps)), float(min(ps)), float(max(ps)), len(ps)))
    print("      ОЦЕНКА, не измерение: precision и recall на борде уже усреднены по двум")
    print("      ферментам и округлены, а макро-среднее отношений не обращается точно.")


if __name__ == "__main__":
    ms, mr = band_table()
    mcc = mcc_table()
    shipped_rates()
    boards(ms, mr, mcc[0])
