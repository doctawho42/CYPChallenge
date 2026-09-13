"""Item 301: estimate the blind test set's label moments WITHOUT probing the leaderboard.

Pre-registered blind in item 300, committed before this ran (commit 2dd5bc5).

Why this exists. Three competing teams now calibrate their submissions onto constants
that were reverse-engineered from leaderboard feedback -- solve R2 = 2*rho*k - k^2 - b^2
from three affine-transformed submissions of one prediction vector and the blind half's
mean and sd fall out. Those constants are public, so using them proves nothing, and they
describe the LIVE HALF only, which is split by chemical series and therefore explicitly
not guaranteed to represent the full 750.

The estimator here touches no leaderboard. It uses only the test SMILES (allowed: computing
from test structures is not de-anonymisation) and our own training labels: for each test
compound, take its nearest training neighbour by Tanimoto on the same Morgan counts the
split uses, and read that neighbour's label. The resulting distribution's mean and sd
estimate the blind moments.

Nearest-neighbour transfer is a shrinkage estimator -- it can only emit labels that already
exist, and averaging over near-duplicates compresses spread -- so its bias is calibrated
leave-one-out on the training set, where the true moments are known. The correction assumes
the bias transfers from training to test; that assumption is stated in the output rather
than hidden, because the test set is an analog expansion and its neighbour structure is
demonstrably different (adlvdl measured nearest-neighbour potency enrichment of 25.9x on
CYP3A4 against 1.6x on CYP2D6).

Runtime: about a minute. Reads committed data only; writes results/logs/k95_moments.json.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES

import json
import numpy as np
import pandas as pd

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]

# Reverse-engineered by team briford / SuperCowPowers from three of their own scored
# submissions and published in jeremycheminf/openadmet_scripts/CYP_Challenge/src/
# cyp_submission/blind_benchmark.py. Obtained BY LEADERBOARD PROBING, and describing the
# LIVE HALF (375 compounds) only -- not the full 750 this estimator covers. Kept here as
# the comparison target, never as an input.
PROBED = {
    "CYP1A2": (4.412, 1.553),
    "CYP2C9": (4.830, 1.101),
    "CYP2D6": (3.107, 1.599),
    "CYP3A4": (4.880, 1.272),
}

# Item 300's usability gate, fixed before the run.
GATE_MEAN = 0.15
GATE_SD = 0.20
# Item 300's corroboration threshold on the means, fixed before the run.
AGREE = 0.30


def tanimoto_binary(A, B, chunk=256):
    """Rows of A against rows of B: |A&B| / (|A|+|B|-|A&B|) on presence/absence.

    Binary rather than count-based MinMax because Butina clustering -- the split this
    project is pinned to -- operates on binary Morgan fingerprints at threshold 0.35, so
    "nearest neighbour" here means nearest in the same metric the folds were built in.
    """
    a = (A > 0).astype(np.float32)
    b = (B > 0).astype(np.float32)
    na = a.sum(1)
    nb = b.sum(1)
    out = np.empty((a.shape[0], b.shape[0]), dtype=np.float32)
    for i in range(0, a.shape[0], chunk):
        j = min(i + chunk, a.shape[0])
        inter = a[i:j] @ b.T
        out[i:j] = inter / np.maximum(na[i:j, None] + nb[None, :] - inter, 1e-9)
    return out


def nn_transfer(Q, Pool, labels):
    """For each query row, the label of its most similar pool row. Returns (labels, sims)."""
    S = tanimoto_binary(Q, Pool)
    k = S.argmax(1)
    return labels[k], S[np.arange(S.shape[0]), k]


def loo_transfer(Pool, labels):
    """Leave-one-out: each pool row takes the label of its nearest OTHER pool row."""
    S = tanimoto_binary(Pool, Pool)
    np.fill_diagonal(S, -1.0)
    k = S.argmax(1)
    return labels[k], S[np.arange(S.shape[0]), k]


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
    tr = tr.set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index()
    FP = np.load(D + "feats.npz")["FP"]
    TFP = np.load(D + "test_feats.npz")["FP"]
    assert FP.shape[0] == len(rows), "feats.npz is not aligned to rows.csv"
    assert TFP.shape[0] == 750, "test_feats.npz does not hold 750 rows"

    print("Оценка моментов закрытого теста БЕЗ зондирования лидерборда")
    print("Пункт 301; предрегистрация -- пункт 300, коммит 2dd5bc5\n")

    print("ШАГ 1. Калибровка смещения оценщика leave-one-out на обучающей выборке.")
    print("       (истинные моменты там известны, поэтому смещение измеримо)\n")
    print("    %-8s %6s | %7s %7s %7s | %7s %7s %7s"
          % ("фермент", "n", "ист.ср", "LOO ср", "сдвиг", "ист.sd", "LOO sd", "отнош."))
    cal, gate_ok = {}, True
    for cyp in CYPS:
        col = cyp + "_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy(float)
        est, sim = loo_transfer(FP[m], y)
        d_mean = float(est.mean() - y.mean())
        r_sd = float(est.std() / y.std())
        cal[cyp] = {"d_mean": d_mean, "r_sd": r_sd, "sim_median": float(np.median(sim))}
        ok = abs(d_mean) <= GATE_MEAN and abs(est.std() - y.std()) <= GATE_SD
        gate_ok &= ok
        print("    %-8s %6d | %7.3f %7.3f %+7.3f | %7.3f %7.3f %7.3f %s"
              % (cyp, m.sum(), y.mean(), est.mean(), d_mean, y.std(), est.std(), r_sd,
                 "" if ok else "<- вне ворот"))

    print("\n    ворота пункта 300: |сдвиг среднего| <= %.2f и |разность sd| <= %.2f"
          % (GATE_MEAN, GATE_SD))
    print("    вердикт по воротам: %s" % ("ГОДЕН" if gate_ok else "НЕ ГОДЕН"))

    print("\n\nШАГ 2. Тот же оценщик на 750 тестовых молекулах.\n")
    print("    %-8s | %7s %7s | %9s %9s | %7s %7s"
          % ("фермент", "сырое", "сырой", "с поправ.", "с поправ.", "зонд", "зонд"))
    print("    %-8s | %7s %7s | %9s %9s | %7s %7s"
          % ("", "ср", "sd", "ср", "sd", "ср", "sd"))
    out = {"gate_ok": bool(gate_ok), "calibration": cal, "estimates": {}}
    for cyp in CYPS:
        col = cyp + "_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy(float)
        est, sim = nn_transfer(TFP, FP[m], y)
        raw_mean, raw_sd = float(est.mean()), float(est.std())
        # additive correction on the mean, multiplicative on the sd, both from LOO
        cor_mean = raw_mean - cal[cyp]["d_mean"]
        cor_sd = raw_sd / cal[cyp]["r_sd"]
        pm, ps = PROBED[cyp]
        out["estimates"][cyp] = {
            "raw_mean": raw_mean, "raw_sd": raw_sd,
            "cor_mean": cor_mean, "cor_sd": cor_sd,
            "probed_mean": pm, "probed_sd": ps,
            "test_sim_median": float(np.median(sim)),
            "train_label_mean": float(y.mean()), "train_label_sd": float(y.std()),
        }
        print("    %-8s | %7.3f %7.3f | %9.3f %9.3f | %7.3f %7.3f"
              % (cyp, raw_mean, raw_sd, cor_mean, cor_sd, pm, ps))

    print("\n    медиана сходства ближайшего соседа (тест -> обучение):")
    for cyp in CYPS:
        e = out["estimates"][cyp]
        print("       %-8s тест %.3f   обучение LOO %.3f"
              % (cyp, e["test_sim_median"], cal[cyp]["sim_median"]))

    print("\n\nШАГ 3. Проверка предсказаний пункта 300, записанных до прогона.\n")
    sd_below = [out["estimates"][c]["raw_sd"] < PROBED[c][1] for c in CYPS]
    print("    (1) каждый оценённый sd НИЖЕ зондированного: %d из 4  -- %s"
          % (sum(sd_below), "подтверждено" if all(sd_below) else "НЕ подтверждено"))
    d2 = out["estimates"]["CYP2D6"]
    dirn = d2["cor_mean"] < d2["train_label_mean"]
    print("    (2) CYP2D6: оценка ниже нашего обучающего среднего (%.3f): %s -- %s"
          % (d2["train_label_mean"], "%.3f" % d2["cor_mean"],
             "подтверждено" if dirn else "НЕ подтверждено"))
    near = {c: abs(out["estimates"][c]["cor_mean"] - PROBED[c][0]) for c in CYPS}
    n_agree = sum(v <= AGREE for v in near.values())
    print("    (3) средние в пределах %.2f от зондированных: %d из 4  (%s)"
          % (AGREE, n_agree, ", ".join("%s %.3f" % (c, near[c]) for c in CYPS)))

    print("\n    ВЕРДИКТ по асимметрии, зафиксированной в пункте 300:")
    if not gate_ok:
        print("       оценщик не прошёл ворота -- числа выше НЕ являются проверкой чужих")
        print("       констант и не цитируются как таковая.")
    elif n_agree == 4:
        print("       зондированные константы ПОДТВЕРЖДЕНЫ источником, не касавшимся борда.")
    else:
        print("       расхождение НЕОДНОЗНАЧНО и НЕ опровергает константы briford:")
        print("       они описывают живую половину (375), а этот оценщик -- все 750, и")
        print("       разбиение по химическим сериям представительность не гарантирует.")

    robustness(tr, FP, TFP, out)

    with open(RES + "logs/k95_moments.json", "w") as f:
        json.dump(out, f, indent=1, ensure_ascii=False)
    print("\n    записано: %slogs/k95_moments.json" % RES)


def tanimoto_minmax(A, B, chunk=128):
    """Generalised (MinMax) Tanimoto on COUNT vectors: sum(min) / sum(max)."""
    out = np.empty((A.shape[0], B.shape[0]), dtype=np.float32)
    for i in range(0, A.shape[0], chunk):
        j = min(i + chunk, A.shape[0])
        a = A[i:j][:, None, :]
        mn = np.minimum(a, B[None, :, :]).sum(-1)
        mx = np.maximum(a, B[None, :, :]).sum(-1)
        out[i:j] = mn / np.maximum(mx, 1e-9)
    return out


def robustness(tr, FP, TFP, out):
    """POST HOC, not pre-registered: does the disagreement survive the arbitrary choices?

    Item 300 fixed k=1 and Tanimoto on the Morgan counts the split uses, which this file
    reads as binary (Butina clusters on binary fingerprints). Two choices were therefore
    mine rather than the pre-registration's: binary versus count MinMax, and k=1 versus
    averaging several neighbours. The second matters most, because averaging neighbours
    compresses sd directly and sd is one of the two quantities in dispute. These arms are
    labelled post hoc and must not be quoted as the pre-registered result.
    """
    print("\n\nШАГ 4 (ПОСТ ФАКТУМ, не предрегистрировано). Устойчивость к моим произвольным")
    print("        выборам: бинарный Танимото против счётного MinMax, и k = 1, 3, 5.\n")

    print("    4a. Смещение САМОГО оценщика при каждом k, LOO на обучении, где истина")
    print("        известна. Это отличает шринкаж от информации: если усреднение соседей")
    print("        уводит оценку вниз и здесь, то уход на тесте -- артефакт, а не сигнал.\n")
    print("    %-8s | %-23s | %-23s" % ("фермент", "сдвиг среднего (LOO)", "отношение sd (LOO)"))
    print("    %-8s | %7s %7s %7s | %7s %7s %7s"
          % ("", "k=1", "k=3", "k=5", "k=1", "k=3", "k=5"))
    loo = {}
    for cyp in CYPS:
        col = cyp + "_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy(float)
        S = tanimoto_binary(FP[m], FP[m])
        np.fill_diagonal(S, -1.0)
        order = np.argsort(-S, axis=1)
        rec = {}
        for k in (1, 3, 5):
            est = y[order[:, :k]].mean(1)
            rec[k] = {"d_mean": float(est.mean() - y.mean()),
                      "r_sd": float(est.std() / y.std())}
        loo[cyp] = rec
        print("    %-8s | %s | %s"
              % (cyp,
                 " ".join("%+7.3f" % rec[k]["d_mean"] for k in (1, 3, 5)),
                 " ".join("%7.3f" % rec[k]["r_sd"] for k in (1, 3, 5))))
    out["loo_by_k"] = {c: {str(k): v for k, v in loo[c].items()} for c in loo}

    print("\n    4b. Те же руки на тесте. СЫРЫЕ, без поправки на смещение -- поправка")
    print("        конфаундила бы именно то, что здесь проверяется.\n")
    print("    %-8s | %-21s | %-21s" % ("фермент", "бинарный: ср (sd)", "MinMax: ср (sd)"))
    print("    %-8s | %6s %6s %6s | %6s %6s %6s"
          % ("", "k=1", "k=3", "k=5", "k=1", "k=3", "k=5"))
    rob = {}
    for cyp in CYPS:
        col = cyp + "_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy(float)
        line, rec = [], {}
        for name, S in (("bin", tanimoto_binary(TFP, FP[m])),
                        ("minmax", tanimoto_minmax(TFP, FP[m]))):
            order = np.argsort(-S, axis=1)
            for k in (1, 3, 5):
                est = y[order[:, :k]].mean(1)
                rec["%s_k%d" % (name, k)] = {"mean": float(est.mean()), "sd": float(est.std())}
                line.append("%6.3f" % est.mean())
        rob[cyp] = rec
        print("    %-8s | %s | %s" % (cyp, " ".join(line[:3]), " ".join(line[3:])))
    print("\n    те же руки по sd:")
    print("    %-8s | %6s %6s %6s | %6s %6s %6s"
          % ("", "k=1", "k=3", "k=5", "k=1", "k=3", "k=5"))
    for cyp in CYPS:
        r = rob[cyp]
        print("    %-8s | %s" % (cyp, " ".join(
            "%6.3f" % r["%s_k%d" % (n, k)]["sd"] for n in ("bin", "minmax") for k in (1, 3, 5))))
    out["robustness_post_hoc"] = rob

    print("\n\nШАГ 5. Арифметика, а не вердикт: что зондированные константы требуют от")
    print("        ЗАКРЫТОЙ половины, если наша оценка полных 750 верна.\n")
    print("    Живая половина + закрытая половина = все 750, поэтому")
    print("    среднее закрытой = 2*(наша оценка 750) - (зондированная живая).\n")
    print("    %-8s %11s %11s %13s" % ("фермент", "наша 750", "зонд живой", "тогда закрытая"))
    imp = {}
    for cyp in CYPS:
        m750 = out["estimates"][cyp]["cor_mean"]
        live = PROBED[cyp][0]
        blind = 2.0 * m750 - live
        imp[cyp] = blind
        print("    %-8s %11.3f %11.3f %13.3f" % (cyp, m750, live, blind))
    out["implied_blinded_half_mean"] = imp
    print("\n    Насколько это требование экстремально -- ПОСЧИТАНО, а не заявлено:")
    print("    %-8s %13s %11s %s" % ("фермент", "тогда закрытая", "доля обуч.", "макс. обуч."))
    print("    %-8s %13s %11s %s" % ("", "среднее", "ниже неё", "метка"))
    for cyp in CYPS:
        col = cyp + "_pIC50_direct_inhibition"
        y = tr[col].dropna().to_numpy(float)
        frac = float((y < imp[cyp]).mean())
        out["estimates"][cyp]["implied_blind_frac_below"] = frac
        print("    %-8s %13.3f %11.4f %12.3f" % (cyp, imp[cyp], frac, y.max()))
    print("\n    Это условное следствие, не измерение: наша оценка 750 сама оценка.")
    print("    Асимметрия пункта 300 остаётся в силе -- это натяжение, не опровержение.")


if __name__ == "__main__":
    main()
