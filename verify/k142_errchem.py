"""What KIND of molecule supplies the error? The one stratification of the ST-RAE numerator nobody
has run.

WHAT IS ALREADY KNOWN, so this does not re-derive it. Item 115 stratified the numerator by POTENCY
and found it U-shaped: 45.5 per cent from the potent quartile, 27.6 from the weak, 26.9 from the two
middle ones. Item 362 stratified by potency AND label noise and found our deficit largest where the
labels are cleanest -- on the strong half of CYP1A2 we take 30 per cent of the achievable rank
against 49 on the weak half. Both cut the data by the LABEL. Neither asks what the molecules are.

WHY THE OBVIOUS VERSION OF THIS QUESTION IS WORTHLESS, and the design exists to avoid it. Comparing
high-error against low-error compounds re-discovers item 115: potent compounds carry more error, and
potency correlates with size and lipophilicity, so the answer would be "the hard ones are big and
greasy" with no information in it. So the comparison is made WITHIN the strong half of each enzyme,
where potency is already restricted, and the ranking statistic is the metric's own per-compound
hinge rather than a squared error.

THE CONTROL IS THE INSTRUMENT. Forty-five properties times four enzymes will produce "findings" by
construction. So the same procedure runs on the thirds assigned AT RANDOM, `--null` times, and the
count of properties clearing the threshold in the real split is read against that null distribution.
A property only counts if the real split beats the null's own 95th percentile. Without that line
this file is a random-number generator with a chemistry vocabulary.

SCAFFOLDS ARE NOT AVAILABLE AS A GROUPING, and that is measured rather than assumed: the journal
records 4520 Bemis-Murcko scaffolds over 4905 molecules, so almost every compound is its own
scaffold and there is nothing to group by. The grouping here is therefore by property, plus the
Butina clusters the split already uses.

PREDICTION, fixed before the run.

    Я ожидаю, что НИЧЕГО не превысит нулевой контроль, и это не кокетство: 348 пунктов
    признаковой работы уже перебрали это представление. Если что-то и выйдет, вероятнее всего
    это размер и липофильность -- и такая находка была бы ПЕРЕОТКРЫТИЕМ, потому что MW и logP
    уже сидят в DESC, то есть модель их видит и всё равно ошибается. Настоящей находкой было бы
    свойство, которого в наших 2295 столбцах НЕТ: заряд при pH 7.4, число ионизуемых центров,
    металл-координирующие группы, стереоцентры. Их я и смотрю отдельно от остальных.

    И вторая половина: если худшие соединения сидят в НЕМНОГИХ кластерах Бутины, это пробел в
    химии; если они рассеяны по кластерам равномерно, это шум метки и никакая новая колонка
    его не уберёт.

    uv run python verify/k142_errchem.py --null 200

Reads results/preds/oof_submitted.json, data/feats.npz, data/{desc,mech}_names.csv.
Writes results/logs/k142_errchem.json. Exit: 0 ran, 2 a control failed, 3 an input is missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES

import argparse
import json

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger
from rdkit.Chem import Crippen, Descriptors, rdMolDescriptors
from scipy.stats import mannwhitneyu

from cypsplit import butina_folds, fold_digest

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
GOLDEN = "2d93c19815e14261"
# Properties computed here rather than taken from DESC, so the NEW ones are unambiguous.
# `нов` marks what our 2295 columns do not carry; those are the only ones a hit would be a
# finding about rather than a re-discovery.
PROPS = [
    ("MW", lambda m: Descriptors.MolWt(m), False),
    ("logP", lambda m: Crippen.MolLogP(m), False),
    ("TPSA", lambda m: rdMolDescriptors.CalcTPSA(m), False),
    ("HBD", lambda m: rdMolDescriptors.CalcNumHBD(m), False),
    ("HBA", lambda m: rdMolDescriptors.CalcNumHBA(m), False),
    ("ароматич_циклов", lambda m: rdMolDescriptors.CalcNumAromaticRings(m), False),
    ("вращ_связей", lambda m: rdMolDescriptors.CalcNumRotatableBonds(m), False),
    ("тяжёлых_атомов", lambda m: m.GetNumHeavyAtoms(), False),
    ("fsp3", lambda m: rdMolDescriptors.CalcFractionCSP3(m), False),
    ("циклов", lambda m: rdMolDescriptors.CalcNumRings(m), False),
    # --- the ones our matrix does not carry ---
    ("заряд_формальный", lambda m: Chem.GetFormalCharge(m), True),
    ("основных_центров", lambda m: sum(1 for a in m.GetAtoms()
                                       if a.GetSymbol() == "N" and a.GetTotalNumHs() >= 0
                                       and not a.GetIsAromatic() and a.GetDegree() <= 3), True),
    ("кислотных_групп", lambda m: len(m.GetSubstructMatches(
        Chem.MolFromSmarts("[CX3](=O)[OX2H1,OX1-]"))), True),
    ("стереоцентров", lambda m: len(Chem.FindMolChiralCenters(m, includeUnassigned=True,
                                                              useLegacyImplementation=False)), True),
    ("металл_коорд", lambda m: len(m.GetSubstructMatches(
        Chem.MolFromSmarts("[#7,#8,#16;!$([#7]C=O)]~[#6]~[#7,#8,#16]"))), True),
    ("галогенов", lambda m: sum(1 for a in m.GetAtoms()
                                if a.GetSymbol() in ("F", "Cl", "Br", "I")), True),
    ("нитро", lambda m: len(m.GetSubstructMatches(Chem.MolFromSmarts("[N+](=O)[O-]"))), True),
    ("сульфонамид", lambda m: len(m.GetSubstructMatches(Chem.MolFromSmarts("S(=O)(=O)N"))), True),
]


def hinge(p, lo, hi):
    """The metric's own per-compound error: zero inside the band, the distance outside it."""
    return np.clip(p - hi, 0, None) + np.clip(lo - p, 0, None)


def smd(a, b):
    """Standardised mean difference, pooled sd. Reported beside the p so a tiny effect on a big
    sample cannot pass as a finding."""
    sa, sb = np.std(a, ddof=1), np.std(b, ddof=1)
    n = np.sqrt(((len(a) - 1) * sa ** 2 + (len(b) - 1) * sb ** 2) / max(1, len(a) + len(b) - 2))
    return float((np.mean(a) - np.mean(b)) / n) if n > 0 else 0.0


def compare(X, worst, best, alpha):
    """How many property columns separate the two groups, and which."""
    hits = []
    for j in range(X.shape[1]):
        a, b = X[worst, j], X[best, j]
        if np.std(a) == 0 and np.std(b) == 0:
            continue
        try:
            u = mannwhitneyu(a, b, alternative="two-sided").pvalue
        except ValueError:
            continue
        if u < alpha:
            hits.append((j, float(u), smd(a, b)))
    return hits


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--null", type=int, default=200, help="сколько случайных расщеплений")
    ap.add_argument("--alpha", type=float, default=0.01)
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)

    for f in ("rows.csv", "feats.npz", "cyp-challenge-TRAIN_inhibition.csv"):
        if not _pl.Path(D + f).exists():
            print("нет %s%s" % (D, f))
            return 3
    op = _pl.Path(RES + "preds/oof_submitted.json")
    if not op.exists():
        print("нет %s" % op)
        return 3

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)
    P = [np.asarray(v, float) for v in json.load(open(op))["P"]]
    for e in range(len(CYPS)):
        if len(P[e]) != int(M[:, e].sum()):
            print("ОТКАЗ: %s предсказаний %d против маски %d"
                  % (CYPS[e], len(P[e]), int(M[:, e].sum())))
            return 2

    fold, ncl = butina_folds(list(rows.SMILES), seed=a.seed)
    dg = fold_digest(fold)
    if a.seed == 0 and dg != GOLDEN:
        print("ОТКАЗ: разбиение сдвинулось, сид 0 дал %s" % dg)
        return 2
    print("  дайджест %s, кластеров %d" % (dg, ncl), flush=True)

    mols = [Chem.MolFromSmiles(s) for s in rows.SMILES]
    if any(m is None for m in mols):
        print("ОТКАЗ: %d SMILES не разобрались" % sum(m is None for m in mols))
        return 2
    names = [p[0] for p in PROPS]
    isnew = np.array([p[2] for p in PROPS])
    X = np.zeros((len(mols), len(PROPS)))
    for j, (_nm, fn, _n) in enumerate(PROPS):
        X[:, j] = [float(fn(m)) for m in mols]
    print("  свойств %d, из них ВНЕ наших 2295 столбцов %d: %s"
          % (len(names), int(isnew.sum()), ", ".join(np.array(names)[isnew])), flush=True)

    rng = np.random.default_rng(a.seed)
    out = {"digest": dg, "alpha": a.alpha, "null_draws": a.null, "props": names,
           "is_new": isnew.tolist(), "per_enzyme": {}}
    print("\n" + "=" * 98)
    print("ЧТО ЭТО ЗА СОЕДИНЕНИЯ: внутри СИЛЬНОЙ половины, худшая треть против лучшей")
    print("=" * 98)
    for e, c in enumerate(CYPS):
        m = M[:, e]
        idx = np.where(m)[0]
        y, lo, hi = Y[m, e], LO[m, e], HI[m, e]
        h = hinge(P[e], lo, hi)
        strong = y > np.median(y)
        si = np.where(strong)[0]
        k = max(5, len(si) // 3)
        order = si[np.argsort(h[si])]
        best, worst = order[:k], order[-k:]
        Xs = X[idx]
        hits = compare(Xs, worst, best, a.alpha)

        # the null: the same split size, assigned at random inside the strong half
        nulls = []
        for _ in range(a.null):
            pr = rng.permutation(si)
            nulls.append(len(compare(Xs, pr[-k:], pr[:k], a.alpha)))
        nulls = np.asarray(nulls)
        q95 = float(np.quantile(nulls, 0.95))
        verdict = "ВЫШЕ нуля" if len(hits) > q95 else "в пределах нуля"
        print("\n  %s: сильных %d, треть %d | хинж сильной половины %.3f против слабой %.3f"
              % (c, len(si), k, float(h[si].mean()), float(h[~strong].mean())))
        print("    свойств прошло alpha=%.3g: %d | нулевой контроль: среднее %.1f, 95-й %.1f -> %s"
              % (a.alpha, len(hits), float(nulls.mean()), q95, verdict))
        for j, u, d in sorted(hits, key=lambda t: t[1])[:6]:
            print("      %-18s p %.2e, SMD %+.2f%s"
                  % (names[j], u, d, "   <- ВНЕ наших столбцов" if isnew[j] else ""))
        # where do the worst sit among the Butina clusters
        cl = fold[idx]
        wc = pd.Series(cl[worst]).value_counts(normalize=True)
        out["per_enzyme"][c] = {
            "n_strong": int(len(si)), "third": int(k),
            "hinge_strong": float(h[si].mean()), "hinge_weak": float(h[~strong].mean()),
            "hits": [{"prop": names[j], "p": u, "smd": d, "new": bool(isnew[j])}
                     for j, u, d in hits],
            "null_mean": float(nulls.mean()), "null_q95": q95,
            "above_null": bool(len(hits) > q95),
            "worst_fold_share": {str(int(kk)): float(vv) for kk, vv in wc.items()}}
        print("    худшая треть по фолдам: %s (равномерно было бы по 0.20)"
              % ", ".join("%d:%.2f" % (int(kk), vv) for kk, vv in sorted(wc.items())))

    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump(out, open(RES + "logs/k142_errchem.json", "w"), ensure_ascii=False, indent=1)
    nnew = sum(1 for c in CYPS for hh in out["per_enzyme"][c]["hits"] if hh["new"])
    above = [c for c in CYPS if out["per_enzyme"][c]["above_null"]]
    print("\n  ферментов, где находок больше нулевого контроля: %s"
          % (", ".join(above) if above else "НИ ОДНОГО"))
    print("  находок по свойствам ВНЕ наших столбцов: %d" % nnew)
    print("  (предсказание в докстроке: ожидалось НИЧЕГО выше нуля)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
