"""Does the per-enzyme screening arm add to the ensemble, or is it the same information again.

Item 177 measured the arm this asks about: a per-enzyme model trained on its own curves plus the
11505 screening measurements for that enzyme is worth **+0.028 of rank on four seeds** over the same
model without them, with every enzyme above its own floor. That is the largest feature-side gain in
the file after the dead zone.

Item 176 is the reason this cannot be assumed to reach the submission. The trunk improved by 0.031
alone and contributed nothing to the five-member ensemble, because it improved by moving toward what
the other members already fit. An ensemble pays for disagreement. So the question here is not how
good the arm is but how *different* it is.

Composition only, from saved predictions, in minutes: both `src/submit.py` (`_combine`) and
`src/abldzens.py` combine members by an unweighted mean, so adding a member is exact arithmetic on
what is already on disk.

    четыре                     the existing ensemble
    четыре + скрининг          the same plus the per-enzyme screening arm as a fifth
    четыре, мёртвая зона       item 162's arm, the current best
    ... + скрининг             and the same with the screening arm added

Pre-registered. The arm scores 0.590 alone against the four-member ensemble's 0.601, so it is a
*weaker* member than the ensemble it joins; by item 92's regime that is not a reason to expect it to
fail, since the Gaussian process is also weaker alone and helps. What decides is correlation with
the existing members' errors, which is printed alongside so a null has an explanation rather than
just a number.

Reads oof_aux.json, oof_aux2.json, oof_aux123.json and the dzens files. Writes nothing.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
ARM = "независимо+скрининг"
# Второй кандидат в члены: аддитивная модель Фри-Вилсона на всех битах. Пункт 191 намерил
# у неё профиль по стратам, ОБРАТНЫЙ профилю GP из пункта 97 --- она вредит на дальних
# соседях и помогает на ближних. Именно противоположность профилей и есть та
# дополнительность, которой ансамблю не хватает по пунктам 176 и 182.
FW = ("oof_fw.json", "Фри-Вилсон, все 2048")


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())

    S = {}
    if os.path.exists(RES + "preds/" + FW[0]):
        _f = json.load(open(RES + "preds/" + FW[0]))["preds"]
        S.update({k: v for k, v in _f.items() if FW[1] in k})
    for fn in ("oof_aux.json", "oof_aux2.json", "oof_aux123.json"):
        p = RES + "preds/" + fn
        if os.path.exists(p):
            S.update(json.load(open(p))["preds"])
    E = {}
    for fn in ("oof_dzens.json", "oof_dzens123.json"):
        p = RES + "preds/" + fn
        if os.path.exists(p):
            E.update(json.load(open(p))["preds"])
    seeds = sorted({int(k.split("|")[0]) for k in E}
                   & {int(k.split("|")[0]) for k in S if ARM in k})
    print(f"сиды, где есть и ансамбль, и скрининговая рука: {seeds}\n")
    if not seeds:
        print("нет пересечения --- нечего складывать")
        return

    rec = []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for ens in ("базовый ансамбль", "мёртвая зона везде"):
            for add in (None, ARM, FW[1]):
                r = {"seed": seed,
                     "состав": ens + ("" if add is None else
                                      " + скрининг" if add == ARM else " + Фри-Вилсон")}
                ok = True
                for c in CYPS:
                    ke = f"{seed}|{ens}|{c}"
                    ks = None if add is None else f"{seed}|{add}|{c}"
                    if ke not in E or (add is not None and ks not in S):
                        ok = False
                        break
                    col = f"{c}_pIC50_direct_inhibition"
                    m = tr[col].notna().to_numpy()
                    y = tr.loc[m, col].to_numpy()
                    lo = tr.loc[m, col + "_conf_low"].to_numpy()
                    hi = tr.loc[m, col + "_conf_high"].to_numpy()
                    a = np.asarray(E[ke], float)
                    if len(a) != int(m.sum()):
                        ok = False
                        break
                    if add is not None:
                        b = np.asarray(S[ks], float)
                        # Ансамбль --- невзвешенное среднее ЧЕТЫРЁХ членов, поэтому
                        # пятый добавляется как (4a + b)/5.
                        p = (4.0 * a + b) / 5.0
                        r[f"{c} corr"] = float(np.corrcoef(a - y, b - y)[0, 1])
                    else:
                        p = a
                    q = fit_apply(p, lo, hi, fold[m], np.ones(len(y)) / len(y))
                    r[f"{c} пара"] = float(strae(y, q, y_true_upper=hi, y_true_lower=lo))
                    r[f"{c} rho"] = float(spearmanr(y, p).statistic)
                if not ok:
                    continue
                for t in ("пара", "rho"):
                    r[f"MACRO {t}"] = float(np.mean([r[f"{c} {t}"] for c in CYPS]))
                rec.append(r)

    df = pd.DataFrame(rec)
    g = df.groupby("состав", sort=False)[["MACRO пара", "MACRO rho"]
                                         + [f"{c} rho" for c in CYPS]].mean()
    print(g.round(4).to_string())
    cc = [c for c in df.columns if c.endswith("corr")]
    if cc:
        print("\nкорреляция ошибок скрининговой руки с ошибками ансамбля:")
        print("  " + "  ".join(f"{c.split()[0][3:]} {df[c].mean():.3f}" for c in cc))
    print(f"\nсидов усреднено: {df.groupby('состав').size().to_dict()}")
    print("""
Как читать. Решает разность строк «+ скрининг» и без него. Корреляция ошибок внизу объясняет
исход: близкая к единице означает, что рука ошибается там же, где ансамбль, и усреднение с ней
ничего не даёт, каким бы хорошим ни был её собственный ранг. Пункт 176 --- ровно этот случай на
стволе.""")


if __name__ == "__main__":
    main()
