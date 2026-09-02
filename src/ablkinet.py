"""The kinetic link: undo the turnover term with pinned physics before fitting, put it back after.

An inhibitor that is also a substrate is emptied out of the EI complex by turnover as well as by
dissociation, so its apparent potency is weaker than its binding affinity. The steady state of

    E + S <-> ES -> E + P            the probe, whose rate we measure
    E + I <-> EI -> E + I'           our compound, if it is also a substrate

gives, for competitive inhibition,

    pIC50_obs = pKi - log10(1 + kcat_I/koff) - log10(1 + [S]/Km)
                      +----- per compound ----+   +-- per enzyme --+

The third term is a constant per enzyme and is already absorbed by `calshift`'s fitted `d_e`
(item 171). The second is per compound, and this file is about it.

**Why subtracting it from the prediction would be arithmetic nonsense, which is the first thing
to get right.** The labels are `pIC50_direct_inhibition` -- the OBSERVED potency, turnover
already inside it. Subtracting delta from a model trained on those labels removes it twice. The
scheme says something else: the latent is pKi and the label is its corrupted image. So the
correction runs the other way,

    train on   y + delta   (the reconstructed latent)
    predict    f(x)
    report     f(x) - delta

and the gain exists exactly when **pKi is a smoother function of structure than pIC50 is** --
when turnover adds to the label a chemically foreign term that pinned physics can take back out.
If the learner were perfect this would change nothing at all, since delta is itself a function of
structure; the whole effect is learnability. That makes it falsifiable rather than tautological.

**Pinning, and why the letter's own constant is wrong.** Bell-Evans-Polanyi says the barrier
tracks the C-H bond strength, and transition-state theory turns a barrier into a rate. But
pinning the slope at 1/(RT ln10) assumes the entire barrier difference reaches the rate, which
over SMARTCyp's 2.5-to-75 kJ/mol span gives deltas up to 9 log units against a pIC50 range of
about five. The BEP slope for hydrogen abstraction is nearer 0.4. So the SHAPE is pinned and one
global amplitude is left:

    delta_i = A * log10(1 + 10^{-(E_i - E0)/w}),  E0 = 49.83 (median barrier), w = 14.3

with `w = RT ln10 / 0.4` from theory and `E0` the median, so that turnover balances dissociation
at the typical compound. `A` is the single free number the physics does not supply.

**And it is swept rather than fitted, because fitting it and reporting the best is the result.**
`A = 0` is the control arm and is the same estimator on the same folds, so the comparison is a
paired one. Item 178 is the standing warning: a more flexible link described the instrument
better and made the model worse, because flexibility finds ways to satisfy a reading without
moving the prediction.

**Pre-registration, from the mechanism rather than from the fit.** The turnover branch is empty
for haem coordinators: azoles, pyridines and pyrimidines sit on the iron and do not turn over. So

    эффект живёт на 1400 НЕкоординаторах, около нуля на 3505 с [nX2]  -> механизм прочитан верно
    эффект равномерен по обеим популяциям                            -> это колонка про размер
    эффект сосредоточен на координаторах                             -> механизм наоборот, закрыто

The split is item 184's, and the counts were re-checked here: 3505 and 1400 of 4905.

E comes from `Score_min` in `data/smartcyp.npz` -- SMARTCyp's activation energy for the most
labile site, in kJ/mol on an absolute inter-molecular scale (item 208 established that the scale
is absolute and not a ranking, which is why ALFABET would add nothing here).

Same folds, learner and metric as src/ablate.py. Writes results/preds/oof_kinet.json.
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
from scipy.stats import spearmanr
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from rdkit import Chem, RDLogger

from cypsplit import butina_folds
from shrinkchoice import fit_apply

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
E0, W = 49.83, 14.3          # медианный барьер; RT ln10 / beta_BEP при beta = 0.4


def shape(E):
    """Форма оборотного члена: убывающая, насыщающаяся по барьеру. Амплитуда снаружи."""
    return np.log10(1.0 + np.power(10.0, -(E - E0) / W))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--amps", default="0,0.3,0.6")
    ap.add_argument("--masks", default="все", help="все|некоорд, через |")
    ap.add_argument("--shuffle", action="store_true", help="контроль: перемешать delta")
    ap.add_argument("--out", default=RES + "preds/oof_kinet.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    amps = [float(x) for x in a.amps.split(",")]
    masks = a.masks.split("|")

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)
    sc = np.load(D + "smartcyp.npz", allow_pickle=True)
    nm = [str(x) for x in sc["names"]]
    E = sc["train"][:, nm.index("Score_min")].astype(float)
    E = np.where(np.isfinite(E), E, np.nanmedian(E))     # одна строка NaN

    patt = Chem.MolFromSmarts("[nX2]")
    coord = np.array([(m := Chem.MolFromSmiles(s)) is not None and m.HasSubstructMatch(patt)
                      for s in rows.SMILES])
    d_raw = shape(E)
    print(f"барьер: [{E.min():.1f}, {E.max():.1f}] кДж/моль, медиана {np.median(E):.1f}")
    print(f"форма delta при A=1: медиана {np.median(d_raw):.3f}, "
          f"90-й проц {np.percentile(d_raw, 90):.3f}, макс {d_raw.max():.3f}")
    print(f"координаторы [nX2]: {int(coord.sum())}, некоординаторы {int((~coord).sum())}\n",
          flush=True)

    Y, LO, HI, M = {}, {}, {}, {}
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        M[c] = tr[col].notna().to_numpy()
        Y[c] = tr[col].to_numpy(float)
        LO[c] = tr[col + "_conf_low"].to_numpy(float)
        HI[c] = tr[col + "_conf_high"].to_numpy(float)

    out, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        rng = np.random.default_rng(1000 + seed)
        for mk in masks:
            base = d_raw.copy()
            if mk == "некоорд":
                # Оборот невозможен там, где молекула сидит на железе: ветка пуста.
                base[coord] = 0.0
            if a.shuffle:
                base = base[rng.permutation(len(base))]
            for A in amps:
                t0 = time.time()
                dlt = A * base
                r = {"seed": seed, "рука": f"A={A} {mk}" + (", перемешан" if a.shuffle else "")}
                for e, c in enumerate(CYPS):
                    m = M[c]
                    y, lo, hi, fi = Y[c][m], LO[c][m], HI[c][m], fold[m]
                    d, Xi = dlt[m], X[m]
                    p = np.zeros(len(y))
                    for f in range(5):
                        trn, te = fi != f, fi == f
                        if te.sum() == 0:
                            continue
                        # Обучаем на ЛАТЕНТНОЙ величине y + delta, предсказываем, снимаем delta.
                        mdl = HistGradientBoostingRegressor(**KW).fit(Xi[trn], (y + d)[trn])
                        p[te] = mdl.predict(Xi[te]) - d[te]
                    q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                    out[f"{seed}|{r['рука']}|{c}"] = p.tolist()
                    r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi,
                                                       y_true_lower=lo)), 4)
                    r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
                for t in ("пара", "rho"):
                    r[f"MACRO {t}"] = round(float(np.mean([r[f"{c} {t}"] for c in CYPS])), 4)
                table.append(r)
                print(f"  сид {seed} {r['рука']:26s} пара {r['MACRO пара']:.4f} "
                      f"rho {r['MACRO rho']:.4f}  "
                      + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                      + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                         + [f"{c} rho" for c in CYPS]].mean().to_string())
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. A = 0 --- контроль: тот же оценщик, те же фолды, поправки нет. Всё сравнение
парное относительно него.

Если ни одно A не поднимает ранг выше поферментного пола (1A2 0.0061, 2C9 0.0071, 2D6 0.0049,
3A4 0.0033), оборотный член не несёт посоединениевого порядка, и кинетическая связь закрыта
тем же нулём, что SMARTCyp-колонка в пункте 179 --- то есть фиксированная форма не оказалась
эффективнее свободной.

Если поднимает --- решает НЕ величина, а разница между «все» и «некоорд». Механизм предсказывает,
что зануление delta на 3505 координаторах гема должно УЛУЧШИТЬ результат: у них ветка оборота
пуста, и ненулевая поправка там есть чистый шум. Если «все» и «некоорд» равны, поправка работает
как ещё одна колонка про размер, а не как оборот.""")


if __name__ == "__main__":
    main()
