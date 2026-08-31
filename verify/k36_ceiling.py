"""What does one screening reading know that the whole structural model does not.

Why this question and not another improvement. Across 130 items no intervention has moved macro
ST-RAE by more than about 0.03, against a noise floor of 0.007. The distance from our model to a
single calibrated screening column is 0.354 macro, and per enzyme:

    фермент   модель+пара   скрининг один   разрыв
    CYP1A2         0.8128          0.3106    0.502
    CYP2C9         0.6559          0.3154    0.340
    CYP2D6         0.9052          0.5992    0.306
    CYP3A4         0.4860          0.2168    0.269

We have been optimising the third decimal of a quantity whose reachable range is in the first. So
before spending the remaining weeks on more arms, it is worth asking what the gap is *made of* --
and this is the one analysis in the queue that can answer "stop".

The design. Take the compounds where the model is badly wrong and the screening column is right.
Those are exactly the compounds carrying the gap. Then ask what they are, on axes chosen in
advance so the answer cannot be read into the data afterwards:

  **chemical**       structural novelty (distance to the training set), size and lipophilicity,
                     mechanistic-block composition, ring systems. If the residue looks like a
                     chemotype, there is something the model cannot see and finding it is worth
                     more than any arm in the queue.
  **non-chemical**   the aggregation flag from item 102, assay position, Emax depth, band width.
                     If the residue is about how a compound *behaves in this assay* rather than
                     what it binds, then prediction from structure has a ceiling near where we
                     already stand.

The second outcome is the valuable one and the one nobody wants. "0.68 is close to the limit of
what structure predicts here" is a stronger statement than another 0.01, and it is only sayable
with a calibrated reference in hand -- which the screen provides and a hypothesis does not. Item
102 walked part of this path with the aggregation suspects and found chemistry rather than
solubility; this is the same question asked of the whole residue instead of one flag.

Both predictions are pre-registered, so a mixed answer is also readable: a residue that is
chemical on CYP1A2 and behavioural on CYP3A4 would say the ceiling is per-enzyme, which is itself
worth knowing.

Reads results/preds/oof.json and the screening table. Writes nothing. ~3 minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from cypsplit import butina_folds
from shrinkchoice import fit_apply

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
FRAC = 0.15        # какую долю берём в «остаток»


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    sc = pd.read_csv(D + "cyp-challenge-single-concentration-TRAIN.csv")
    em = pd.read_csv(D + "cyp-challenge-TRAIN_Emax.csv").set_index("Molecule_Name")
    z = np.load(D + "feats.npz")
    dn = [l.strip() for l in open(D + "desc_names.csv")]
    mn = [l.strip() for l in open(D + "mech_names.csv")]
    O = json.load(open(RES + "preds/oof.json"))
    fold, _ = butina_folds(list(rows.SMILES))

    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]
    nn = np.zeros(len(rows))
    for f in range(5):
        te, trn = np.where(fold == f)[0], np.where(fold != f)[0]
        ref = [fps[j] for j in trn]
        for i in te:
            nn[i] = max(DataStructs.BulkTanimotoSimilarity(fps[i], ref))

    logp = z["DESC"][:, dn.index("MolLogP")].astype(float)
    mw = z["DESC"][:, dn.index("MolWt")].astype(float) if "MolWt" in dn else np.zeros(len(rows))
    narom = z["DESC"][:, dn.index("NumAromaticRings")].astype(float) \
        if "NumAromaticRings" in dn else np.zeros(len(rows))
    nbas = z["MECH"][:, mn.index("n_basicN_ali")].astype(float)

    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy()
        lo = tr.loc[m, col + "_conf_low"].to_numpy()
        hi = tr.loc[m, col + "_conf_high"].to_numpy()
        fi, u = fold[m], np.ones(int(m.sum())) / int(m.sum())
        p = fit_apply(np.asarray(O[f"FP+DESC+MECH|{c}"], float), lo, hi, fi, u)

        sub = sc[sc.enzyme == c].set_index("Molecule_Name")
        l2 = rows.set_index("Molecule_Name").join(sub["log2fc_estimate"])["log2fc_estimate"] \
            .to_numpy(float)[m]
        ok = ~np.isnan(l2)
        q = np.full(len(y), np.nan)
        for f in range(5):
            te, trn = (fi == f) & ok, (fi != f) & ok
            if te.sum() == 0 or trn.sum() < 50:
                continue
            q[te] = IsotonicRegression(increasing=False, out_of_bounds="clip").fit(
                l2[trn], y[trn]).predict(l2[te])
        g = ~np.isnan(q)

        e_mod = np.maximum(0.0, np.maximum(lo - p, p - hi))
        e_scr = np.maximum(0.0, np.maximum(lo - q, q - hi))
        adv = np.where(g, e_mod - e_scr, -np.inf)      # насколько скрининг лучше модели
        k = int(FRAC * g.sum())
        sel = np.zeros(len(y), bool)
        sel[np.argsort(-adv)[:k]] = True
        rest = g & ~sel

        emx = em[f"{c}_EmaxVsPosCtrl_direct_inhibition"].reindex(
            rows.Molecule_Name[m]).to_numpy(float)
        w = hi - lo
        print(f"\n=== {c}: остаток --- {k} соединений, где скрининг выигрывает больше всего "
              f"(средний выигрыш {adv[sel].mean():.3f} pIC50) ===")
        print(f"{'ось':26s} {'остаток':>10s} {'прочие':>10s} {'разница':>9s}")
        for nm, v in (("сходство с обучением", nn[m]), ("logP", logp[m]),
                      ("молекулярная масса", mw[m]), ("ароматических колец", narom[m]),
                      ("основных N", nbas[m]), ("ширина полосы", w),
                      ("Emax", emx), ("метка pIC50", y)):
            a, b = np.nanmean(v[sel]), np.nanmean(v[rest])
            print(f"{nm:26s} {a:10.3f} {b:10.3f} {a-b:+9.3f}")

    print("""
Как читать. Оси разделены заранее на две группы, и вывод читается по тому, какая группа
разошлась.

Химические --- сходство с обучением, размер, липофильность, кольца, основания. Если остаток
отличается по ним, значит есть хемотип, которого модель не видит, и его поиск стоит больше
любой руки в очереди.

Поведенческие --- ширина полосы, Emax, метка. Если остаток отличается только по ним, то он
про то, как соединение ведёт себя В ЭТОМ анализе, а не про его сродство, и предсказание из
структуры упирается в потолок примерно там, где мы стоим. Это и есть исход, ради которого
файл написан: сказать «дальше некуда» можно только с калиброванным эталоном в руках.""")


if __name__ == "__main__":
    main()
