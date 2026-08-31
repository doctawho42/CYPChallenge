"""Does a large pretrained representation carry anything our matrix does not.

The question is open despite item 68, and the reason is in `src/embed.py`'s own docstring. That
script chose chemprop `rdkit2d` deliberately, a checkpoint pretrained to predict the RDKit
descriptors that make up our DESC block. It is the right choice for "does a learned
representation beat handing the same information over directly" and the wrong one for "does a
pretrained representation know anything we do not", because an encoder trained to reproduce our
features has little room to add to them. Item 122 then measured that its concatenated arm loses
under the affine pair, and that closed the narrow question, not this one.

`src/embed2.py` supplies the other kind: transformers pretrained on raw structure at scale, with
no descriptor supervision anywhere in the objective.

Four arms per checkpoint, and the ordering of the first three is the whole design:

  база                 FP+DESC+MECH, 2295 columns. Must reproduce 0.7150 / 0.5651 at seed 0 or
                       nothing below is comparable to anything above.
  эмбеддинг один       the substitution arm. Expected to lose, and item 101 records this as the
                       third time a representation lost as a replacement and won as an addition,
                       so a loss here says nothing about the arm that follows.
  база + эмбеддинг     the concatenation. This is the arm the question rests on.
  база + эмбеддинг-64  the same, reduced to 64 components. 768 extra columns against about 1400
                       labelled rows per enzyme is a great deal of width, and item 68 already
                       found width to matter on this data; without this arm a loss in the third
                       could be dimensionality rather than content.

The reduction is fitted inside each fold on the training rows only. It is unsupervised, so the
leak would be mild, but "mild" is not a standard this repository uses.

Same folds, masks, metric and learner settings as src/ablate.py. Reads data/feats.npz,
data/rows.csv and an embedding written by src/embed2.py. Writes results/preds/oof_pre_<slug>.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import re
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.decomposition import PCA
from sklearn.ensemble import HistGradientBoostingRegressor
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
          l2_regularization=1.0, random_state=0)
NCOMP = 64


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--emb", required=True, help="npz от src/embed2.py")
    ap.add_argument("--seeds", default="0")
    ap.add_argument("--out", default="")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    slug = re.sub(r"[^a-z0-9]+", "_", _pl.Path(a.emb).stem.replace("emb_", "").lower()).strip("_")
    out = a.out or RES + f"preds/oof_pre_{slug}.json"

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)
    ez = np.load(a.emb, allow_pickle=True)
    E = ez["train"].astype(np.float32)
    if len(E) != len(X):
        raise SystemExit(f"эмбеддинг не той длины: {len(E)} против {len(X)}")
    print(f"эмбеддинг {a.emb}: {E.shape}, модель {ez['model']}, пулинг {ez['pool']}")
    print(f"база {X.shape}\n")

    out_p, table = {}, []
    for seed in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        for nm in ("база", "эмбеддинг один", "база + эмбеддинг", f"база + эмбеддинг-{NCOMP}"):
            t0 = time.time()
            r = {"seed": seed, "рука": nm}
            for c in CYPS:
                col = f"{c}_pIC50_direct_inhibition"
                m = tr[col].notna().to_numpy()
                y = tr.loc[m, col].to_numpy()
                lo = tr.loc[m, col + "_conf_low"].to_numpy()
                hi = tr.loc[m, col + "_conf_high"].to_numpy()
                Xi, Ei, fi = X[m], E[m], fold[m]
                p = np.zeros_like(y)
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    if nm.endswith(str(NCOMP)):
                        # Понижение фитится на обучающих строках фолда: без надзора, но
                        # «слабая утечка» --- не тот стандарт, которым здесь пользуются.
                        pca = PCA(n_components=NCOMP, random_state=0).fit(Ei[trn])
                        Etr, Ete = pca.transform(Ei[trn]), pca.transform(Ei[te])
                        A, B = np.hstack([Xi[trn], Etr]), np.hstack([Xi[te], Ete])
                    elif nm == "эмбеддинг один":
                        A, B = Ei[trn], Ei[te]
                    elif nm == "база + эмбеддинг":
                        A, B = np.hstack([Xi[trn], Ei[trn]]), np.hstack([Xi[te], Ei[te]])
                    else:
                        A, B = Xi[trn], Xi[te]
                    p[te] = HistGradientBoostingRegressor(**KW).fit(A, y[trn]).predict(B)
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out_p[f"{seed}|{nm}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for tag in ("пара", "rho"):
                r[f"MACRO {tag}"] = round(float(np.mean([r[f"{c} {tag}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {nm:22s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    print()
    print(df.groupby("рука", sort=False)[["MACRO пара", "MACRO rho"]
                                         + [f"{c} rho" for c in CYPS]].mean().to_string())
    json.dump({"table": table, "preds": out_p}, open(out, "w"))
    print(f"\nсохранено: {out}")
    print("""
Как читать. Первая рука обязана дать 0.7150 / 0.5651 на сиде 0.

Смотреть на ранг. Пункт 77 и здесь решает: прирост сырого ST-RAE аффинная пара перепишет, и
именно так закрылась предыдущая попытка с энкодером (пункт 122) --- сырое улучшалось, ранг
падал, после пары становилось хуже.

Вторая рука почти наверняка проиграет, и это ничего не значит: замена проигрывала уже трижды
(энкодер, FCFP, механистический блок) и дважды из трёх добавка выигрывала. Решает третья.

Если третья проигрывает, а четвёртая выигрывает --- дело было в ширине, и вывод касается
размерности, а не предобучения. Если проигрывают обе --- представление действительно ничего
не добавляет к нашей матрице, и это первый честный ответ на этот вопрос в журнале.""")


if __name__ == "__main__":
    main()
