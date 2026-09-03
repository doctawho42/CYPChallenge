"""Active-site blocks for CYP1A2 and CYP3A4, the two enzymes that have never had one.

The observation this rests on, and it is arithmetic rather than analogy. Item 167 read every saved
ablation per enzyme against that enzyme's own floor from item 165:

    вмешательство            сидов      1A2      2C9      2D6      3A4     пол
    механистический блок         4  -0.0019  +0.0148  +0.0454  +0.0008
    пулирование                  4  +0.0157  +0.0071  +0.0374  -0.0057
    блок формы                   4  +0.0060  -0.0003  +0.0087  -0.0008
    блок кислот                  4  +0.0012  -0.0023  -0.0025  -0.0006
    пол этого фермента                0.0061   0.0071   0.0049   0.0033

**CYP3A4 clears its floor on nothing at all; CYP1A2 clears it on pooling alone.** And `src/feats.py`
opens with a docstring reading "Mechanistic feature block for CYP2D6" -- the one enzyme with a
deliberately built site block is the one enzyme every intervention helps. CYP2C9's +0.0148 from the
same block is an accident of `n_acid` and `ph74_n_anion` meeting Arg108.

Item 139 makes the misallocation plain: the distance to a single calibrated screening column is
0.502 on CYP1A2, 0.340 on CYP2C9, 0.306 on CYP2D6 and 0.289 on CYP3A4. Effort has gone to the
enzyme with the second smallest gap.

What each block encodes, and why these columns and not others.

**CYP1A2 -- a narrow planar slot.** Its cavity is about 375 cubic angstroms, the smallest of the
four, lined with phenylalanines and shaped like a slit. It takes flat fused aromatics and refuses
anything thick. This is not imported from the literature: item 156 measured it here. Of nine
fragments decoded from the fingerprint, the **tertiary aliphatic amine costs -0.58 of pIC50 on
CYP1A2 after controlling for size, lipophilicity, aromaticity, TPSA, acceptors and sp3 fraction --
and does nothing on the other three enzymes.** sp3 basicity being rejected by a slit is a feature
specification, and the columns follow from it: thinness, flatness, fused-aromatic fraction, sp3
branching, and the presence of the offending amine.

**CYP3A4 -- a large flexible cavity.** About 1400 cubic angstroms, the largest, able to hold two
ligands at once, and famous for it. The columns are volume against that cavity, surface, radius of
gyration, flexibility, and the number of separate aromatic systems a molecule presents, since a
two-moiety ligand is what a two-site cavity is unusual for.

**CYP2C9 and CYP2D6 are controls, not proposals.** CYP2D6's block exists and must reproduce its
+0.045, or the harness is broken. CYP2C9 gets the acid columns that item 167 already measured at
-0.0023, so if they come back positive here something is wrong with this script rather than with
item 167.

Arms, and the third is the one that can kill the idea cheaply:

    база                    FP+DESC+MECH, поферментно
    +сайт                   плюс блок СВОЕГО фермента
    +сайт перемешанный      те же колонки, перемешанные между молекулами. Ширина матрицы,
                            маргиналы и число колонок те же; разорвана связь с молекулой.
    +форма целиком          все шестнадцать колонок shape3d каждому ферменту без разбора.
                            Если это не хуже целевого блока, то целиться было незачем и
                            ответ --- "добавьте форму всем", а не "стройте блоки центров".

Pre-registered per enzyme against its own floor from item 165, because item 165 exists precisely so
that this is not judged on the macro: CYP1A2 must clear **0.0061** and CYP3A4 **0.0033**. The
permutation arm must clear neither. A gain on CYP2D6 from its own block that is smaller than
+0.045 means the harness differs from the ablation grid and the numbers are not comparable to it.

The learner is the plain-tree booster pinned to `src/ablpairloss.py`: 200 trees, depth 5,
max_features 0.3, which item 140 measured at 0.5687 macro against HistGB's 0.5651 and sixty times
faster. Four seeds are affordable here only because of that.

Reads data/feats.npz, data/rows.csv, data/shape3d.npz, data/acid.npz.
Writes results/preds/oof_site.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time

import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from sklearn.tree import DecisionTreeRegressor
from rdkit import Chem, RDLogger
from rdkit.Chem import Descriptors
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply

RDLogger.DisableLog("rdApp.*")
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
NTREE, LR, DEPTH, MAXFEAT = 200, 0.06, 5, 0.3      # пины из src/ablpairloss.py
FLOOR = {"CYP1A2": 0.0061, "CYP2C9": 0.0071, "CYP2D6": 0.0049, "CYP3A4": 0.0033}
CAV3A4 = 1400.0        # объём полости CYP3A4, кубические ангстремы
TERT_AMINE = Chem.MolFromSmarts("[NX3;H0;!$(N[#6]=[O,N,S]);!$(N[a]);!$(N[SX4](=O)=O);!$(N#*);"
                                "!$([N+]);!$(N=*)]([CX4])([CX4])[CX4]")


def arom_systems(mol):
    """Число РАЗДЕЛЬНЫХ ароматических систем и доля тяжёлых атомов в крупнейшей.

    Слитые кольца считаются одной системой: объединяются те, что делят атом.
    Для CYP3A4 важно число систем (двухлигандная полость), для CYP1A2 --- доля в крупнейшей
    (плоское слитое ядро против болтающихся заместителей).
    """
    ri = mol.GetRingInfo()
    rings = [set(r) for r in ri.AtomRings()
             if all(mol.GetAtomWithIdx(i).GetIsAromatic() for i in r)]
    comps = []
    for r in rings:
        hit = [c for c in comps if c & r]
        merged = set(r)
        for c in hit:
            merged |= c
            comps.remove(c)
        comps.append(merged)
    n_heavy = mol.GetNumHeavyAtoms()
    big = max((len(c) for c in comps), default=0)
    return len(comps), (big / n_heavy if n_heavy else 0.0)


def site_blocks(smiles, shape, sh_names, acid):
    """Поферментные колонки активного центра. Возвращает dict фермент -> матрица."""
    idx = {n: i for i, n in enumerate(sh_names)}
    pbf = shape[:, idx["PBF"]]
    npr1, npr2 = shape[:, idx["NPR1"]], shape[:, idx["NPR2"]]
    pmi1, pmi3 = shape[:, idx["PMI1"]], shape[:, idx["PMI3"]]
    asph, spher = shape[:, idx["Asphericity"]], shape[:, idx["SpherocityIndex"]]
    vol, sasa = shape[:, idx["vol_vdw"]], shape[:, idx["sasa"]]
    rgyr, nrot = shape[:, idx["RadiusOfGyration"]], shape[:, idx["n_rot_frac"]]

    n_sys, frac_big = [], []
    n_branch, n_tert, fsp3 = [], [], []
    for s in smiles:
        m = Chem.MolFromSmiles(s)
        if m is None:
            n_sys.append(0); frac_big.append(0.0); n_branch.append(0)
            n_tert.append(0); fsp3.append(0.0); continue
        a, b = arom_systems(m)
        n_sys.append(a); frac_big.append(b)
        n_branch.append(sum(1 for at in m.GetAtoms()
                            if at.GetSymbol() == "C" and at.GetHybridization().name == "SP3"
                            and at.GetDegree() >= 3))
        n_tert.append(len(m.GetSubstructMatches(TERT_AMINE)) if TERT_AMINE else 0)
        fsp3.append(Descriptors.FractionCSP3(m))
    n_sys = np.array(n_sys, float); frac_big = np.array(frac_big, float)
    n_branch = np.array(n_branch, float); n_tert = np.array(n_tert, float)
    fsp3 = np.array(fsp3, float)
    thin = pmi1 / np.maximum(pmi3, 1e-9)

    return {
        # Узкая плоская щель: тонкость, плоскостность, слитое ядро, sp3-помеха, амин из 156.
        "CYP1A2": (np.column_stack([pbf, npr1, npr2, thin, asph, spher,
                                    frac_big, n_branch, n_tert, fsp3]).astype(np.float32),
                   ["PBF", "NPR1", "NPR2", "PMI1/PMI3", "Asphericity", "Spherocity",
                    "доля в крупнейшем ядре", "sp3-ветвлений", "трет. алиф. аминов", "FractionCSP3"]),
        # Большая гибкая полость: объём против 1400 A^3, поверхность, гибкость, две системы.
        "CYP3A4": (np.column_stack([vol, vol / CAV3A4, sasa, rgyr, nrot,
                                    n_sys, frac_big]).astype(np.float32),
                   ["vol_vdw", "vol/1400", "SASA", "RadiusOfGyration", "n_rot_frac",
                    "ароматических систем", "доля в крупнейшем ядре"]),
        # Контроль: пункт 167 намерил эти колонки в -0.0023 на 2C9.
        "CYP2C9": (acid.astype(np.float32), [f"acid{i}" for i in range(acid.shape[1])]),
        # Контроль: блок 2D6 уже в MECH; сюда идут его трёхмерные углы.
        "CYP2D6": (shape[:, [idx["bN_arom_ang_min"], idx["bN_arom_ang_mean"]]].astype(np.float32),
                   ["bN_arom_ang_min", "bN_arom_ang_mean"]),
    }


def boost(Xtr, ytr, Xte, seed):
    rng = np.random.default_rng(seed)
    s = np.full(len(ytr), float(ytr.mean()))
    pred = np.full(len(Xte), s[0])
    for _ in range(NTREE):
        t = DecisionTreeRegressor(max_depth=DEPTH, max_features=MAXFEAT,
                                  random_state=int(rng.integers(1 << 30)))
        t.fit(Xtr, ytr - s)
        s = s + LR * t.predict(Xtr)
        pred = pred + LR * t.predict(Xte)
    return pred


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--arms", default="база|+сайт|+сайт перемешанный|+форма целиком")
    ap.add_argument("--out", default=RES + "preds/oof_site.json")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",")]
    arms = a.arms.split("|")

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]]).astype(np.float32)
    shape = np.load(D + "shape3d.npz")["train"].astype(np.float64)
    QUANT = None
    if _pl.Path(D + "quantum.npz").exists():
        _q = np.load(D + "quantum.npz", allow_pickle=True)
        QUANT = np.nan_to_num(_q["train"], nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
    SCYP = None
    if _pl.Path(D + "smartcyp.npz").exists():
        _z = np.load(D + "smartcyp.npz", allow_pickle=True)
        SCYP = np.nan_to_num(_z["train"], nan=0.0, posinf=0.0, neginf=0.0).astype(np.float32)
    sh_names = [l.strip() for l in open(D + "shape_names.csv")]
    acid = np.load(D + "acid.npz")["A"].astype(np.float64)
    SB = site_blocks(list(rows.SMILES), shape, sh_names, acid)
    SHAPE_ALL = np.nan_to_num(shape, posinf=0.0, neginf=0.0).astype(np.float32)

    print(f"база {X.shape[1]} колонок; блоки центров:")
    for c in CYPS:
        B, nm = SB[c]
        print(f"  {c}: {B.shape[1]:2d} колонок  {', '.join(nm)}")
    print(f"\nпол по ферментам (пункт 165): " + "  ".join(f"{c[3:]} {FLOOR[c]:.4f}" for c in CYPS))
    print()

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
        for arm in arms:
            t0 = time.time()
            r = {"seed": seed, "рука": arm}
            for c in CYPS:
                m = M[c]
                y, lo, hi, fi = Y[c][m], LO[c][m], HI[c][m], fold[m]
                if arm == "база":
                    Xi = X[m]
                elif arm == "+форма целиком":
                    Xi = np.hstack([X[m], SHAPE_ALL[m]])
                elif arm.startswith("+квант"):
                    # Электронная структура из GFN2-xTB: HOMO, LUMO, щель, диполь, заряд на
                    # основном и ароматическом азоте, функция Фукуи f-, свободный конус.
                    # Проходит ворота пункта 166: ни одна из этих величин не выводится из
                    # 217 дескрипторов RDKit. И пункт 156 даёт острую предрегистрацию ---
                    # все выжившие фрагменты фингерпринта суть sp2-азот, координирующий
                    # железо гема, а f- и заряд на азоте это ровно та электроника, которая
                    # должна предсказывать силу координации.
                    if QUANT is None:
                        raise SystemExit("нет data/quantum.npz")
                    B = QUANT[m]
                    if "перемешанный" in arm:
                        B = B[np.random.default_rng(seed * 23 + 7).permutation(len(B))]
                    Xi = np.hstack([X[m], B])
                elif arm.startswith("+SMARTCyp"):
                    if SCYP is None:
                        raise SystemExit("нет data/smartcyp.npz --- сначала src/smartcyp.py")
                    B = SCYP[m]
                    if "перемешанный" in arm:
                        B = B[np.random.default_rng(seed * 13 + 5).permutation(len(B))]
                    Xi = np.hstack([X[m], B])
                else:
                    B = np.nan_to_num(SB[c][0], posinf=0.0, neginf=0.0)[m]
                    if "перемешанный" in arm:
                        B = B[np.random.default_rng(seed * 7 + 1).permutation(len(B))]
                    Xi = np.hstack([X[m], B])
                p = np.zeros(len(y))
                for f in range(5):
                    trn, te = fi != f, fi == f
                    if te.sum() == 0:
                        continue
                    p[te] = boost(Xi[trn], y[trn], Xi[te], seed * 10 + f)
                q = fit_apply(p, lo, hi, fi, np.ones(len(y)) / len(y))
                out[f"{seed}|{arm}|{c}"] = p.tolist()
                r[f"{c} пара"] = round(float(strae(y, q, y_true_upper=hi, y_true_lower=lo)), 4)
                r[f"{c} rho"] = round(float(spearmanr(y, p).statistic), 4)
            for tag in ("пара", "rho"):
                r[f"MACRO {tag}"] = round(float(np.mean([r[f"{c} {tag}"] for c in CYPS])), 4)
            table.append(r)
            print(f"  сид {seed} {arm:22s} пара {r['MACRO пара']:.4f} rho {r['MACRO rho']:.4f}  "
                  + " ".join(f"{c[3:]} {r[f'{c} rho']:.3f}" for c in CYPS)
                  + f"  ({time.time()-t0:.0f} с)", flush=True)

    df = pd.DataFrame(table)
    g = df.groupby("рука", sort=False)[["MACRO rho"] + [f"{c} rho" for c in CYPS]].mean()
    print("\n" + g.round(4).to_string())
    if "база" in g.index:
        print(f"\n{'рука':22s} " + "".join(f"{c[3:]+' Δ':>12s}" for c in CYPS))
        for arm in g.index:
            if arm == "база":
                continue
            line = f"{arm:22s} "
            for c in CYPS:
                d = g.loc[arm, f"{c} rho"] - g.loc["база", f"{c} rho"]
                mark = "*" if abs(d) > FLOOR[c] else " "
                line += f"{d:+11.4f}{mark}"
            print(line)
        print("  * --- больше СОБСТВЕННОГО пола фермента (пункт 165), не макро-пола")
    json.dump({"table": table, "preds": out}, open(a.out, "w"))
    print(f"\nсохранено: {a.out}")
    print("""
Как читать. Решают две строки на своих ферментах: «+сайт» на CYP1A2 против 0.0061 и на CYP3A4
против 0.0033. CYP2D6 --- проверка стенда: его углы обязаны что-то дать, иначе сравнивать не с
чем. CYP2C9 обязан дать около нуля, иначе стенд противоречит пункту 167.

«+сайт перемешанный» отделяет содержание колонок от их числа: ширина матрицы и маргиналы те
же, разорвана связь с молекулой.

«+форма целиком» --- самый дешёвый соперник. Если шестнадцать колонок формы всем без разбора
дают столько же, сколько целевые блоки, то целиться было незачем, и вывод в том, что форму
надо было просто добавить, а не строить под каждый центр.""")


if __name__ == "__main__":
    main()
