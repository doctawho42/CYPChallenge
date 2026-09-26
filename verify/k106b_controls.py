"""Three controls on verify/k106_probe.py's +0.0138, each built so it COULD refute it.

A gain of +0.0138 of macro rank is 3.8x this project's floor and larger than anything the
ensemble has bought since item 164. CLAUDE.md's rule is that a measurement is paired with a
control that must succeed AND that exercises the suspect clause. Three clauses are suspect:

  S  SHUFFLE  -- "the gain is the embedding" could instead be "adding a sixth vector to a
                 four-or-fewer-member mean helps whatever the vector is". Permuting the
                 embedding's ROWS keeps every marginal, every column statistic and the whole
                 RidgeCV path identical and destroys only the molecule correspondence. If the
                 gain survives a shuffle it was never about chemprop.
  L  LIBRARY  -- "a pretrained encoder" could instead be "a RidgeCV head with a wide alpha
                 grid, standardised, is just a good sixth member". Same head, same alphas,
                 same standardisation, OUR OWN DESC+MECH block instead of the embedding.
                 This is the arm item 316 would have run if it had run one.
  N  NEIGHBOUR-- the corpus the checkpoint was pretrained on cannot be audited (the 750 blind
                 compounds must never be looked up externally, and neither may our training
                 set). If the encoder MEMORISED per-molecule CYP activity, the probe would not
                 need a near neighbour in the training folds and our models would; the probe's
                 advantage would then be FLAT or GROW as similarity to the training folds
                 falls. If the encoder learned transferable chemistry, the advantage should
                 behave like a model's and shrink. This is suggestive, not decisive -- stated
                 as such -- but it is the only handle on that risk that queries no database.

Everything is computed from local files and from the test/training SMILES only.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from rdkit import Chem, DataStructs, RDLogger
from rdkit.Chem import rdFingerprintGenerator

from cypsplit import butina_folds, fold_digest
import submit as S
from submit import CYPS, _keep, DESC_MECH
from k106_probe import (probe_oof, members, load_emb, combine_keep, combine_add, score,
                        MEMBERS, ALPHAS)

RDLogger.DisableLog("rdApp.*")
LOG = RES + "logs/k106b_controls.json"


def setup():
    rows = pd.read_csv(D + "rows.csv")
    z = np.load(D + "feats.npz")
    XFULL = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)
    S.LO, S.HI = LO, HI
    return rows, XFULL, y, LO, HI, mask


def probe_all(E, y, mask, fold):
    return [probe_oof(E, y, mask, fold, e)[0] for e in range(len(CYPS))]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    rows, XFULL, y, LO, HI, mask = setup()
    E = load_emb(rows)
    DM = np.nan_to_num(XFULL[:, -DESC_MECH:].astype(np.float64), posinf=0.0, neginf=0.0)
    print(f"  DESC+MECH для контроля L: {DM.shape}", flush=True)

    out = json.load(open(LOG)) if _pl.Path(LOG).exists() else {}
    out.setdefault("seeds", {})

    for s in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        print(f"\n  сид {s}: дайджест {fold_digest(fold)}", flush=True)
        pdz = members(s, XFULL, y, mask, fold, LO, HI)
        A = score(combine_keep(pdz), y, LO, HI, mask, fold)

        pr = probe_all(E, y, mask, fold)
        Dtrue = score(combine_add(pdz, pr), y, LO, HI, mask, fold)

        # --- S: the same embedding, rows permuted -------------------------------------------
        perm = np.random.default_rng(1000 + s).permutation(len(E))
        # Control on the control: the permutation must actually move rows.
        if int((perm == np.arange(len(E))).sum()) > len(E) // 100:
            raise SystemExit("перестановка почти тождественна")
        prS = probe_all(E[perm], y, mask, fold)
        DS = score(combine_add(pdz, prS), y, LO, HI, mask, fold)
        soloS = score(prS, y, LO, HI, mask, fold)

        # --- L: identical head, our own DESC+MECH -------------------------------------------
        prL = probe_all(DM, y, mask, fold)
        DL = score(combine_add(pdz, prL), y, LO, HI, mask, fold)
        soloL = score(prL, y, LO, HI, mask, fold)

        rec = {"A": A, "D": Dtrue, "D_shuffle": DS, "D_library": DL,
               "solo_probe": score(pr, y, LO, HI, mask, fold),
               "solo_shuffle": soloS, "solo_library": soloL}
        out["seeds"][str(s)] = rec
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
        print(f"    A {A[0]:.4f} | D {Dtrue[0]:+.4f}->{Dtrue[0]-A[0]:+.4f} | "
              f"S {DS[0]-A[0]:+.4f} | L {DL[0]-A[0]:+.4f}", flush=True)

    # --- N: similarity stratification, seed 0 only ------------------------------------------
    if 0 in seeds:
        neighbour_control(rows, XFULL, y, LO, HI, mask, E, out)

    report(out)


def neighbour_control(rows, XFULL, y, LO, HI, mask, E, out):
    """Per-compound max Tanimoto to compounds in the OTHER folds -- i.e. to the rows the model
    that predicted it was actually trained on -- then the probe's advantage over GP by quartile.

    Same construction as verify/f12_cvhard.py line-for-line in spirit: Morgan r=2, 2048 bits,
    maximum over the training folds, which is what a held-out prediction can lean on.
    """
    print("\n  === контроль N: преимущество зонда против сходства с обучающими фолдами ===",
          flush=True)
    fold, _ = butina_folds(list(rows.SMILES), seed=0)
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    fps = [gen.GetFingerprint(Chem.MolFromSmiles(s)) for s in rows.SMILES]
    nn = np.empty(len(fps))
    for i in range(len(fps)):
        other = np.where(fold != fold[i])[0]
        nn[i] = max(DataStructs.BulkTanimotoSimilarity(fps[i], [fps[j] for j in other]))
    # Control: the statistic must be capable of separating. A compound compared against its OWN
    # fold-mates including itself would read 1.000 everywhere; it does not.
    print(f"    сходство с ближайшим из ДРУГИХ фолдов: медиана {np.median(nn):.3f}, "
          f"доля 1.000 {np.mean(nn > 0.999):.4f} (контроль: должно быть мало)", flush=True)

    pdz = members(0, XFULL, y, mask, fold, LO, HI)
    pr = probe_all(E, y, mask, fold)
    GP = dict(pdz)["GP"]
    ship = [np.mean([P[e] for k, P in pdz if _keep(e, k)], axis=0) for e in range(len(CYPS))]

    edges = np.percentile(nn, [0, 25, 50, 75, 100])
    print(f"\n    {'квартиль сходства':>22s} {'n':>6s} " +
          " ".join(f"{c:>17s}" for c in CYPS), flush=True)
    print(f"    {'':>22s} {'':>6s} " + " ".join(f"{'зонд/GP/состав':>17s}" for _ in CYPS),
          flush=True)
    rowsout = []
    for q in range(4):
        lab = f"{edges[q]:.2f}-{edges[q+1]:.2f}"
        cells, ns, adv = [], [], []
        for e in range(len(CYPS)):
            m = mask[:, e]
            sel = (nn[m] >= edges[q]) & (nn[m] <= edges[q + 1] if q == 3 else nn[m] < edges[q + 1])
            yy = y[m, e][sel]
            if sel.sum() < 30 or len(np.unique(yy)) < 5:
                cells.append(f"{'--':>17s}")
                continue
            rp = float(spearmanr(yy, pr[e][sel]).statistic)
            rg = float(spearmanr(yy, GP[e][sel]).statistic)
            rs = float(spearmanr(yy, ship[e][sel]).statistic)
            cells.append(f"{rp:5.3f}/{rg:5.3f}/{rs:5.3f}")
            ns.append(int(sel.sum()))
            adv.append(rp - rs)
        print(f"    {lab:>22s} {int(np.mean(ns)) if ns else 0:6d} " + " ".join(cells), flush=True)
        rowsout.append({"bin": lab, "adv_vs_shipped": adv})
    out["neighbour"] = {"edges": [float(x) for x in edges], "bins": rowsout}
    json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
    print("\n    Читать так: если преимущество зонда над составом РАСТЁТ, когда сходство"
          "\n    падает, это подпись запоминания (зонду не нужен близкий сосед, нашим моделям"
          "\n    нужен). Если падает --- зонд ведёт себя как модель. Свидетельство"
          "\n    наводящее, не решающее.", flush=True)


def report(out):
    ks = sorted([k for k in out["seeds"] if "D_shuffle" in out["seeds"][k]], key=int)
    if not ks:
        return
    print(f"\n\n  === КОНТРОЛИ, сидов {len(ks)} ({', '.join(ks)}) ===", flush=True)
    print(f"\n  {'':>40s} {'соло ранг':>10s} {'Δ над A':>10s} {'знак':>7s}", flush=True)
    for arm, solo, lab in [("D", "solo_probe", "зонд chemprop_medium (замер)"),
                           ("D_shuffle", "solo_shuffle", "S: тот же зонд, строки переставлены"),
                           ("D_library", "solo_library", "L: та же голова на DESC+MECH")]:
        d = np.array([out["seeds"][k][arm][0] - out["seeds"][k]["A"][0] for k in ks])
        r = np.array([out["seeds"][k][solo][0] for k in ks])
        print(f"  {lab:>40s} {r.mean():10.4f} {d.mean():+10.4f} "
              f"{int((d > 0).sum()):4d}/{len(d):<2d}", flush=True)
    print(f"\n  Контроль S должен ОБНУЛИТЬ или обратить прирост; контроль L отделяет"
          f"\n  «предобученный кодировщик» от «ещё одна гребневая голова».", flush=True)


if __name__ == "__main__":
    main()
