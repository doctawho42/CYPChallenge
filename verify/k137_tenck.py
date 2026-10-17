"""The ten organiser checkpoints items 349 and 350 left unmeasured, each as a member and then as a
STACK --- because item 350 predicted each would fail and item 346 showed that failures can add.

THE TWO PREDICTIONS THIS TESTS, both already in writing. Item 350: "The refined criterion predicts
the ten remaining fail, now on FIVE measured descriptor cases spanning 2D, 3D, topological and shape
families and target widths from 60 to 2048." And item 346's lesson, which cuts the other way: three
members each consistently positive and consistently UNDER the 0.0036 floor --- item 335's gbm head
at +0.0023, item 342's graph net at +0.0006, item 345's chemeleon at +0.0033 --- summed to +0.0047
and cleared it, because "the project's rule judges a candidate against the composition, never
candidates against one another, so a threshold that is right for one attempt discards a SET of
attempts whose sum clears it."

So the question item 350 did not ask is not whether each of the ten helps. It is whether the ten
TOGETHER do. Ten sub-floor positives would be worth more than any one of them, and ten negatives
settle the criterion on its own prediction.

    atompair  ecfp  erg  jazzy  osmordred  osmordred_pca80/90/95  pubchem
    and minimol, which item 349 ALREADY measured at -0.0013 at sign 0/4 --- kept in the list
    deliberately, as the reproduction control for this whole sweep.

WHAT IS MEASURED. Item 324's device throughout: a RidgeCV probe on the frozen embedding, our five
Butina folds, standardisation and the dead-column filter on each fold's training rows, the
prediction entering as a sixth member beside the five `_keep` admits, scored after the per-fold
affine pair. Then every positive candidate added together, and the whole set added together, since
item 346's result says the sum is the question.

THE CONTROLS. `minimol` must reproduce item 349's -0.0013 at 0/4, or the device here is not the
device there. Seed 0's Butina digest must be the golden `2d93c19815e14261`. Each embedding is built
by `src/embed.py` out of process and carries no `molecule_name`, so each gets item 345's alignment
control: the probe refitted on the array ROLLED BY ONE ROW must collapse to near zero, which is
asserted rather than printed. And no embedding may have a dead column count that differs from what
`src/embed.py` reported when it wrote the file.

    uv run python verify/k137_tenck.py --seeds 0,1

Writes results/logs/k137_tenck.json.
Exit: 0 ran, 2 a control failed, 3 no embedding is present.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import glob
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from cypsplit import butina_folds, fold_digest
from submit import CYPS, _keep

import k106_probe as K
from k106_probe import score, probe_oof

LOG = RES + "logs/k137_tenck.json"
GOLDEN = "2d93c19815e14261"
MEMBERS = ("поферментно", "пул", "GP", "гребневая", "ствол")
FLOOR, FLOOR_PAIR = 0.0036, 0.007
ITEM349_MINIMOL = -0.0013
ROLL_MAX = 0.05            # the rolled-row control must stay under this in |rho|


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    embs = sorted(glob.glob(D + "emb_ck_*.npz"))
    if not embs:
        print("нет ни одного data/emb_ck_*.npz", file=_sys.stderr)
        return 3
    tags = [_pl.Path(p).stem.replace("emb_ck_", "") for p in embs]
    print("кандидатов на диске: %d --- %s" % (len(tags), ", ".join(tags)), flush=True)

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)
    K.S.LO, K.S.HI = LO, HI

    out = {"candidates": tags, "seeds": seeds, "per_seed": {}}
    for s in seeds:
        fold, _n = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        if s == 0 and dg != GOLDEN:
            print("ОТКАЗ: сид 0 дал %s" % dg, file=_sys.stderr)
            return 2
        mp = _pl.Path(RES + f"preds/members_seed{s}.json")
        if not mp.exists():
            print("  сид %d: нет %s --- пропуск" % (s, mp.name))
            continue
        mem = json.load(open(mp))
        pdz = [(k, [np.asarray(mem[k][e], float) for e in range(len(CYPS))]) for k in MEMBERS]
        kept = [[Q[e] for k, Q in pdz if _keep(e, k)] for e in range(len(CYPS))]
        five = score([np.mean(kept[e], 0) for e in range(len(CYPS))], Y, LO, HI, M, fold)
        cell = {"digest": dg, "пять": five, "arms": {}}
        print("\n  сид %d: пять = ранг %.4f, пара %.4f" % (s, five[0], five[1]), flush=True)
        P = {}
        for p_, tag in zip(embs, tags):
            E = np.load(p_, allow_pickle=True)["train"].astype(np.float64)
            if E.shape[0] != len(rows):
                print("ОТКАЗ: %s даёт %d строк против %d" % (tag, E.shape[0], len(rows)),
                      file=_sys.stderr)
                return 2
            pe = [probe_oof(E, Y, M, fold, e)[0] for e in range(len(CYPS))]
            # item 345's alignment control, asserted: nameless file, so the order is exercised
            rc = probe_oof(np.roll(E, 1, axis=0), Y, M, fold, 0)[0]
            m0 = M[:, 0]
            rho_roll = float(spearmanr(Y[m0, 0], rc).statistic)
            rho_ok = float(spearmanr(Y[m0, 0], pe[0]).statistic)
            if abs(rho_roll) > ROLL_MAX or rho_ok <= abs(rho_roll):
                print("КОНТРОЛЬ ВЫРАВНИВАНИЯ ПРОВАЛЕН для %s: сдвинутый ро %.4f против %.4f"
                      % (tag, rho_roll, rho_ok), file=_sys.stderr)
                return 2
            P[tag] = pe
            six = score([np.mean(kept[e] + [pe[e]], 0) for e in range(len(CYPS))],
                        Y, LO, HI, M, fold)
            solo = score(pe, Y, LO, HI, M, fold)
            cell["arms"][tag] = {"six": six, "solo": solo, "d_rank": six[0] - five[0],
                                 "d_pair": six[1] - five[1], "roll_rho": rho_roll}
            print("    %-18s шестым %+.4f ранга, %+.4f пары | одиночн. %.4f | сдвиг-контроль %.4f"
                  % (tag, six[0] - five[0], six[1] - five[1], solo[0], rho_roll), flush=True)
        # --- the stacks item 346's lesson asks for ---
        pos = [t for t in tags if cell["arms"][t]["d_rank"] > 0]
        for nm, sub in (("все %d" % len(tags), tags), ("положительные (%d)" % len(pos), pos)):
            if not sub:
                continue
            Q = [np.mean(kept[e] + [P[t][e] for t in sub], 0) for e in range(len(CYPS))]
            sc = score(Q, Y, LO, HI, M, fold)
            cell["arms"]["СТЕК " + nm] = {"six": sc, "solo": None,
                                          "d_rank": sc[0] - five[0], "d_pair": sc[1] - five[1],
                                          "members": sub}
            print("    СТЕК %-13s      %+.4f ранга, %+.4f пары   (%s)"
                  % (nm, sc[0] - five[0], sc[1] - five[1], ", ".join(sub)), flush=True)
        out["per_seed"][str(s)] = cell
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)

    have = [s for s in seeds if str(s) in out["per_seed"]]
    if not have:
        print("ни одного сида", file=_sys.stderr)
        return 3
    names = sorted({k for s in have for k in out["per_seed"][str(s)]["arms"]},
                   key=lambda k: (k.startswith("СТЕК"), k))
    print("\n" + "=" * 92)
    print("%-20s %9s %8s %7s %9s %9s" % ("кандидат", "Δранг", "sd", "знак", "Δпара", "одиночн."))
    print("=" * 92)
    out["summary"] = {}
    for n in names:
        have_n = [s for s in have if n in out["per_seed"][str(s)]["arms"]]
        v = np.asarray([out["per_seed"][str(s)]["arms"][n]["d_rank"] for s in have_n], float)
        pr = np.asarray([out["per_seed"][str(s)]["arms"][n]["d_pair"] for s in have_n], float)
        so = [out["per_seed"][str(s)]["arms"][n]["solo"] for s in have_n]
        som = float(np.mean([x[0] for x in so])) if so[0] is not None else float("nan")
        sd = float(np.std(v, ddof=1)) if len(v) > 1 else float("nan")
        out["summary"][n] = {"d_rank": float(v.mean()), "sd": sd, "d_pair": float(pr.mean()),
                             "solo": som, "sign": [int((v > 0).sum()), len(v)]}
        print("%-20s %+9.4f %8.4f %5d/%d %+9.4f %9s"
              % (n, v.mean(), sd, int((v > 0).sum()), len(v), pr.mean(),
                 "—" if np.isnan(som) else "%.4f" % som))
    print("\nполы: макро ранг %.4f, ST-RAE %.3f" % (FLOOR, FLOOR_PAIR))
    if "minimol" in out["summary"]:
        mm = out["summary"]["minimol"]["d_rank"]
        print("контроль воспроизведения: minimol %+.4f против %+.4f пункта 349 (|Δ| %.4f)"
              % (mm, ITEM349_MINIMOL, abs(mm - ITEM349_MINIMOL)))
        if abs(mm - ITEM349_MINIMOL) > 0.004:
            print("КОНТРОЛЬ ПРОВАЛЕН: устройство не воспроизводит пункт 349", file=_sys.stderr)
            return 2
    json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
    print("\nзапись: %s" % LOG)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
