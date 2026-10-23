"""The union corpus for a multi-head fine-tune: our labels as the primary task, the qHTS panel as
auxiliary supervision, on ONE molecular graph set and ONE fold vector.

WHY THIS AND NOT THE FIVE THINGS ALREADY CLOSED. Every previous use of external data produced a
MEMBER -- predictions averaged into the composition -- and every one died over the ensemble: items
176, 182, 191, 194 and 356, the last being this very panel at -0.0092 as a seventh member with a
solo rank of 0.5728 against item 348's 0.60 screen. Item 356's own words put the closure "on the
ceiling rather than on the source". A multi-head fine-tune does not make a member. It makes a better
ENCODER, and a frozen head on a pretrained encoder is the only intervention in nine uploads that
moved the board's ORDERING (items 324, 344). That is the mechanism difference, and it is why this is
not the sixth attempt at a closed thing.

Item 341 sizes the prize off `briford`'s published shape -- four nets averaged over
challenge + ChEMBL + qHTS -- at >= +0.048 of rho on our folds. Item 343 closed the NAIVE form: one
head fine-tuned from a warm start scores 0.5613 against the frozen probe's 0.6261, because 6525
labels overfit 0.82 M encoder parameters. The auxiliary heads are what stops that, and the ratio is
the argument: ~68.5 k auxiliary labels against 6.5 k primary.

WHAT IS IN AND WHAT IS NOT. ChEMBL is NOT here. `data/ext_thirdhead.npz` holds 8004 rows of 247
FEATURES and no SMILES, and nothing local can invert them, so those rows cannot be turned into
molecular graphs. The union is therefore challenge + panel, two of the three sources item 341 names,
and saying so is better than implying a completeness we do not have.

THE PANEL'S SHAPE, measured rather than assumed. 54177 rows are (molecule, enzyme) pairs over FIVE
enzymes -- CYP2C19 is there and is kept, since more chemistry costs nothing and conflicts with
nothing -- covering 13126 unique canonical SMILES. `censored=1` means pAC50 is absent ENTIRELY (0 of
31671 have one) and those rows carry score ~0: they are the inactives. So the panel gives a binary
label on 54177 pairs and a continuous pAC50 on 22409 of them, left-truncated at 4.10 by the assay's
detection limit, with means 5.18--5.36 against our 4.10--4.96.

That truncation is why item 144 refused the panel for calibration and why items 363/364 found a
single source does not transfer ordering. It is harmless HERE: the auxiliary head learns "active, and
if so how potent", the primary heads keep our own scale, and the scale mismatch never touches a
number we submit.

THE LEAKAGE TRAP, which is the reason this is a script with controls and not a notebook cell. 182 of
the panel's molecules ARE our training molecules. If a panel row for a fold-3 molecule enters the
training set while fold 3 is held out, the encoder has seen that molecule's activity and the
out-of-fold estimate is worthless. So a shared molecule INHERITS its challenge fold, and the control
asserts it for all 182 rather than trusting the join.

    uv run python src/union.py --seed 0

Writes data/union_seed0.npz. Exit: 0 ran, 2 a control failed, 3 an input is missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES

import argparse
import json

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

from cypsplit import butina_folds, fold_digest

RDLogger.DisableLog("rdApp.*")

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
PANEL_ENZ = ["CYP1A2", "CYP2C9", "CYP2C19", "CYP2D6", "CYP3A4"]
GOLDEN = "2d93c19815e14261"


def canon(s):
    m = Chem.MolFromSmiles(s) if isinstance(s, str) else None
    return Chem.MolToSmiles(m) if m is not None else None


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    out = a.out or (D + "union_seed%d.npz" % a.seed)

    for f in ("rows.csv", "cyp-challenge-TRAIN_inhibition.csv",
              "cyp-challenge-TEST-BLINDED.csv", "ncgc/panel.csv"):
        if not _pl.Path(D + f).exists():
            print("нет %s%s" % (D, f))
            return 3

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    panel = pd.read_csv(D + "ncgc/panel.csv")

    fold, ncl = butina_folds(list(rows.SMILES), seed=a.seed)
    dg = fold_digest(fold)
    print("  сид %d: дайджест %s, кластеров %d" % (a.seed, dg, ncl), flush=True)
    if a.seed == 0 and dg != GOLDEN:
        print("ОТКАЗ: разбиение сдвинулось, сид 0 дал %s" % dg)
        return 2

    # ---- canonical keys -------------------------------------------------------------
    ch_can = [canon(s) for s in rows.SMILES]
    if any(c is None for c in ch_can):
        print("ОТКАЗ: %d обучающих SMILES не разобрались" % sum(c is None for c in ch_can))
        return 2
    te_can = set(filter(None, (canon(s) for s in te.SMILES)))
    panel = panel.copy()
    panel["can"] = [canon(s) for s in panel.smiles]
    panel = panel[panel.can.notna()]
    print("  панель: строк %d, уникальных молекул %d" % (len(panel), panel.can.nunique()),
          flush=True)

    # CONTROL: nothing from the panel may be a blind-test molecule. Computed from the test
    # SMILES we already hold -- no database was queried and none may be.
    leak = te_can & set(panel.can)
    print("      КОНТРОЛЬ утечки: панель ∩ слепой тест = %d из %d" % (len(leak), len(te)),
          flush=True)
    if leak:
        print("ОТКАЗ: панель содержит тестовые молекулы — обучать на них нельзя")
        return 2
    # ...and the same matcher must FIND the training overlap, or it proves nothing.
    ch_set = set(ch_can)
    shared = ch_set & set(panel.can)
    print("      КОНТРОЛЬ матчера: панель ∩ обучающие = %d (должно быть > 0)" % len(shared),
          flush=True)
    if not shared:
        print("ОТКАЗ: ноль пересечения с обучающими — сравнение канонических не работает")
        return 2

    # ---- the molecule table: challenge rows first, in rows.csv order ----------------
    ch_fold = {c: int(f) for c, f in zip(ch_can, fold)}
    smiles = list(ch_can)
    folds = [int(f) for f in fold]
    source = [0] * len(ch_can)                     # 0 = challenge, 1 = panel only
    seen = dict((c, i) for i, c in enumerate(ch_can))
    for c in panel.can.drop_duplicates():
        if c in seen:
            continue
        seen[c] = len(smiles)
        smiles.append(c)
        folds.append(-1)                           # panel-only: always trainable
        source.append(1)
    n = len(smiles)
    folds = np.asarray(folds, np.int64)
    source = np.asarray(source, np.int64)
    print("  объединение: %d молекул (%d челлендж + %d только панель)"
          % (n, int((source == 0).sum()), int((source == 1).sum())), flush=True)

    # CONTROL: a shared molecule must carry its CHALLENGE fold, not -1.
    bad = [c for c in shared if folds[seen[c]] != ch_fold[c]]
    print("      КОНТРОЛЬ наследования фолда: общих %d, с неверным фолдом %d"
          % (len(shared), len(bad)), flush=True)
    if bad:
        print("ОТКАЗ: общая молекула не унаследовала фолд — это утечка через фолд")
        return 2
    if (folds[source == 1] != -1).any():
        print("ОТКАЗ: молекула только из панели получила реальный фолд")
        return 2

    # ---- primary targets: our pIC50 and the band, masked per enzyme -----------------
    Y = np.full((n, len(CYPS)), np.nan, np.float32)
    LO = np.full((n, len(CYPS)), np.nan, np.float32)
    HI = np.full((n, len(CYPS)), np.nan, np.float32)
    for e, c in enumerate(CYPS):
        col = c + "_pIC50_direct_inhibition"
        for i, cc in enumerate(ch_can):
            Y[seen[cc], e] = tr[col].iloc[i]
            LO[seen[cc], e] = tr[col + "_conf_low"].iloc[i]
            HI[seen[cc], e] = tr[col + "_conf_high"].iloc[i]

    # ---- auxiliary targets: binary outcome and pAC50, five enzymes ------------------
    # `luc` flags a luciferase counter-screen artefact (472 rows of 54177). Masked rather
    # than kept: a known assay artefact is not a label, and dropping it is recorded here.
    use = panel[panel.luc == 0]

    # REPLICATES, and the first version of this file got the control wrong by ignoring them.
    # 3845 of 48685 (molecule, enzyme) pairs carry more than one row -- 225 molecules map
    # several CIDs onto one canonical structure -- so a cell-by-cell assignment silently
    # overwrites and the transfer count came out 3894 short. They are aggregated on purpose:
    #   pAC50   -> MEAN of the replicates, the standard summary of repeated measurement;
    #   outcome -> masked when the replicates CONFLICT (320 pairs of 3663 carry both Active
    #              and Inactive), because a conflict is missing information and "any active"
    #              would bias the label toward activity.
    # The replicate spread is reported rather than buried: its median is 0.400 of a log unit,
    # which is wider than our own labels' whole credible band (half-width 0.136 to 0.263), and
    # that is the quantity that should set how hard the auxiliary heads are weighted.
    det = use[use.outcome.isin(["Active", "Inactive"])].copy()
    det["act"] = (det.outcome == "Active").astype(float)
    gb = det.groupby(["can", "enzyme"])["act"]
    bin_mean, bin_nuniq = gb.mean(), gb.nunique()
    bin_ok = bin_mean[bin_nuniq == 1]                       # unanimous replicates only
    reg = use[use.pAC50.notna()].groupby(["can", "enzyme"])["pAC50"]
    reg_mean, reg_n = reg.mean(), reg.count()
    spread = (reg.max() - reg.min())[reg_n > 1]

    AB = np.full((n, len(PANEL_ENZ)), np.nan, np.float32)      # 1 active, 0 inactive
    AR = np.full((n, len(PANEL_ENZ)), np.nan, np.float32)      # pAC50 where measured
    ecol = {c: e for e, c in enumerate(PANEL_ENZ)}
    for (cc, en), v in bin_ok.items():
        AB[seen[cc], ecol[en]] = float(v)
    for (cc, en), v in reg_mean.items():
        AR[seen[cc], ecol[en]] = float(v)

    # CONTROL: the filled cells must equal the aggregated pairs, counted independently by
    # groupby rather than by the loop that wrote them. The reference is the number of PAIRS,
    # not of rows -- which is the correction the first version needed.
    got_bin, got_reg = int(np.isfinite(AB).sum()), int(np.isfinite(AR).sum())
    want_bin, want_reg = len(bin_ok), len(reg_mean)
    # Counted straight off the groupby rather than inferred: the first draft of this line
    # mixed two indexes and printed 874 / 0.8 per cent, which is not any quantity at all.
    allp = use.groupby(["can", "enzyme"]).size()
    print("      повторы: пар всего %d, с повторами %d (%.1f%%), строк в них %d"
          % (len(allp), int((allp > 1).sum()), 100.0 * float((allp > 1).mean()),
             int(allp[allp > 1].sum())), flush=True)
    print("      с определённым исходом %d пар, из них противоречивых %d (%.1f%%) --- замаскированы"
          % (len(bin_mean), int((bin_nuniq > 1).sum()),
             100.0 * float((bin_nuniq > 1).mean())), flush=True)
    print("      разброс pAC50 в повторах: пар %d, медиана %.3f, макс %.3f"
          % (len(spread), float(spread.median()) if len(spread) else float("nan"),
             float(spread.max()) if len(spread) else float("nan")), flush=True)
    print("      КОНТРОЛЬ переноса: бинарных %d против %d сведённых пар%s"
          % (got_bin, want_bin, "" if got_bin == want_bin else "  РАСХОЖДЕНИЕ"), flush=True)
    print("      КОНТРОЛЬ переноса: pAC50 %d против %d сведённых пар%s"
          % (got_reg, want_reg, "" if got_reg == want_reg else "  РАСХОЖДЕНИЕ"), flush=True)
    if got_bin != want_bin or got_reg != want_reg:
        print("ОТКАЗ: перенос метки потерял или размножил пары")
        return 2
    prim = int(np.isfinite(Y).sum())
    print("      основных меток %d, вспомогательных %d (отношение 1 : %.1f)"
          % (prim, got_bin + got_reg, (got_bin + got_reg) / prim), flush=True)

    np.savez_compressed(out, smiles=np.array(smiles, dtype=object), fold=folds,
                        source=source, Y=Y, LO=LO, HI=HI, AB=AB, AR=AR,
                        cyps=np.array(CYPS, dtype=object),
                        panel_enz=np.array(PANEL_ENZ, dtype=object),
                        digest=dg, seed=a.seed)
    print("\nзаписано: %s" % out)
    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump({"seed": a.seed, "digest": dg, "n": n,
               "n_challenge": int((source == 0).sum()),
               "n_panel_only": int((source == 1).sum()),
               "shared": len(shared), "leak_vs_test": len(leak),
               "primary_labels": prim, "aux_binary": got_bin, "aux_pac50": got_reg,
               "luc_dropped": int((panel.luc == 1).sum())},
              open(RES + "logs/union_seed%d.json" % a.seed, "w"), ensure_ascii=False, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
