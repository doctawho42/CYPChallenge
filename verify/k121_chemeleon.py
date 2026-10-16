"""The SECOND frozen checkpoint as a seventh member --- rejected once, on an instrument since
measured non-predictive.

WHY IT IS BEING REOPENED, and the reason is in our own board history rather than in an argument.
Across nine uploads our board Spearman has taken exactly two values: 0.6935 for the five submissions
before 25 September and 0.7062 for the four after. ST-RAE moved 0.8149 -> 0.5357 over the same
period, which is 103rd place to 53rd. **Ordering moved once, and the thing that moved it was the
frozen-encoder probe of items 324/325.** Item 327 recorded the transfer: the out-of-fold rank gain
reached the board at 0.92, against 0.52 and 0.437 for items 314 and 317, which were post-processing.
So the apparatus measures rank correctly and we moved rank once in nine tries.

The probe's device --- a frozen public encoder, a RidgeCV head, the prediction entering the mean as
its own member --- is therefore the only lever in this project with a demonstrated effect on the
quantity the board actually ranks by. And there is a second checkpoint sitting in the same public
repository that this device was never run on.

HOW IT WAS REJECTED. Item 324, verbatim: "We extracted from `chemprop_medium.pt` only --- 7.7 MB,
sha256 `b8137a0d...` --- because the larger checkpoint's residuals correlate 0.91 with it on two
enzymes." That is the whole record. A grep of `verify/README.md` for `chemprop_chemeleon` returns
ZERO against a control of five for `chemprop_medium`: it was never fitted, never scored, never
admitted as a member. It was refused on a residual correlation measured on two of four enzymes ---
by the same error-correlation screen item 335 then measured as non-predictive, having passed four
candidates whose seventh-member outcomes span -0.0048 to +0.0138 in no particular order. Its own
author, whose file it is, calls it "a genuinely decorrelated second cluster from (2), not a
duplicate" and uses both.

WHAT IS MEASURED HERE. `probe_oof` from `k106_probe`, unchanged --- the same `ALPHAS`, the same
per-fold standardisation and dead-column filter --- on the new embedding, then `score` and `errcorr`
exactly as items 324, 335, 341 and 342 used them. The device is identical to the shipped probe's so
the comparison is about the checkpoint and nothing else.

Three references, all on these folds and already in the journal:

    зонд на chemprop_medium (отгружен)       0.6261 одиночного ранга, +0.0138 седьмым
    состав из шести (с зондом)               0.6568
    полоса корреляции ошибок пяти членов     0.966--0.982, сам зонд 0.901

The deciding number is the gain as a SEVENTH member beside the five and the existing probe, because
item 335 established that the screen does not predict it and only the direct measurement does. The
rejection figure is also reproduced properly: the two probes' error correlation on all four enzymes,
not two.

    uv run python verify/k121_chemeleon.py

Needs data/emb_chemprop_chemeleon.npz, written out of process by
`src/embed_medium.py --ckpt <...>/chemprop_chemeleon.pt --rows data/rows.csv`, whose own controls
reported `bn=Identity`, both extraction routes identical at 0.0, and order preserved under a
1024-molecule permutation.

Writes results/logs/k121_chemeleon.json. Exit: 0 ran, 2 a control failed, 3 the embedding is absent.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from cypsplit import butina_folds, fold_digest
from submit import CYPS

import k106_probe as K
from k106_probe import score, errcorr, members, probe_oof, combine_add, load_emb

LOG = RES + "logs/k121_chemeleon.json"
NEW = D + "emb_chemprop_chemeleon.npz"
FLOOR = 0.0036
GOLDEN = "2d93c19815e14261"
# Item 324 and item 335, measured on these folds.
PROBE_SOLO, PROBE_GAIN, PROBE_CORR = 0.6261, 0.0138, 0.901
BAND = (0.966, 0.982)


def load_any(path, rows):
    """`k106_probe.load_emb`'s controls, for an arbitrary path. Returns (E, named).

    The order test is shown to be CAPABLE of failing by repeating it on the same names rolled by
    one position: a comparison that cannot fail is not a control. The shape is checked against the
    table rather than trusted, because the test-side file in this family stores its array under the
    key `train` as well and a name that lies cannot be caught by reading the name.

    NOT EVERY EMBEDDING CARRIES NAMES. `src/embed.py` writes `emb_chemprop_rdkit2d.npz` by
    embedding `list(rows.SMILES)` in order and checking only the length, and the three
    language-model files do the same, so those four store `train` with no `molecule_name`. For
    them this function CANNOT check the order and says so by returning named=False; the caller
    must then run an alignment control that can fail. Returning a quiet True here would be the
    exact defect this file is otherwise built to avoid.
    """
    z = np.load(path, allow_pickle=True)
    E = z["train"].astype(np.float64)
    r = rows.Molecule_Name.to_numpy(object)
    if E.shape[0] != len(r):
        print("ОТКАЗ: %d строк вложения против %d в rows.csv" % (E.shape[0], len(r)))
        raise SystemExit(2)
    if not np.isfinite(E).all():
        print("ОТКАЗ: во вложении есть не-конечные значения")
        raise SystemExit(2)
    if "molecule_name" not in z.files:
        print("      вложение %s, имён в файле НЕТ --- порядок проверяется выравниванием ниже"
              % (E.shape,), flush=True)
        return E, False
    nm = np.asarray(z["molecule_name"], object)
    ok = bool((nm == r).all())
    rolled = bool((np.roll(nm, 1) == r).all())
    print("      вложение %s, порядок совпал: %s; сдвинутый контроль (ожидается False): %s"
          % (E.shape, ok, rolled), flush=True)
    if not ok or rolled:
        print("ОТКАЗ: порядок вложения не совпал с rows.csv (или контроль не сработал)")
        raise SystemExit(2)
    return E, True


def align_control(E, Y, M, fold, e=0):
    """For a nameless embedding: the probe must collapse when the rows are rolled by one.

    This is the control the name match would have been. It exercises the ORDER clause directly ---
    a misaligned feature matrix carries no signal about the labels --- and it must succeed, since
    the aligned array is known to score well. A control that only proved `probe_oof` runs would
    say nothing about alignment.
    """
    # `score` always loops over all four enzymes, so one enzyme is ranked directly here.
    yy = Y[M[:, e], e]
    a = float(spearmanr(yy, probe_oof(E, Y, M, fold, e)[0]).statistic)
    b = float(spearmanr(yy, probe_oof(np.roll(E, 1, axis=0), Y, M, fold, e)[0]).statistic)
    print("      ВЫРАВНИВАНИЕ (%s): как есть %.4f, сдвиг на строку %.4f -> %s"
          % (CYPS[e], a, b, "ПОДТВЕРЖДЕНО" if a - b > 0.05 else "НЕ ПОДТВЕРЖДЕНО"), flush=True)
    return a - b > 0.05


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seed", type=int, default=0)
    # Any frozen checkpoint can go through this device. The defaults reproduce the chemeleon
    # run bit for bit, so adding these flags cannot move a committed number; item 345's
    # `rdkit2d` negative is this same code with --emb/--tag/--out, which is what makes that
    # figure reproducible rather than a one-off heredoc.
    ap.add_argument("--emb", default=NEW, help="вложение-кандидат, .npz с train")
    ap.add_argument("--tag", default="chemeleon", help="как называть кандидата в выводе")
    ap.add_argument("--out", default=None, help="куда писать (по умолчанию %s)" % LOG)
    a = ap.parse_args(argv)
    dest = a.out or LOG

    if not _pl.Path(a.emb).exists():
        print("нет %s --- сначала извлеките вложение в окружении энкодера" % a.emb)
        return 3

    rows = pd.read_csv(D + "rows.csv")
    z = np.load(D + "feats.npz")
    XFULL = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)
    K.S.LO, K.S.HI = LO, HI

    fold, ncl = butina_folds(list(rows.SMILES), seed=a.seed)
    dg = fold_digest(fold)
    print("\n  сид %d: дайджест %s, кластеров %d" % (a.seed, dg, ncl), flush=True)
    if a.seed == 0 and dg != GOLDEN:
        print("ОТКАЗ: разбиение сдвинулось, сид 0 дал %s" % dg)
        return 2

    Eold = load_emb(rows)                 # chemprop_medium, the shipped probe's embedding
    Enew, named = load_any(a.emb, rows)   # the candidate; chemprop_chemeleon by default
    print("      столбцов: medium %d, %s %d"
          % (Eold.shape[1], a.tag, Enew.shape[1]), flush=True)
    if not named and not align_control(Enew, Y, M, fold):
        print("ОТКАЗ: выравнивание безымянного вложения не подтверждено")
        return 2

    pdz = members(a.seed, XFULL, Y, M, fold, LO, HI)
    old = [probe_oof(Eold, Y, M, fold, e)[0] for e in range(len(CYPS))]
    new = [probe_oof(Enew, Y, M, fold, e)[0] for e in range(len(CYPS))]

    # ---- controls that must succeed -------------------------------------------------
    s_old = score(old, Y, LO, HI, M, fold)
    ref = json.load(open(RES + "logs/k106_probe.json"))["seeds"].get(str(a.seed))
    if ref and abs(s_old[0] - ref["solo_probe"][0]) > 1e-6:
        print("ОТКАЗ: отгруженный зонд не воспроизводится: %.4f против %.4f"
              % (s_old[0], ref["solo_probe"][0]))
        return 2
    print("      КОНТРОЛЬ зонд воспроизведён: %.4f" % s_old[0], flush=True)

    five = [np.mean([Q[e] for k, Q in pdz if K._keep(e, k)], 0) for e in range(len(CYPS))]
    six = combine_add(pdz, old)
    s_five, s_six = score(five, Y, LO, HI, M, fold), score(six, Y, LO, HI, M, fold)
    print("      КОНТРОЛЬ состав: пять %.4f -> шесть %.4f, +%.4f (пункт 324: +%.4f)"
          % (s_five[0], s_six[0], s_six[0] - s_five[0], PROBE_GAIN), flush=True)
    if s_six[0] - s_five[0] < 0.005:
        print("ОТКАЗ: добавление зонда не воспроизводит выигрыш пункта 324")
        return 2
    # The two embeddings must not be the same array under two names.
    dmin = min(float(np.min(np.abs(new[e] - old[e]))) for e in range(len(CYPS)))
    if dmin < 1e-12:
        print("ОТКАЗ: предсказания двух зондов совпадают --- это, вероятно, один файл")
        return 2

    # ---- the measurement ------------------------------------------------------------
    s_new = score(new, Y, LO, HI, M, fold)
    _pm, vm_new, _vs, cal = errcorr(new, pdz, Y, M)
    _pm2, vm_old, _vs2, _c2 = errcorr(old, pdz, Y, M)

    # The figure item 324 rejected it on, reproduced on ALL FOUR enzymes rather than two.
    pair_corr = [float(np.corrcoef(new[e] - Y[M[:, e], e], old[e] - Y[M[:, e], e])[0, 1])
                 for e in range(len(CYPS))]

    seven = [np.mean([Q[e] for k, Q in pdz if K._keep(e, k)] + [old[e], new[e]], 0)
             for e in range(len(CYPS))]
    s7 = score(seven, Y, LO, HI, M, fold)
    # And the swap: the new checkpoint INSTEAD of the shipped one, which is the other question
    # a second encoder raises and costs nothing to ask here.
    swap = combine_add(pdz, new)
    s_swap = score(swap, Y, LO, HI, M, fold)

    out = json.load(open(dest)) if _pl.Path(dest).exists() else {}
    out.setdefault("meta", {}).update(floor=FLOOR, probe_solo=PROBE_SOLO, band=list(BAND),
                                      emb=str(a.emb), tag=a.tag)
    out.setdefault("seeds", {})[str(a.seed)] = {
        "solo_new": s_new, "solo_old": s_old, "five": s_five, "six": s_six,
        "seven": s7, "swap": s_swap,
        "gain_seventh": s7[0] - s_six[0], "gain_pair_seventh": s7[1] - s_six[1],
        "gain_swap": s_swap[0] - s_six[0],
        "errcorr_new_vs_five": vm_new, "errcorr_old_vs_five": vm_old,
        "errcorr_new_vs_old": pair_corr,
        "member_band": {k: float(np.mean(v)) for k, v in cal.items()},
        "cols": int(Enew.shape[1]),
    }
    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump(out, open(dest, "w"), ensure_ascii=False, indent=1)

    print("\n" + "=" * 94)
    print("ЗАМОРОЖЕННЫЙ ЧЕКПОЙНТ %s, сид %d" % (a.tag, a.seed))
    print("=" * 94)
    print("  %-30s %9s %9s %9s" % ("", "одиночн.", "пара", "corr к пяти"))
    print("  %-30s %9.4f %9.4f %9.3f" % ("зонд medium (отгружен)", s_old[0], s_old[1],
                                         float(np.mean(vm_old))))
    print("  %-30s %9.4f %9.4f %9.3f" % ("зонд %s (новый)" % a.tag, s_new[0], s_new[1],
                                         float(np.mean(vm_new))))
    print("  %-30s %9.4f %9s %9s" % ("полоса пяти членов", float(np.mean(
        [np.mean(v) for v in cal.values()])), "-", "-"))
    print()
    print("  корреляция ОШИБОК двух зондов по ферментам: %s   среднее %.3f"
          % (" ".join("%.3f" % v for v in pair_corr), float(np.mean(pair_corr))))
    print("  (пункт 324 отверг его по 0.91 на ДВУХ ферментах; здесь все четыре)")
    print()
    print("  состав: пять %.4f -> шесть %.4f -> СЕМЬ %.4f  (%+.4f ранга, пара %+.4f)"
          % (s_five[0], s_six[0], s7[0], s7[0] - s_six[0], s7[1] - s_six[1]))
    print("  замена: вместо medium взять %s -> %.4f  (%+.4f к шести)"
          % (a.tag, s_swap[0], s_swap[0] - s_six[0]))
    print("  пол %.4f" % FLOOR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
