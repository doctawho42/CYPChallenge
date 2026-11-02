"""Score the multi-head union fine-tune, in the PINNED environment.

The split of labour is the one `src/embed.py` established and item 343 followed: chemprop resolves
scikit-learn to 1.9.1 and CLAUDE.md records that 1.9.0 moves every prediction here by ~0.135 pIC50,
so the GPU box trains and writes a file, and every number that reaches the journal is computed here
under 1.3.2 by the same `score` items 324, 343, 345 and 350 used.

TWO ARMS, because the question is about the ENCODER and not about this model's own head.

  `голова`  -- the fine-tuned model's own out-of-fold predictions. Directly comparable to item
               343's 0.5613, which is what a NAIVE single-head fine-tune of this same checkpoint
               scored, and to the frozen probe's 0.6261.
  `зонд`    -- the shipped RidgeCV probe run on the fine-tuned encoder's embedding, nested: for
               fold f the embedding comes from the encoder that never saw fold f, and the probe is
               fitted on that same embedding over the other four folds. This is the device that
               measured 0.6261 on the FROZEN encoder, so it is the one that answers whether
               fine-tuning improved the representation or damaged it.

THE REFERENCE NUMBERS, all on these folds and already in the journal:

    замороженный зонд (пункт 324, 8 сидов)                 0.6261
    наивное одноголовое дообучение (пункт 343)             0.5613
    случайные веса той же формы (пункт 343, rand600)       0.4942
    слабейший член отгруженного состава                    0.5778
    отгруженный состав из шести                            0.6568
    экран пункта 348: член помогает выше ~0.60

THE CONTROL THAT MAKES THE NUMBERS READABLE. The five shipped members plus the FROZEN probe must
reproduce item 324's +0.0138 on this machine before any fine-tuned arm is scored. If that fails the
device has moved and nothing below compares to anything.

    uv run python verify/k141_ftunion.py --arms unionA,noaux,randaux

Reads data/ft/ft_union_<tag>.npz. Writes results/logs/k141_ftunion.json.
Exit: 0 ran, 2 a control failed, 3 no arm found.
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
from sklearn.linear_model import RidgeCV

from cypsplit import butina_folds, fold_digest
from submit import CYPS

import k106_probe as K
from k106_probe import score, errcorr, members, probe_oof, combine_add, load_emb

LOG = RES + "logs/k141_ftunion.json"
FLOOR = 0.0036
GOLDEN = "2d93c19815e14261"
REF = {"замороженный зонд (п.324)": 0.6261, "наивное дообучение (п.343)": 0.5613,
       "случайные веса (п.343)": 0.4942, "слабейший член": 0.5778,
       "состав из шести": 0.6568}


def nested_probe(emb, y, mask, fold):
    """RidgeCV on a PER-FOLD embedding, nested so the encoder never saw the rows it predicts.

    `emb[f]` is the fold-f encoder applied to all 4905 rows. For fold f the probe is fitted on
    that same embedding over the other four folds and predicts fold f, so the representation and
    the head are both out of fold. Standardisation is on the fitting rows only --- fitting the
    scaler on the union of a fold's halves is item 202's defect 2.
    """
    out = []
    for e in range(len(CYPS)):
        m = mask[:, e]
        p = np.full(int(m.sum()), np.nan)
        ym, fm = y[m, e], fold[m]
        for f in range(5):
            trn, te = fm != f, fm == f
            if te.sum() == 0:
                continue
            E = emb[f][m]
            mu, sd = E[trn].mean(0), E[trn].std(0)
            keep = sd > 1e-6
            Z = (E[:, keep] - mu[keep]) / sd[keep]
            p[te] = RidgeCV(alphas=K.ALPHAS).fit(Z[trn], ym[trn]).predict(Z[te])
        if not np.isfinite(p).all():
            raise SystemExit("%s: вложенный зонд оставил незаполненные строки" % CYPS[e])
        out.append(p)
    return out


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--arms", default="unionA")
    ap.add_argument("--seed", type=int, default=0)
    a = ap.parse_args(argv)
    tags = [t.strip() for t in a.arms.split(",") if t.strip()]

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
    print("  сид %d: дайджест %s, кластеров %d" % (a.seed, dg, ncl), flush=True)
    if a.seed == 0 and dg != GOLDEN:
        print("ОТКАЗ: разбиение сдвинулось, сид 0 дал %s" % dg)
        return 2

    # ---- the control: the shipped six must reproduce item 324 before anything else -------
    pdz = members(a.seed, XFULL, Y, M, fold, LO, HI)
    frozen = [probe_oof(load_emb(rows), Y, M, fold, e)[0] for e in range(len(CYPS))]
    five = [np.mean([Q[e] for k, Q in pdz if K._keep(e, k)], 0) for e in range(len(CYPS))]
    six = combine_add(pdz, frozen)
    s_five, s_six = score(five, Y, LO, HI, M, fold), score(six, Y, LO, HI, M, fold)
    s_frozen = score(frozen, Y, LO, HI, M, fold)
    print("      КОНТРОЛЬ состав: пять %.4f -> шесть %.4f, +%.4f (пункт 324: +0.0138)"
          % (s_five[0], s_six[0], s_six[0] - s_five[0]), flush=True)
    if s_six[0] - s_five[0] < 0.005:
        print("ОТКАЗ: замороженный зонд не воспроизводит выигрыш пункта 324")
        return 2
    print("      КОНТРОЛЬ: одиночный ранг ЗАМОРОЖЕННОГО зонда %.4f (пункт 324: 0.6261)"
          % s_frozen[0], flush=True)

    out = json.load(open(LOG)) if _pl.Path(LOG).exists() else {}
    out.setdefault("meta", {}).update(floor=FLOOR, refs=REF, frozen_here=s_frozen[0])
    res = out.setdefault("arms", {})

    for tag in tags:
        p = _pl.Path(D + "ft/ft_union_%s.npz" % tag)
        if not p.exists():
            print("      %-9s файла нет, пропуск" % tag, flush=True)
            continue
        f = np.load(p, allow_pickle=True)
        if str(f["digest"]) != dg:
            print("ОТКАЗ: %s собран на разбиении %s, а не %s" % (tag, str(f["digest"]), dg))
            return 2
        meta = json.loads(str(f["meta"]))
        P = f["oof"]
        if not np.isfinite(P).all():
            nm = int((~np.isfinite(P)).any(1).sum())
            print("ОТКАЗ: %s оставил %d строк незаполненными" % (tag, nm))
            return 2

        head = [P[M[:, e], e] for e in range(len(CYPS))]
        s_head = score(head, Y, LO, HI, M, fold)
        emb = [f["fold%d" % i].astype(np.float64) for i in range(5)]
        pr = nested_probe(emb, Y, M, fold)
        s_pr = score(pr, Y, LO, HI, M, fold)

        # as a further member on top of the shipped six, both ways
        sev_h = [np.mean([Q[e] for k, Q in pdz if K._keep(e, k)] + [frozen[e], head[e]], 0)
                 for e in range(len(CYPS))]
        sev_p = [np.mean([Q[e] for k, Q in pdz if K._keep(e, k)] + [frozen[e], pr[e]], 0)
                 for e in range(len(CYPS))]
        g_h = score(sev_h, Y, LO, HI, M, fold)
        g_p = score(sev_p, Y, LO, HI, M, fold)
        _pm, vm_h, _vs, _c = errcorr(head, pdz, Y, M)
        _pm2, vm_p, _vs2, _c2 = errcorr(pr, pdz, Y, M)

        res[tag] = {"aux": meta["aux"], "lam_bin": meta["lam_bin"], "lam_reg": meta["lam_reg"],
                    "epochs": [d["epochs"] for d in meta["log"]],
                    "head_solo": s_head, "probe_solo": s_pr,
                    "head_seventh": g_h[0] - s_six[0], "probe_seventh": g_p[0] - s_six[0],
                    "head_seventh_pair": g_h[1] - s_six[1],
                    "probe_seventh_pair": g_p[1] - s_six[1],
                    "head_corr_vs_five": float(np.mean(vm_h)),
                    "probe_corr_vs_five": float(np.mean(vm_p))}
        print("      %-9s голова %.4f | зонд %.4f | Δ7-м голова %+.4f зонд %+.4f | эпох %s"
              % (tag, s_head[0], s_pr[0], g_h[0] - s_six[0], g_p[0] - s_six[0],
                 "/".join(str(d["epochs"]) for d in meta["log"])), flush=True)

    if not res:
        print("ни одного арма не найдено")
        return 3

    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)

    print("\n" + "=" * 96)
    print("МНОГОГОЛОВОЕ ДООБУЧЕНИЕ НА ОБЪЕДИНЕНИИ, сид %d" % a.seed)
    print("=" * 96)
    print("  %-30s %9s %10s %10s" % ("ориентир", "одиночн.", "", ""))
    for k, v in REF.items():
        print("  %-30s %9.4f" % (k, v))
    print("  %-30s %9.4f   <- измерено здесь" % ("замороженный зонд, этот прогон", s_frozen[0]))
    print("  " + "-" * 90)
    print("  %-14s %4s %9s %9s | %10s %10s | %9s"
          % ("арм", "aux", "голова", "зонд", "Δ7-м гол.", "Δ7-м зонд", "corr зонда"))
    for tag, v in res.items():
        print("  %-14s %4s %9.4f %9.4f | %+10.4f %+10.4f | %9.3f"
              % (tag, "да" if v["aux"] else "нет", v["head_solo"][0], v["probe_solo"][0],
                 v["head_seventh"], v["probe_seventh"], v["probe_corr_vs_five"]))
    print("\n  пол макро-ранга %.4f; экран пункта 348 около 0.60" % FLOOR)
    print("  РЕШАЕТ: превосходит ли зонд на ДООБУЧЕННОМ энкодере замороженные 0.6261.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
