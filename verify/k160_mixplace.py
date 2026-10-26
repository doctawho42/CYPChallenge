"""Place CYP2D6's predictions as a two-component MIXTURE -- the top fraction at the actives' level,
the rest where the wide band makes them free -- and measure it out of fold against the affine pair
the project ships, which is the thing item 316 says will win.

Why it is worth measuring despite that. Item 384 found that the credible band is seven to ten times
wider for low-activity compounds than for actives, so ST-RAE is decided almost entirely by the
actives, and the affine family cannot exploit that: it moves everything together. A mixture
placement can put the actives precisely and dump the rest anywhere low. Item 316 measured monotone
piecewise-linear placement as 0.0011 to 0.0014 WORSE than the pair out of fold at every bin count
from 8 to 128, sign 1/4, and item 31 closed out-of-fold isotonic before it -- so this arm is inside
a family that has lost twice, and the ONLY thing that distinguishes it is item 357's scoping: that
floor was computed on OUR label distribution, "where there is no mismatch by construction", and says
nothing about correcting a mismatch with the test.

So the honest form of the question is: does a mixture placement beat the pair when the mixture is
handed the TRUE active fraction and the TRUE pools, out of fold? If it loses even then, no derived
fraction can rescue it and the candidate is dead before it is uploaded.

CYP2D6 is the target because item 316 measured where the deficit is: of the 0.0721 of metric that
separates us from the best team at our own ordering, "CYP2D6 alone is 0.0449 -- sixty-two per cent",
and CYP2D6 is also the one isoform whose rank TRANSFERS WORSE to the board than our CV suggests
(-0.0424 against +0.1642 and +0.1016 on CYP1A2 and CYP2C9).

THE MEASUREMENT. Five cached members per seed, composed by the unweighted mean the project ships,
four seeds. The baseline is `submit.fit_shrinkage` itself -- the shipped code, called with its own
globals set -- not a re-implementation. The candidate's pools and fraction are built FOLD-WISE from
the other folds' labels, never from the labels being scored, because item 316's 0.0223 in-sample
gain was "pure selection contamination" and this file must not repeat it.

The score is the ST-RAE numerator the pair itself minimises, max(q-hi,0)+max(lo-q,0) averaged over
compounds, so candidate and baseline are commensurable by construction.

CONTROLS:
  * the affine baseline is the shipped function, and its own optimum is reported per enzyme;
  * an ORACLE arm places every compound at the centre of its own band, which bounds what any
    placement can reach and proves the scorer responds;
  * the candidate is also run with the pools taken IN SAMPLE, which must look better than the
    out-of-fold version -- if it does not, the fold machinery is not doing anything.

Part two writes the submission regardless of the verdict, so the decision stays the user's:
`activity_submission_mix2d6.csv` is `rankmap` with CYP2D6 alone re-placed, every other isoform bit
for bit identical.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import json

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

import submit
from submit import CYPS
from cypsplit import butina_folds

FLOOR = 4.0
SEEDS = [0, 1, 2, 3]
MEMBERS = ["поферментно", "пул", "GP", "гребневая", "ствол"]


def mixplace(p, lo_pool, hi_pool, q):
    """Monotone by construction: the top q fraction takes the ACTIVE pool's quantiles, the rest the
    inactive pool's, and every active value is >= FLOOR > every inactive one."""
    n = len(p)
    k = int(round(q * n))
    k = max(0, min(n, k))
    out = np.empty(n)
    order = np.argsort(p, kind="stable")
    if n - k > 0:
        out[order[: n - k]] = np.quantile(lo_pool, np.linspace(0, 1, n - k))
    if k > 0:
        out[order[n - k:]] = np.quantile(hi_pool, np.linspace(0, 1, k))
    return out


def pen(q, lo, hi):
    return float(np.mean(np.maximum(q - hi, 0.0) + np.maximum(lo - q, 0.0)))


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr["%s_pIC50_direct_inhibition" % c].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr["%s_pIC50_direct_inhibition_conf_low" % c].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr["%s_pIC50_direct_inhibition_conf_high" % c].to_numpy(float) for c in CYPS], 1)
    MASK = np.isfinite(Y)
    submit.LO, submit.HI = LO, HI          # the shipped fitter reads these as module globals
    fold, ncl = butina_folds(list(rows.SMILES))
    fold = np.asarray(fold)
    print("  сплит: %d кластеров; размеченных по ферментам %s"
          % (ncl, "/".join(str(int(MASK[:, e].sum())) for e in range(4))))

    tab = {c: {"affine": [], "mix_oof": [], "mix_insample": [], "oracle": []} for c in CYPS}
    for s in SEEDS:
        M = json.load(open(RES + "preds/members_seed%d.json" % s))
        P = [np.mean([np.asarray(M[m][e], float) for m in MEMBERS], 0) for e in range(4)]
        print("\n  --- сид %d, аффинная пара отгружаемым кодом ---" % s, flush=True)
        pair = submit.fit_shrinkage(P, Y, MASK)
        for e, c in enumerate(CYPS):
            m = MASK[:, e]
            idx = np.where(m)[0]
            y, lo, hi = Y[m, e], LO[m, e], HI[m, e]
            L, mu, sh = pair[e]
            tab[c]["affine"].append(pen(mu * (1 - L) + sh + L * P[e], lo, hi))
            tab[c]["oracle"].append(pen((lo + hi) / 2, lo, hi))
            # in-sample pools and fraction
            act, ina = y[y >= FLOOR], y[y < FLOOR]
            tab[c]["mix_insample"].append(
                pen(mixplace(P[e], ina, act, float((y >= FLOOR).mean())), lo, hi))
            # fold-wise: pools and fraction from the OTHER folds only
            q = np.empty(len(y))
            f = fold[idx]
            for k in np.unique(f):
                te, trn = f == k, f != k
                yt = y[trn]
                a, i2 = yt[yt >= FLOOR], yt[yt < FLOOR]
                if len(a) < 20 or len(i2) < 20:
                    q[te] = mu * (1 - L) + sh + L * P[e][te]      # not enough pool: leave as pair
                    continue
                q[te] = mixplace(P[e][te], i2, a, float((yt >= FLOOR).mean()))
            tab[c]["mix_oof"].append(pen(q, lo, hi))

    print("\n  ШТРАФ ST-RAE (числитель, меньше — лучше), среднее по %d сидам" % len(SEEDS))
    print("  %-8s %10s %12s %12s %10s %12s"
          % ("фермент", "аффинная", "смесь ВНЕ", "смесь в выб.", "оракул", "Δ вне-афф."))
    deltas = []
    for c in CYPS:
        a = np.mean(tab[c]["affine"]); o = np.mean(tab[c]["mix_oof"])
        i = np.mean(tab[c]["mix_insample"]); orc = np.mean(tab[c]["oracle"])
        sgn = sum(1 for k in range(len(SEEDS)) if tab[c]["mix_oof"][k] < tab[c]["affine"][k])
        deltas.append(o - a)
        print("  %-8s %10.4f %12.4f %12.4f %10.4f %+12.4f  знак %d/%d"
              % (c, a, o, i, orc, o - a, sgn, len(SEEDS)))
    print("  %-8s %10.4f %12.4f %12.4f %10.4f %+12.4f"
          % ("МАКРО", np.mean([np.mean(tab[c]["affine"]) for c in CYPS]),
             np.mean([np.mean(tab[c]["mix_oof"]) for c in CYPS]),
             np.mean([np.mean(tab[c]["mix_insample"]) for c in CYPS]),
             np.mean([np.mean(tab[c]["oracle"]) for c in CYPS]), np.mean(deltas)))

    orc_ok = all(np.mean(tab[c]["oracle"]) < np.mean(tab[c]["affine"]) for c in CYPS)
    ins_ok = all(np.mean(tab[c]["mix_insample"]) <= np.mean(tab[c]["mix_oof"]) + 1e-12 for c in CYPS)
    print("\n  КОНТРОЛЬ оракула (центр своей полосы должен бить всё): %s"
          % ("пройден" if orc_ok else "ПРОВАЛЕН — шкала не реагирует"))
    print("  КОНТРОЛЬ фолдов (в выборке должно быть не хуже, чем вне): %s"
          % ("пройден" if ins_ok else "ПРОВАЛЕН — фолды ничего не делают"))
    verdict = "ВЫИГРЫВАЕТ" if np.mean(deltas) < 0 else "ПРОИГРЫВАЕТ"
    print("  ВЕРДИКТ вне фолда: смесь %s аффинной паре на %+.4f макро"
          % (verdict, np.mean(deltas)))

    # ---------------- part two: build the submission regardless ----------------
    print("\n  --- сборка кандидата: rankmap с переставленным ТОЛЬКО CYP2D6 ---")
    k158 = json.load(open(RES + "logs/k158_derive.json"))
    qd = float(k158["derived"]["CYP2D6"]["p_active_mean"])
    src = pd.read_csv(RES + "submission/activity_submission_rankmap.csv")
    out = src.copy()
    col = "CYP2D6_pIC50_direct_inhibition"
    y = Y[MASK[:, 2], 2]
    act, ina = y[y >= FLOOR], y[y < FLOOR]
    v = mixplace(src[col].to_numpy(float), ina, act, qd)
    rho = float(spearmanr(src[col].to_numpy(float), v).statistic)
    if abs(rho - 1.0) > 1e-12:
        raise SystemExit("порядок CYP2D6 не сохранён, ро = %.9f" % rho)
    out[col] = v
    for c in CYPS:
        cc = "%s_pIC50_direct_inhibition" % c
        same = np.array_equal(out[cc].to_numpy(float), src[cc].to_numpy(float))
        if (c == "CYP2D6") == same:
            raise SystemExit("%s: ожидалось %s, получено обратное"
                             % (c, "изменение" if c == "CYP2D6" else "побитовое совпадение"))
    print("  КОНТРОЛЬ: ро CYP2D6 с источником %.6f; три прочих фермента побитово те же" % rho)
    print("  CYP2D6: было ср %.3f sd %.3f -> стало ср %.3f sd %.3f (доля активных %.3f)"
          % (src[col].mean(), src[col].std(ddof=1), v.mean(), v.std(ddof=1), qd))
    path = RES + "submission/activity_submission_mix2d6.csv"
    out.to_csv(path, index=False)
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    from validation.activity_validation import validate_activity_submission
    ok, errs = validate_activity_submission(path, expected_ids=set(te.Molecule_Name))
    bad, _ = validate_activity_submission(path, expected_ids={"нет-такой"})
    print("  валидатор: %s; КОНТРОЛЬ валидатора (неверные id): %s"
          % ("принято" if ok else "ОТКЛОНЕНО %s" % errs,
             "отказал как и ожидалось" if not bad else "ПРИНЯЛ — проверка не работает"))
    print("  записано %s" % path)
    json.dump({"table": tab, "macro_delta": float(np.mean(deltas)), "q_derived": qd},
              open(RES + "logs/k160_mixplace.json", "w"), ensure_ascii=False, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
