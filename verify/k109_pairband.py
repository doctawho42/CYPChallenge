"""The PAIRED reveal band on macro rank --- the instrument item 294 named and did not build.

Item 294 measured the reveal band on the ABSOLUTE macro rank: half-width 0.0270 at n=375, the live
board's size. Item 327 then read the probe's board move of +0.0127 against it and had to say the
board confirms the gain without establishing it, because the move is less than half that band.

That comparison is the wrong one, and `k92_rankband.py` says so in its own closing note: the board
scores every submission on the SAME molecules, so per-compound noise cancels in a DIFFERENCE, and
the project's gains must be compared against a paired floor rather than against the spread of one
absolute score. Item 327 recorded the paired band as unmeasured without crediting k92 for having
named it. This file measures it.

WHAT IS AND IS NOT COMPUTABLE, said first because the obvious reading of "compute the paired band"
is not available. A paired band on the BOARD's own difference cannot be computed at all: it would
need the identity of the 375 revealed compounds AND their labels, and the labels are the blind half.
What is computable is the paired band under the SAME procedure item 294 used -- subsample our own
out-of-fold predictions without replacement down to 375 -- which answers the question actually at
issue: how much of the +-0.0270 is shared-sample noise that a difference cancels.

THE CONTROL WITHOUT WHICH THIS FILE IS DECORATION. The procedure must be item 294's, not a
lookalike. So the run first reproduces k92's own logged macro band on `oof_submitted.json` from
`results/logs/k92_rankband.json`, to the fourth decimal, using the same rng seeding
(`BOOT_SEED + n`). If that misses, the numbers below are measuring something else and the run says
so instead of printing them.

PAIRING IS THE WHOLE POINT, so the draws are shared: one set of subsample indices per enzyme per
resample, evaluated under BOTH arms. Drawing separately would reproduce the unpaired band twice and
call their difference paired, which is the defect this file exists to avoid.

    uv run python verify/k109_pairband.py                 # seed 0, the analogue of k92
    uv run python verify/k109_pairband.py --seeds 0,1,2,3,4,5,6,7

Writes results/logs/k109_pairband.json. Runtime: about a minute per seed.
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

from cypsplit import butina_folds
import submit as S
from submit import CYPS, _combine
from k92_rankband import BOOT_SEED, N_BOOT, rho, summarise          # one definition, not a copy
from k106_probe import combine_add, load_emb, members, probe_oof

SIZES = [750, 375]
LOG = RES + "logs/k109_pairband.json"

# item 294's own logged figures, to reproduce rather than to trust
K92_LOG = RES + "logs/k92_rankband.json"


def draws(n, N_by_e):
    """One shared set of subsample indices per enzyme per resample.

    Seeded exactly as k92 (`BOOT_SEED + n`, one stream per enzyme in CYPS order), so the
    subsamples here ARE k92's subsamples and the control below can be exact rather than close.
    """
    out = []
    for e in range(len(CYPS)):
        rng = np.random.default_rng(BOOT_SEED + n)
        # k92 makes a fresh generator per enzyme with this same seed, so each enzyme's stream
        # restarts; reproduced rather than tidied, because tidying it would move the numbers.
        out.append([rng.choice(N_by_e[e], n, replace=False) for _ in range(N_BOOT)])
    return out


def band_of(vecs, Y, mask, idx, n):
    """Macro rank across one set of shared draws. `vecs` is per-enzyme predictions."""
    per_e = []
    for e in range(len(CYPS)):
        y, p = Y[mask[:, e], e], np.asarray(vecs[e], float)
        per_e.append(np.array([rho(y[i], p[i]) for i in idx[e]]))
    return np.mean(per_e, axis=0)


def control(Y, mask):
    """Reproduce item 294's macro band on the arm it was measured on, or refuse."""
    want = json.load(open(K92_LOG))["макро полоса n=375"]
    P = [np.asarray(v, float) for v in json.load(open(RES + "preds/oof_submitted.json"))["P"]]
    N_by_e = [int(mask[:, e].sum()) for e in range(len(CYPS))]
    got = summarise(band_of(P, Y, mask, draws(375, N_by_e), 375))
    hw_w = (want["hi95"] - want["lo95"]) / 2
    hw_g = (got["hi95"] - got["lo95"]) / 2
    ok = all(abs(got[k] - want[k]) < 5e-5 for k in ("lo95", "hi95", "sd"))
    print("КОНТРОЛЬ: воспроизводится ли полоса пункта 294 этой же процедурой\n")
    print(f"    пункт 294 (лог k92)   [{want['lo95']:+.4f}, {want['hi95']:+.4f}]  "
          f"sd {want['sd']:.4f}  полуширина {hw_w:.4f}")
    print(f"    здесь                 [{got['lo95']:+.4f}, {got['hi95']:+.4f}]  "
          f"sd {got['sd']:.4f}  полуширина {hw_g:.4f}")
    print(f"    -> {'СОВПАЛО, процедура та же' if ok else 'НЕ СОВПАЛО'}\n")
    return ok, hw_w


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    rows = pd.read_csv(D + "rows.csv")
    z = np.load(D + "feats.npz")
    XFULL = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(Y)
    S.LO, S.HI = LO, HI
    N_by_e = [int(mask[:, e].sum()) for e in range(len(CYPS))]

    ok, hw_unpaired_294 = control(Y, mask)
    if not ok:
        raise SystemExit("процедура расходится с пунктом 294 --- числа ниже мерили бы другое")

    E = load_emb(rows)
    out = {"meta": {"n_boot": N_BOOT, "seed_rng": BOOT_SEED, "sizes": SIZES,
                    "полуширина_непарной_294": hw_unpaired_294}, "seeds": {}}

    for s in seeds:
        print(f"сид {s}", flush=True)
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        pdz = members(s, XFULL, Y, mask, fold, LO, HI)
        probe = [probe_oof(E, Y, mask, fold, e)[0] for e in range(len(CYPS))]
        A = [_combine(pdz, e) for e in range(len(CYPS))]
        Dm = combine_add(pdz, probe)
        rec = {}
        for n in SIZES:
            idx = draws(n, N_by_e)
            a_ = band_of(A, Y, mask, idx, n)
            d_ = band_of(Dm, Y, mask, idx, n)
            dif = d_ - a_
            sa, sd_, sdif = summarise(a_), summarise(d_), summarise(dif)
            hw = lambda t: (t["hi95"] - t["lo95"]) / 2
            rec[str(n)] = {"A": sa, "D": sd_, "D-A": sdif,
                           "полуширина_A": hw(sa), "полуширина_D": hw(sd_),
                           "полуширина_разности": hw(sdif),
                           # Stored, not only printed: the share of draws on which the probe
                           # wins is the sharpest single number here, and a log that cannot be
                           # re-read for it makes the claim rest on a terminal scrollback.
                           "доля_положительных": float((dif > 0).mean())}
            print(f"    n={n}")
            print(f"        арм A  среднее {sa['mean']:+.4f}  полуширина {hw(sa):.4f}")
            print(f"        арм D  среднее {sd_['mean']:+.4f}  полуширина {hw(sd_):.4f}")
            print(f"        D - A  среднее {sdif['mean']:+.4f}  sd {sdif['sd']:.4f}  "
                  f"полуширина {hw(sdif):.4f}  "
                  f"(уже непарной в {hw(sa)/hw(sdif):.1f} раза)")
            print(f"        доля розыгрышей, где разность > 0: "
                  f"{float((dif > 0).mean()):.3f}")
        out["seeds"][str(s)] = rec

    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
    print(f"\nсохранено: {LOG}")
    print("""
Как читать. Полоса РАЗНОСТИ --- это то, с чем надо сравнивать наши приросты по рангу, и именно её
отсутствие заставило пункт 327 сказать «доска подтверждает, но не устанавливает». Полоса
АБСОЛЮТНОГО счёта (пункт 294) осталась той же и по-прежнему нужна для другого: сказать, лежит ли
наш абсолютный ранг там, где мы его ждали.

Чего здесь НЕТ, и это не мелочь. Полоса посчитана подвыборкой НАШИХ вневыборочных предсказаний,
то есть на нашем распределении. Тест лежит там, где у кросс-валидации почти нет массы (пункт 123,
потолок ESS 26.1%; пункт 129, медиана ближайшего соседа 0.587 против <=0.450). Полоса на ПАРНОЙ
разности самой доски не вычислима вовсе: для неё нужны и список раскрытых 375 соединений, и их
метки, а метки --- это и есть слепая половина.""")


if __name__ == "__main__":
    raise SystemExit(main())
