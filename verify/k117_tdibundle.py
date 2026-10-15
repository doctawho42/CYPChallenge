"""The frozen-encoder head against the arm that ACTUALLY SHIPS, which items 337/338 did not do.

THE DEFECT THIS FILE REPAIRS. Items 337 and 338 measured the frozen `chemprop_medium` head against
`oof_gbm` on FP+DESC+MECH -- a plain direct classifier on `is_TDI`. That is the `--no-bundle` branch
of `src/submit.py`, the behaviour from BEFORE item 250. What ships is the conjunction: `--no-bundle`
is the opt-OUT, `results/submission/submission.meta.json` records the shipped arm as
"связка: ворота(ансамбль+мёртвая зона) * сдвиг", and `verify/k91_mccband.py`'s own docstring says it
in as many words -- `tdi_probs.json` holds "прямой классификатор, то есть плечо, которое больше не
подаётся". Every control inside items 337/338 held bit for bit, and none of them could catch this,
because they faithfully compared arm C with arm A. A control that the baseline reproduces is only
worth something when the baseline is the right one.

The size of the error, indicatively: at seed 0 the shipped bundle scores macro MCC 0.2430 against arm
A's 0.2246, so roughly a third of the +0.0557 that item 338 reported for arm C was a gain over a
weaker baseline. Indicatively only -- the bundle's rows are not arm A's rows (2334 against 2346 on
CYP3A4) and cross-script comparison is what this file exists to stop doing.

WHAT IS MEASURED HERE. The shipped conjunction, reproduced by the same functions item 250 used ---
`member_oof` and `dz_refit` from `k75_armens`, `platt_oof`, `oof_clf` and `decide` from `k76_bundle`,
assembled as `k91_mccband` assembles it --- and then three variants of adding the encoder, all on ONE
row mask and ONE fold vector:

    S     подаваемая связка            gp * dp, затем decide
    S+E   связка плюс вложение         среднее калиброванных, затем decide
    E     вложение одно                через то же устройство
    S+D   РАЗЛИЧАЮЩИЙ КОНТРОЛЬ         то же, но голова на НАШЕМ DESC+MECH

`S+D` is what decides what a win by `S+E` means. If the linear head on our own block also lifts the
bundle, the gain is the head and not the frozen weights --- the test item 324 built and items 337/338
passed on the wrong baseline.

THE CONTROL THAT MAKES THE BASELINE TRUSTWORTHY THIS TIME. `results/preds/bundle_oof.json` is the
committed per-compound output of the shipped arm at seed 0, written by `k91_mccband`. The bundle
reproduced here must match it --- same n, and MCC equal to 0.3578 on CYP3A4 and 0.1282 on CYP2D6. A
mismatch means this file is measuring something else again, and it exits rather than reporting.

    uv run python verify/k117_tdibundle.py --seeds 0

The four ensemble members with their dead-zone refits are the cost, about an hour a seed. Everything
is cached under the fold digest, so adding a variant later costs seconds.

Writes results/logs/k117_tdibundle.json. Exit: 0 ran, 2 a control failed, 3 inputs missing.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json
import time

import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef, roc_auc_score

from cypsplit import butina_folds, fold_digest
from k75_armens import member_oof, dz_refit, TDI, GATE, L
from k76_bundle import platt_oof, oof_clf, decide      # same functions as item 250
from k106_probe import load_emb
# The head comes from the SHIPPED side, so the arm measured here and the arm
# src/submit.py builds under --tdiemb are one piece of code, not two copies.
from submit import tdiemb_oof as oof_logit

KINDS = ["поферментно", "пул", "GP", "гребневая"]      # k91's four, the shipped arm
LOG = RES + "logs/k117_tdibundle.json"
CACHE = RES + "preds/k117_bundle_seed%d.json"
# Item 235's floors, measured on the shipped path over four seeds.
FLOOR = {"CYP3A4": 0.0281, "CYP2D6": 0.0419, "макро": 0.0076}
# results/preds/bundle_oof.json, committed, seed 0, the shipped arm per compound.
EXPECT = {"CYP3A4": (2334, 0.3578), "CYP2D6": (1493, 0.1282)}


def build(seed, digest, Xs, ys, folds, LOs, HIs, shifts, gates):
    """`gp` and `dp` of the shipped conjunction, cached under the fold digest.

    This is the hour. Caching it means a later variant of the COMBINATION -- which is a question
    about how the encoder enters, not about the bundle -- costs seconds. Items 335 and 338 each
    paid the whole bill twice for exactly that distinction.
    """
    path = CACHE % seed
    store = json.load(open(path)) if _pl.Path(path).exists() else {}
    if store.get("digest") != digest:
        store = {"digest": digest}
    if "gp" in store and "dp" in store:
        print("      связка из кэша", flush=True)
        return ([np.asarray(v, float) for v in store["gp"]],
                [np.asarray(v, float) for v in store["dp"]])
    dz = {}
    for k in KINDS:
        t0 = time.time()
        P = member_oof(k, Xs, ys, folds)
        tg = [np.clip(P[e], LOs[e], HIs[e]) for e in range(len(TDI))]
        dz[k] = dz_refit(k, Xs, tg, folds)
        print(f"      плечо, член {k}: {time.time() - t0:.0f} с", flush=True)
    arm = [np.mean([dz[k][e] for k in KINDS], 0) for e in range(len(TDI))]
    gp, dp = [], []
    for e, c in enumerate(TDI):
        dp.append(oof_clf(Xs[e], shifts[e] > L, folds[e], f"{c} сдвиг"))
        gp.append(platt_oof(arm[e] - GATE, gates[e], folds[e], raw=True))
    store["gp"] = [[float(x) for x in v] for v in gp]
    store["dp"] = [[float(x) for x in v] for v in dp]
    _pl.Path(RES + "preds").mkdir(parents=True, exist_ok=True)
    json.dump(store, open(path, "w"))
    return gp, dp


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    for f in ("emb_chemprop_medium.npz", "feats.npz", "rows.csv"):
        if not _pl.Path(D + f).exists():
            print(f"нет data/{f}")
            return 3

    rows = pd.read_csv(D + "rows.csv")
    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    XOWN = np.hstack([z["DESC"], z["MECH"]])
    E = load_emb(rows)
    tdi = (pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv")
             .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    inh = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
             .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())

    out = json.load(open(LOG)) if _pl.Path(LOG).exists() else {}
    out.setdefault("meta", {}).update(floors=FLOOR, kinds=KINDS, expect=EXPECT)
    out.setdefault("seeds", {})

    for s in seeds:
        fold, ncl = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        print(f"\n  сид {s}: дайджест {dg}, кластеров {ncl}", flush=True)
        if s == 0 and dg != "2d93c19815e14261":
            print(f"ОТКАЗ: разбиение сдвинулось, сид 0 дал {dg}")
            return 2

        # The bundle's own row mask: it needs the pre-incubation arm and its band, which is
        # STRICTER than `lab.notna()` and is why items 337/338's row set was not this one.
        Xs, ys, folds, LOs, HIs, labs, shifts, gates, Es, Os = ([] for _ in range(10))
        for c in TDI:
            lab = tdi[f"{c}_is_TDI"]
            ta = tdi[f"{c}_pIC50_TDI_condition"].to_numpy(float)
            dr = inh[f"{c}_pIC50_direct_inhibition"].to_numpy(float)
            lo = tdi[f"{c}_pIC50_TDI_condition_conf_low"].to_numpy(float)
            hi = tdi[f"{c}_pIC50_TDI_condition_conf_high"].to_numpy(float)
            m = (lab.notna() & np.isfinite(ta) & np.isfinite(dr)
                 & np.isfinite(lo) & np.isfinite(hi)).to_numpy()
            Xs.append(X[m]); ys.append(ta[m]); folds.append(fold[m])
            LOs.append(lo[m]); HIs.append(hi[m])
            labs.append(lab[m].astype(int).to_numpy())
            shifts.append((ta - dr)[m]); gates.append((ta > GATE)[m].astype(int))
            Es.append(E[m]); Os.append(XOWN[m])
            print(f"    {c}: n {int(m.sum())}, доля метки "
                  f"{lab[m].astype(bool).mean():.3f}", flush=True)
            if int(m.sum()) != EXPECT[c][0]:
                print(f"ОТКАЗ: {c} дал {int(m.sum())} строк против {EXPECT[c][0]} "
                      f"в bundle_oof.json --- это не та маска")
                return 2

        gp, dp = build(s, dg, Xs, ys, folds, LOs, HIs, shifts, gates)

        res = {}
        for e, c in enumerate(TDI):
            y, fi = labs[e], folds[e]
            dec_s, cal_s = decide(gp[e] * dp[e], y, fi)
            mcc_s = float(matthews_corrcoef(y, dec_s))

            # CONTROL: the bundle reproduced here must BE the committed shipped arm.
            if s == 0 and abs(mcc_s - EXPECT[c][1]) > 5e-4:
                print(f"ОТКАЗ: {c} связка дала MCC {mcc_s:.4f} против {EXPECT[c][1]:.4f} "
                      f"в bundle_oof.json --- измеряется не отгруженное плечо")
                return 2

            pe = oof_logit(Es[e], y, fi)             # frozen 600 dimensions
            po = oof_logit(Os[e], y, fi)             # our DESC+MECH, the control
            cal_e = platt_oof(pe, y, fi)
            cal_o = platt_oof(po, y, fi)
            arms = {"S связка": (gp[e] * dp[e], dec_s, mcc_s),
                    "S+E связка+вложение": ((cal_s + cal_e) / 2, None, None),
                    "E вложение одно": (pe, None, None),
                    "S+D КОНТРОЛЬ DESC+MECH": ((cal_s + cal_o) / 2, None, None)}
            cell = {}
            for nm, (score, dcached, mcached) in arms.items():
                if dcached is None:
                    d, _ = decide(score, y, fi)
                    mv = float(matthews_corrcoef(y, d))
                else:
                    d, mv = dcached, mcached
                cell[nm] = {"mcc": mv, "auc": float(roc_auc_score(y, score)),
                            "rate": float(d.mean())}
                print("      %-26s MCC %+.4f  Δ к S %+.4f  AUC %.4f  доля реш. %.3f"
                      % (nm, mv, mv - mcc_s, cell[nm]["auc"], cell[nm]["rate"]), flush=True)
            cell["истинная_доля"] = float(y.mean())
            cell["n"] = int(len(y))
            res[c] = cell

        out["seeds"][str(s)] = res
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
        print(f"      записано: {LOG}", flush=True)

    # ---- the table ---------------------------------------------------------------------
    ARMS = ["S связка", "S+E связка+вложение", "E вложение одно", "S+D КОНТРОЛЬ DESC+MECH"]
    print("\n" + "=" * 94)
    print("ПРОТИВ ОТГРУЖЕННОЙ СВЯЗКИ (не против прямого классификатора, как в пунктах 337/338)")
    print("=" * 94)
    print("  %-26s %9s %9s   %9s %9s   %9s" % ("арм", "CYP3A4", "Δ к S", "CYP2D6", "Δ к S",
                                               "макро Δ"))
    got = {a: {} for a in ARMS}
    for c in TDI:
        for a in ARMS:
            got[a][c] = float(np.mean([out["seeds"][str(s)][c][a]["mcc"] for s in seeds]))
    for a in ARMS:
        d3 = got[a]["CYP3A4"] - got["S связка"]["CYP3A4"]
        d2 = got[a]["CYP2D6"] - got["S связка"]["CYP2D6"]
        print("  %-26s %9.4f %+9.4f   %9.4f %+9.4f   %+9.4f"
              % (a, got[a]["CYP3A4"], d3, got[a]["CYP2D6"], d2, (d3 + d2) / 2))
    print("  %-26s %9.4f %9s   %9.4f %9s   %9.4f"
          % ("пол (пункт 235)", FLOOR["CYP3A4"], "", FLOOR["CYP2D6"], "", FLOOR["макро"]))
    print(f"\n  сиды: {seeds}. Решает S+E против S при провале S+D.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
