"""Score the union of the three endpoint-pretrained embeddings under the in-context learner, and
first prove that the LOCAL device reproduces the GPU box.

Item 376 registered the prediction before this ran: `U` -- `chemprop_medium` (600) concatenated with
`chemeleon` (2048) and `surrogate_adme` (2048), 4696 columns capped at PCA 100 -- should score above
both `E` (0.6409) and `S` (0.6387) alone, in 0.645 to 0.660, and is the first candidate with a route
to item 373's bar of 0.6480. If it gets there it replaces TWO slots with ONE, which the 1/n law says
is strictly better. The named risk was that PCA 100 over 4696 dimensions retains mostly the dominant
block's variance, leaving U at about max(E, S).

THE DEVICE CONTROL, WHICH MUST SUCCEED BEFORE U IS READ AT ALL. The GPU box became unreachable --
`vsb.fbb.msu.ru` stopped resolving, SERVFAIL from our resolver and nothing from three public ones,
while the control query (`huggingface.co`) resolved -- so these arms were produced locally on CPU
under torch 2.14.1 instead of the box's 2.13.0+cu130. A different torch build and a different
machine is exactly the kind of change that moves numbers quietly, so `E` and `S` were re-run here
and must reproduce the committed files from the box. If they do not, nothing in this file is
comparable with items 372, 373 or 374 and the gate says so instead of printing a table.

    uv run python verify/k154_union.py --loc <dir written by src/tabicl_run.py --device cpu>

Reads that directory's tabicl_oof.json plus the committed box files. Writes
results/logs/k154_union.json. Exit: 0 ran, 2 the device control failed.
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
from k106_probe import score, errcorr, members, probe_oof, load_emb

BOX = {"E": RES + "preds/k147_tabicl_E100.json", "S": RES + "preds/k147_tabicl_S100.json"}
FLOOR, BAR = 0.0036, 0.6480
PRED = (0.645, 0.660)


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--loc", required=True, help="каталог с tabicl_oof.json локального прогона")
    ap.add_argument("--tol", type=float, default=0.02,
                    help="допуск контроля устройства по max|Δ| предсказания")
    a = ap.parse_args(argv)
    L = json.load(open(a.loc.rstrip("/") + "/tabicl_oof.json"))

    rows = pd.read_csv(D + "rows.csv")
    z = np.load(D + "feats.npz")
    XF = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr["%s_pIC50_direct_inhibition" % c].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr["%s_pIC50_direct_inhibition_conf_low" % c].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr["%s_pIC50_direct_inhibition_conf_high" % c].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)
    K.S.LO, K.S.HI = LO, HI

    # ---- the device control -------------------------------------------------------------
    print("  КОНТРОЛЬ УСТРОЙСТВА: локальный CPU против карты, на тех же армах и сидах")
    worst = 0.0
    dev = {}
    for repk, path in BOX.items():
        B = json.load(open(path))
        for s in (0, 1, 2, 3):
            kb = "%s_seed%d" % (repk, s)
            if kb not in L["preds"]:
                raise SystemExit("в локальном прогоне нет %s" % kb)
            if str(B["digests"][s]) != str(L["digests"][s]):
                raise SystemExit("%s: дайджест %s локально против %s на карте"
                                 % (kb, L["digests"][s], B["digests"][s]))
            lb = [np.asarray(v, float) for v in B["preds"][kb]]
            ll = [np.asarray(v, float) for v in L["preds"][kb]]
            d = max(float(np.max(np.abs(x - y))) for x, y in zip(lb, ll))
            r = float(np.mean([spearmanr(x, y).statistic for x, y in zip(lb, ll)]))
            dev[kb] = {"max_abs": d, "spearman": r}
            worst = max(worst, d)
    byrep = {k: np.mean([v["max_abs"] for kk, v in dev.items() if kk.startswith(k)])
             for k in BOX}
    byrho = {k: np.mean([v["spearman"] for kk, v in dev.items() if kk.startswith(k)])
             for k in BOX}
    for k in BOX:
        print("    %s: средний max|Δ| %.4f, Spearman локального с картой %.6f"
              % (k, byrep[k], byrho[k]))
    print("    худший max|Δ| по всем восьми арм-сидам: %.4f (допуск %.4f)" % (worst, a.tol))
    if worst > a.tol:
        print("ОТКАЗ: локальное устройство не воспроизводит карту. Числа ниже не сравнимы с "
              "пунктами 372-374, и таблица не печатается.")
        return 2
    print("    ПРОЙДЕН\n")

    # ---- the measurement ----------------------------------------------------------------
    acc = {}
    for s in (0, 1, 2, 3):
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        assert str(L["digests"][s]) == dg, (L["digests"][s], dg)
        pdz = members(s, XF, Y, M, fold, LO, HI)
        probe = [probe_oof(load_emb(rows), Y, M, fold, e)[0] for e in range(len(CYPS))]
        chem = [np.asarray(v, float) for v in
                json.load(open(RES + "preds/k122_chemeleon_probe_seed%d.json" % s))["P"]]
        kept = [[Q[e] for k, Q in pdz if K._keep(e, k)] for e in range(len(CYPS))]
        arm = {k: [np.asarray(v, float) for v in L["preds"]["%s_seed%d" % (k, s)]]
               for k in ("E", "S", "U")}

        def put(name, P, base):
            sc = score(P, Y, LO, HI, M, fold)
            b = score(base, Y, LO, HI, M, fold)
            acc.setdefault(name, []).append({"rank": sc[0], "pair": sc[1],
                                             "d_rank": sc[0] - b[0], "d_pair": sc[1] - b[1]})

        six = [np.mean(kept[e] + [probe[e]], 0) for e in range(len(CYPS))]
        seven = [np.mean(kept[e] + [arm["E"][e], arm["S"][e]], 0) for e in range(len(CYPS))]
        put("шесть (отгружено)", six, six)
        put("семь: пять+E+S (п.374)", seven, six)
        put("шесть: пять+U", [np.mean(kept[e] + [arm["U"][e]], 0) for e in range(len(CYPS))], six)
        put("семь: пять+E+U", [np.mean(kept[e] + [arm["E"][e], arm["U"][e]], 0)
                               for e in range(len(CYPS))], six)
        put("восемь: пять+E+S+U", [np.mean(kept[e] + [arm["E"][e], arm["S"][e], arm["U"][e]], 0)
                                   for e in range(len(CYPS))], six)
        put("семь+chem: пять+chem+U", [np.mean(kept[e] + [chem[e], arm["U"][e]], 0)
                                       for e in range(len(CYPS))], six)
        for k in ("E", "S", "U"):
            put("одиночный %s" % k, arm[k], six)
            _p, vm, _v, _c = errcorr(arm[k], pdz, Y, M)
            acc.setdefault("_corr5_%s" % k, []).append(float(np.mean(vm)))
        for k in ("E", "S"):
            acc.setdefault("_corrU_%s" % k, []).append(
                float(np.mean([abs(np.corrcoef(Y[M[:, e], e] - arm["U"][e],
                                               Y[M[:, e], e] - arm[k][e])[0, 1])
                               for e in range(len(CYPS))])))

    print("  %-26s %9s %8s %10s %7s %10s" % ("арм", "ранг", "sd", "Δк шести", "знак", "ΔST-RAE"))
    out = {}
    for k, v in acc.items():
        if k.startswith("_"):
            continue
        r = [x["rank"] for x in v]
        g = [x["d_rank"] for x in v]
        out[k] = {"rank": float(np.mean(r)), "sd": float(np.std(r, ddof=1)),
                  "d_rank": float(np.mean(g)), "sign": int(sum(1 for x in g if x > 0)),
                  "d_pair": float(np.mean([x["d_pair"] for x in v]))}
        print("  %-26s %9.4f %8.4f %+10.4f %5d/%d %+10.4f"
              % (k, out[k]["rank"], out[k]["sd"], out[k]["d_rank"], out[k]["sign"], len(g),
                 out[k]["d_pair"]))
    corr = {k[1:]: float(np.mean(v)) for k, v in acc.items() if k.startswith("_")}
    print("\n  корреляция ошибок: U с пятью %.3f, E с пятью %.3f, S с пятью %.3f;"
          % (corr["corr5_U"], corr["corr5_E"], corr["corr5_S"]))
    print("                     U с E %.3f, U с S %.3f" % (corr["corrU_E"], corr["corrU_S"]))

    u = out["одиночный U"]["rank"]
    e_, s_ = out["одиночный E"]["rank"], out["одиночный S"]["rank"]
    print("\n  ЗАРЕГИСТРИРОВАННОЕ ПРЕДСКАЗАНИЕ (пункт 376, до прогона): U в %.3f--%.3f и выше "
          "обоих" % PRED)
    print("    U %.4f; E %.4f; S %.4f" % (u, e_, s_))
    ok_band = PRED[0] <= u <= PRED[1]
    ok_above = u > max(e_, s_)
    print("    в полосе: %s; выше обоих: %s; выше планки %.4f: %s"
          % ("ДА" if ok_band else "НЕТ", "ДА" if ok_above else "НЕТ", BAR,
             "ДА" if u > BAR else "НЕТ"))
    print("  пол ранга %.4f, пол ST-RAE 0.007" % FLOOR)
    json.dump({"device": dev, "device_by_rep": {"max_abs": byrep, "spearman": byrho},
               "arms": out, "corr": corr,
               "prediction": {"band": PRED, "bar": BAR, "in_band": ok_band,
                              "above_both": ok_above, "above_bar": bool(u > BAR)}},
              open(RES + "logs/k154_union.json", "w"), ensure_ascii=False, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
