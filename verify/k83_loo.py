"""Реестр вкладов, перемеренный ВЫБИВАНИЕМ из подаваемой конфигурации. Пункт 272.

Четыре слагаемых табло мерялись как ДОБАВЛЕНИЕ к той базе, какая была на момент предложения.
Пункт 269 показал, почему сумме таких чисел верить нельзя: +0.0108 над одним членом, +0.0067
над более сильным и -0.0007 над ансамблем --- не от разбавления, а от избыточности. Значит
размер каждой строки реестра зависит от базы и от порядка добавления, а выбить компонент из
готовой модели не пробовал никто.

Эталон --- то, что подаётся: oof_members(mode="ансамбль5", dead=True), dz_pass и _keep, то есть
SOLO = {"CYP3A4": ("GP",)} действует. Выбивания K1, K2, K3, K5, K6 суть рекомбинации ОДНОЙ
сборки членов и стоят сверх неё ноль; K4 требует второй сборки на FP+DESC.

K4 --- ЧАСТИЧНОЕ выбивание и помечено так: ствол это torch-модель на DESC+MECH, закоммиченная
предсказаниями, свой механистический блок он сохраняет. K4 занижает цену MECH.

Условия и прогноз записаны в пункте 272 до прогона. Пишет results/logs/k83_loo.json после
каждого сида.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse, json, time
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds
from shrinkchoice import fit_apply
import submit as S
from submit import CYPS, oof_members, dz_pass, _keep


def combine(parts, mask, drop=None):
    """Среднее членов по правилу подачи, с выбитым членом. Если после выбивания на ферменте
    не остаётся ни одного члена (SOLO на CYP3A4 при drop='GP'), берутся все прочие."""
    out, fell = [], []
    for e in range(len(CYPS)):
        keep = [P[e] for k, P in parts if k != drop and _keep(e, k)]
        if not keep:
            keep = [P[e] for k, P in parts if k != drop]
            fell.append(CYPS[e])
        out.append(np.mean(keep, axis=0))
    return out, fell


def score(P, y, LO, HI, mask, fold):
    rk, pr = [], []
    for e in range(len(CYPS)):
        m = mask[:, e]
        yy, lo, hi = y[m, e], LO[m, e], HI[m, e]
        u = np.ones(len(yy)) / len(yy)
        rk.append(float(spearmanr(yy, P[e]).statistic))
        pr.append(float(strae(yy, fit_apply(P[e], lo, hi, fold[m], u),
                              y_true_upper=hi, y_true_lower=lo)))
    return float(np.mean(rk)), float(np.mean(pr)), rk, pr


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    ap.add_argument("--no-mech", dest="mech", action="store_false",
                    help="пропустить K4 (вторая сборка на FP+DESC), вдвое дешевле")
    a = ap.parse_args()

    z = np.load(D + "feats.npz")
    XFULL = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    XNOMECH = np.hstack([z["FP"], z["DESC"]])
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)
    S.LO, S.HI = LO, HI                       # submit держит их в модульных глобалах

    dst = RES + "logs/k83_loo.json"
    out = json.load(open(dst)) if _pl.Path(dst).exists() else {}
    print(f"  {'сид':>3s} {'конфигурация':>26s} {'макро-ранг':>11s} {'макро-пара':>11s} "
          f"{'Δранг':>9s}", flush=True)
    for s in [int(x) for x in a.seeds.split(",") if x.strip()]:
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        res = out.get(str(s), {})
        t0 = time.time()
        parts = oof_members(XFULL, y, mask, fold, "ансамбль5", seed=s, dead=True)
        print(f"      сборка членов: {time.time() - t0:.0f} c", flush=True)
        t0 = time.time()
        pdz = dz_pass(parts, XFULL, mask, fold, (LO, HI))
        print(f"      проход мёртвой зоны: {time.time() - t0:.0f} c", flush=True)

        ref, _ = combine(pdz, mask)
        res["эталон"] = score(ref, y, LO, HI, mask, fold)
        base = res["эталон"][0]
        print(f"  {s:3d} {'эталон (подаётся)':>26s} {res['эталон'][0]:11.4f} "
              f"{res['эталон'][1]:11.4f} {'':>9s}", flush=True)

        cfgs = [("K1 без мёртвой зоны", parts, None),
                ("K2 без пула", pdz, "пул"),
                ("K3 без ствола", pdz, "ствол"),
                ("K5 без GP", pdz, "GP"),
                ("K6 без гребневой", pdz, "гребневая")]
        for nm, pp, drop in cfgs:
            P, fell = combine(pp, mask, drop)
            sc = score(P, y, LO, HI, mask, fold)
            res[nm] = sc
            tag = f"  (SOLO пуст на {', '.join(fell)}, взяты все прочие)" if fell else ""
            print(f"  {s:3d} {nm:>26s} {sc[0]:11.4f} {sc[1]:11.4f} {sc[0]-base:+9.4f}{tag}",
                  flush=True)

        if a.mech:
            t0 = time.time()
            p2 = oof_members(XNOMECH, y, mask, fold, "ансамбль5", seed=s, dead=True)
            p2 = dz_pass(p2, XNOMECH, mask, fold, (LO, HI))
            P, _ = combine(p2, mask)
            sc = score(P, y, LO, HI, mask, fold)
            res["K4 без MECH (частичное)"] = sc
            print(f"  {s:3d} {'K4 без MECH (частичное)':>26s} {sc[0]:11.4f} {sc[1]:11.4f} "
                  f"{sc[0]-base:+9.4f}   ({time.time() - t0:.0f} c)", flush=True)

        out[str(s)] = res
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(dst, "w"), ensure_ascii=False, indent=1)

    names = [k for k in out[list(out)[0]] if k != "эталон"]
    print(f"\n  {'выбивание':>26s} {'Δранг ср.':>10s} {'sd':>8s} {'знак<0':>8s} {'Δпара':>9s}"
          f"   реестр", flush=True)
    LEDGER = {"K1 без мёртвой зоны": "+0.0197", "K4 без MECH (частичное)": "+0.0163",
              "K2 без пула": "+0.0141", "K3 без ствола": "+0.0054"}
    tot = 0.0
    for nm in names:
        ks = [k for k in out if nm in out[k] and "эталон" in out[k]]
        if not ks:
            continue
        d = np.array([out[k][nm][0] - out[k]["эталон"][0] for k in ks])
        q = np.array([out[k][nm][1] - out[k]["эталон"][1] for k in ks])
        if nm in LEDGER:
            tot += -d.mean()
        print(f"  {nm:>26s} {d.mean():+10.4f} {d.std(ddof=1) if len(d) > 1 else 0:8.4f} "
              f"{int((d < 0).sum()):5d}/{len(d):<2d} {q.mean():+9.4f}   {LEDGER.get(nm, '---')}",
              flush=True)
    print(f"\n  сумма четырёх выбиваний реестра: {tot:+.4f}   против суммы добавлений 0.0555 "
          f"и траектории 0.0579", flush=True)


if __name__ == "__main__":
    main()
