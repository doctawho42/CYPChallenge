"""Item 316's closure, put to the test it names in its own caveat: a RidgeCV probe on FROZEN
chemprop_medium embeddings, entered as a SIXTH member and scored BY RANK over the shipped
ensemble after the affine pair, nested, four seeds.

WHY THIS IS NOT ALREADY CLOSED. Items 61/117/154 closed frozen embeddings as CONCATENATED
COLUMNS into the same learner, and item 61 deliberately used the organisers' rdkit2d
checkpoint -- pretrained on the very descriptors our DESC block holds, i.e. the encoder's
worst case as a decorrelated member. This is a SEPARATE head whose PREDICTION enters the
mean. Item 316 measured a sixth member at any fitted weight as exactly +0.0000 and then
wrote its own limit: "this is about OUR library --- five learners over one feature matrix ---
not about a library of architecturally distinct models". This is that member.

ARMS, scored over the ENSEMBLE, never over the member:

  A  REFERENCE : the shipped composition (oof_members ансамбль5 + dz_pass + _keep/SOLO),
                 from the committed members_seed{s}.json cache.
  D  ARM       : A with the probe added as a sixth member at EQUAL weight. The probe gets no
                 dead-zone pass, exactly as ствол gets none -- it is the other member that is
                 not refitted here.
  Dz SENSITIVITY: as D, with the probe dead-zone-passed on its own design (same RidgeCV, same
                 alphas, target = clip(p_oof, lo, hi)). Reported because the dz pass is worth
                 +0.0167 of rank on the shipped five (item 164) and leaving it off could
                 understate the probe. NOT the deciding arm.
  W  ITEM 316's ARM: the probe at a NESTED FITTED weight per (enzyme, outer fold), chosen on
                 the other four folds by Spearman and applied held out.

THE ALPHA GRID IS CARRIED IN, NOT SEARCHED. JacksonBurns's committed notes record that small
alphas overfit badly ("-0.04 garbage"); the grid is [100, 1000, 10000] and sklearn's default
0.1/1/10 is not used. RidgeCV's own generalised-CV picks within that grid on TRAINING rows
only, per fold -- so it is nested, and no alpha is chosen against the outcome.

THE GATE, pre-registered and not widened: >= +0.0031 of macro Spearman over the ensemble with
the sign holding on at least 3 of 4 seeds. One board place costs 0.0050 of macro ST-RAE
(read directly off results/leaderboard_2026-09-25T0644Z.json: we are 0.5535, the entry above
is 0.5485); converting that to rank needs the slope, which item 322 remeasured at -1.62 and
item 323 qualified. At -1.62 a place is +0.0031 of rank -- BELOW the 0.0036 macro floor -- so
anything in [0.0031, 0.0036] is INSIDE THE FLOOR and is not a pass. Item 323 further shows
the windows nearest our own rank fit shallower (-0.97 at ranks 43-73, -0.58 at 48-68), so
+0.0031 is a LOWER BOUND on the price of a place, not the price; the report prints the
conversion at all three slopes rather than quoting one.
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
from sklearn.linear_model import RidgeCV
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae

from cypsplit import butina_folds, fold_digest
from shrinkchoice import fit_apply
import submit as S
from submit import CYPS, oof_members, dz_pass, _keep, _oof_trunk, SOLO

EMB = D + "emb_chemprop_medium.npz"
LOG = RES + "logs/k106_probe.json"
ALPHAS = np.array([100.0, 1000.0, 10000.0])
WGRID = np.round(np.arange(0.0, 1.0001, 0.05), 4)
GATE = 0.0031          # one board place at the causal slope -1.62 (items 322, 323)
FLOOR = 0.0036         # project macro-rank noise floor
PLACE_STRAE = 0.0050   # measured off the committed board, not quoted from memory
MEMBERS = ("поферментно", "пул", "GP", "гребневая", "ствол")


def score(P, y, LO, HI, mask, fold):
    """verify/k103_thirdhead.py's score(), unchanged: rank = mean Spearman
    (affine-pair-invariant), pair = mean ST-RAE AFTER the per-fold affine pair."""
    rk, pr = [], []
    for e in range(len(CYPS)):
        m = mask[:, e]
        yy, lo, hi = y[m, e], LO[m, e], HI[m, e]
        u = np.ones(len(yy)) / len(yy)
        rk.append(float(spearmanr(yy, P[e]).statistic))
        pr.append(float(strae(yy, fit_apply(P[e], lo, hi, fold[m], u),
                              y_true_upper=hi, y_true_lower=lo)))
    return float(np.mean(rk)), float(np.mean(pr)), rk, pr


def combine_keep(parts):
    return [np.mean([P[e] for k, P in parts if _keep(e, k)], axis=0) for e in range(len(CYPS))]


def combine_add(parts, probe):
    """The shipped mean with the probe added as one more equally weighted member."""
    return [np.mean([P[e] for k, P in parts if _keep(e, k)] + [probe[e]], axis=0)
            for e in range(len(CYPS))]


def combine_fitted(parts, probe, y, mask, fold):
    """Item 316's arm: a per-(enzyme, outer fold) weight on the probe, chosen NESTED.

    The weight is chosen on the four training folds by Spearman against the labels there and
    applied to the held-out fold. The member predictions are already out of fold, so no label
    of the held-out fold touches the choice.
    """
    out, chosen = [], []
    for e in range(len(CYPS)):
        m = mask[:, e]
        yy, fm = y[m, e], fold[m]
        base = np.mean([P[e] for k, P in parts if _keep(e, k)], axis=0)
        pr = probe[e]
        blend = np.zeros(len(yy))
        for f in range(5):
            te, trn = fm == f, fm != f
            if te.sum() == 0:
                continue
            sc = [float(spearmanr(yy[trn], (1 - a) * base[trn] + a * pr[trn]).statistic)
                  for a in WGRID]
            a = float(WGRID[int(np.argmax(sc))])
            blend[te] = (1 - a) * base[te] + a * pr[te]
            chosen.append({"enzyme": CYPS[e], "fold": int(f), "a": a,
                           "rho_at_a": sc[int(np.argmax(sc))], "rho_at_0": sc[0]})
        out.append(blend)
    return out, chosen


def probe_oof(E, y, mask, fold, e):
    """RidgeCV on the frozen embedding, out of fold on OUR Butina folds, one enzyme.

    Standardisation and the dead-column filter are fitted on the TRAINING rows of each fold
    and applied to the held-out rows -- never on the union, which is item 202's defect 2.
    """
    m = mask[:, e]
    Xm, ym, fm = E[m], y[m, e], fold[m]
    p = np.full(int(m.sum()), np.nan)
    picked, ncols = [], []
    for f in range(5):
        trn, te = fm != f, fm == f
        if te.sum() == 0:
            continue
        mu, sd = Xm[trn].mean(0), Xm[trn].std(0)
        keep = sd > 1e-6
        Z = (Xm[:, keep] - mu[keep]) / sd[keep]
        r = RidgeCV(alphas=ALPHAS).fit(Z[trn], ym[trn])
        p[te] = r.predict(Z[te])
        picked.append(float(r.alpha_))
        ncols.append(int(keep.sum()))
    if not np.isfinite(p).all():
        raise SystemExit(f"{CYPS[e]}: зонд оставил не заполненные строки")
    return p, picked, ncols


def probe_dz(E, y, mask, fold, e, p_oof, LO, HI):
    """The dead-zone pass for the probe: same design, target = clip(p_oof, lo, hi)."""
    m = mask[:, e]
    Xm, fm = E[m], fold[m]
    tgt = np.clip(p_oof, LO[m, e], HI[m, e])
    q = np.full(len(tgt), np.nan)
    for f in range(5):
        trn, te = fm != f, fm == f
        if te.sum() == 0:
            continue
        mu, sd = Xm[trn].mean(0), Xm[trn].std(0)
        keep = sd > 1e-6
        Z = (Xm[:, keep] - mu[keep]) / sd[keep]
        q[te] = RidgeCV(alphas=ALPHAS).fit(Z[trn], tgt[trn]).predict(Z[te])
    return q


def members(s, XFULL, y, mask, fold, LO, HI):
    """k103_thirdhead.py's members(), unchanged, including its bit-for-bit trunk control."""
    p = RES + f"preds/members_seed{s}.json"
    if not _pl.Path(p).exists():
        parts = oof_members(XFULL, y, mask, fold, "ансамбль5", seed=s, dead=True)
        return dz_pass(parts, XFULL, mask, fold, (LO, HI))
    mem = json.load(open(p))
    if set(mem) != set(MEMBERS):
        raise SystemExit(f"кэш членов сида {s} содержит {sorted(mem)}, а не {sorted(MEMBERS)}")
    pdz = [(k, [np.asarray(mem[k][e], float) for e in range(len(CYPS))]) for k in MEMBERS]
    fresh = _oof_trunk(y, mask, fold, s, dead=True)
    cached = dict(pdz)["ствол"]
    d = max(float(np.max(np.abs(cached[e] - fresh[e]))) for e in range(len(CYPS)))
    if d != 0.0:
        raise SystemExit(f"сид {s}: кэш членов расходится со стволом на {d:.3e}; пересоберите")
    print(f"      члены из кэша (контроль: ствол совпал побитово, {d:.1e})", flush=True)
    return pdz


def load_emb(rows):
    """The embedding, realigned against data/rows.csv rather than trusted.

    The equality test is shown to be CAPABLE of failing by running it again on the same array
    rolled by one position -- a comparison that cannot fail is not a control (CLAUDE.md).
    """
    z = np.load(EMB, allow_pickle=True)
    E, nm = z["train"].astype(np.float64), np.asarray(z["molecule_name"], object)
    r = rows.Molecule_Name.to_numpy(object)
    ok = bool((nm == r).all())
    rolled = bool((np.roll(nm, 1) == r).all())
    print(f"  вложение {E.shape}, порядок совпал: {ok}; сдвинутый контроль (ожидается False): "
          f"{rolled}", flush=True)
    if not ok or rolled:
        raise SystemExit("порядок вложения не совпал с rows.csv (или контроль не сработал)")
    if not np.isfinite(E).all():
        raise SystemExit("во вложении есть не-конечные значения")
    return E


def errcorr(probe, parts, y, mask):
    """Item 289's screen: Pearson correlation of out-of-fold ERROR.

    Two references, because the gate and the task ask different questions: item 289's own
    number is against the FIVE-MEMBER MEAN's error (verify/k88_family.py), and the competitor's
    0.67-0.79 was against each other member individually.
    """
    per_member, vs_mean, vs_ship = {k: [] for k in MEMBERS}, [], []
    calib = {k: [] for k in MEMBERS}
    for e in range(len(CYPS)):
        m = mask[:, e]
        yy = y[m, e]
        M = {k: P[e] for k, P in parts}
        five = np.mean([M[k] for k in MEMBERS], axis=0)
        ship = np.mean([M[k] for k in MEMBERS if _keep(e, k)], axis=0)
        ep = probe[e] - yy
        for k in MEMBERS:
            per_member[k].append(float(np.corrcoef(ep, M[k] - yy)[0, 1]))
            calib[k].append(float(np.corrcoef(M[k] - yy, five - yy)[0, 1]))
        vs_mean.append(float(np.corrcoef(ep, five - yy)[0, 1]))
        vs_ship.append(float(np.corrcoef(ep, ship - yy)[0, 1]))
    return per_member, vs_mean, vs_ship, calib


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0,1,2,3")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

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
    E = load_emb(rows)
    print(f"  состав по ферментам: " + ", ".join(
        f"{c}={'+'.join(SOLO[c]) if c in SOLO else 'все пять'}" for c in CYPS), flush=True)

    out = json.load(open(LOG)) if _pl.Path(LOG).exists() else {}
    out.setdefault("meta", {}).update(alphas=list(ALPHAS), gate=GATE, floor=FLOOR,
                                      wgrid=[float(x) for x in WGRID], emb=EMB)
    out.setdefault("seeds", {})

    for s in seeds:
        t0 = time.time()
        fold, ncl = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        print(f"\n  сид {s}: дайджест {dg}, кластеров {ncl}", flush=True)
        if s == 0 and dg != "2d93c19815e14261":
            raise SystemExit(f"разбиение сдвинулось: сид 0 дал {dg}")

        probe, dzprobe, alph, ncols = [], [], [], []
        for e in range(len(CYPS)):
            p, pk, nc = probe_oof(E, y, mask, fold, e)
            probe.append(p)
            dzprobe.append(probe_dz(E, y, mask, fold, e, p, LO, HI))
            alph.append(pk)
            ncols.append(nc)
        print(f"      alpha по фолдам: " + "; ".join(
            f"{CYPS[e]} {[int(x) for x in alph[e]]}" for e in range(len(CYPS))), flush=True)
        print(f"      живых колонок: {sorted(set(sum(ncols, [])))} из {E.shape[1]}", flush=True)

        pdz = members(s, XFULL, y, mask, fold, LO, HI)

        res = out["seeds"].get(str(s), {})
        res["A"] = score(combine_keep(pdz), y, LO, HI, mask, fold)
        res["D"] = score(combine_add(pdz, probe), y, LO, HI, mask, fold)
        res["Dz"] = score(combine_add(pdz, dzprobe), y, LO, HI, mask, fold)
        W, chosen = combine_fitted(pdz, probe, y, mask, fold)
        res["W"] = score(W, y, LO, HI, mask, fold)
        res["solo_probe"] = score(probe, y, LO, HI, mask, fold)
        res["solo_probe_dz"] = score(dzprobe, y, LO, HI, mask, fold)
        for k, P in pdz:
            res[f"solo_{k}"] = score(P, y, LO, HI, mask, fold)
        res["weights"] = chosen
        res["alphas"] = alph

        pm, vm, vs, cal = errcorr(probe, pdz, y, mask)
        res["errcorr"] = {"per_member": pm, "vs_five_mean": vm, "vs_shipped": vs,
                          "calibration_members_vs_five": cal}

        # Controls that must succeed AND exercise the suspect clause ------------------------
        # 1. the probe is not a constant and is not a copy of an existing member;
        # 2. the fitted-weight grid REACHES a>0 and the chooser is not merely declining;
        # 3. a=1.0 reproduces the probe exactly, i.e. the grid can express "all probe".
        dmin = min(float(np.min(np.abs(probe[e] - dict(pdz)[k][e]))) for e in range(len(CYPS))
                   for k in MEMBERS)
        spread = [float(np.std(probe[e])) for e in range(len(CYPS))]
        npos = sum(1 for c in chosen if c["a"] > 0)
        amax = max(c["a"] for c in chosen)
        one = [np.mean([P[e] for k, P in pdz if _keep(e, k)], axis=0) for e in range(len(CYPS))]
        blend1 = [(1 - 1.0) * one[e] + 1.0 * probe[e] for e in range(len(CYPS))]
        d1 = max(float(np.max(np.abs(blend1[e] - probe[e]))) for e in range(len(CYPS)))
        res["controls"] = {"min_abs_probe_minus_member": dmin, "probe_sd": spread,
                           "weight_cells_positive": npos, "weight_cells": len(chosen),
                           "max_a": amax, "a1_reproduces_probe_maxabs": d1,
                           "probe_range": [[float(np.min(p)), float(np.max(p))] for p in probe]}
        if min(spread) < 1e-6:
            raise SystemExit("контроль не сработал: зонд постоянен на каком-то ферменте")
        if d1 != 0.0:
            raise SystemExit("контроль не сработал: a=1 не воспроизводит зонд")
        print(f"      контроль: sd зонда {['%.3f' % v for v in spread]}, "
              f"a>0 в {npos}/{len(chosen)} ячейках, max a {amax}, a=1 даёт зонд ({d1:.1e})",
              flush=True)

        out["seeds"][str(s)] = res
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
        print(f"      сид {s}: {time.time()-t0:.0f} c, записано в {LOG}", flush=True)

    report(out)


def report(out):
    ks = sorted(out["seeds"], key=int)
    ks = [k for k in ks if all(x in out["seeds"][k] for x in ("A", "D", "Dz", "W"))]
    if not ks:
        return
    print(f"\n\n  === сидов: {len(ks)} ({', '.join(ks)}) ===", flush=True)

    print(f"\n  ЗОНД САМ ПО СЕБЕ (анатомия)", flush=True)
    print(f"  {'':>28s} " + " ".join(f"{c:>8s}" for c in CYPS) + f" {'макро':>8s}", flush=True)
    for nm, lab in [("solo_probe", "зонд chemprop_medium"),
                    ("solo_probe_dz", "зонд + мёртвая зона"),
                    ("solo_поферментно", "поферментно"), ("solo_пул", "пул"),
                    ("solo_GP", "GP"), ("solo_гребневая", "гребневая"),
                    ("solo_ствол", "ствол")]:
        kk = [k for k in ks if nm in out["seeds"][k]]
        if not kk:
            continue
        per = np.array([out["seeds"][k][nm][2] for k in kk]).mean(0)
        mac = np.array([out["seeds"][k][nm][0] for k in kk]).mean()
        print(f"  {lab:>28s} " + " ".join(f"{v:8.4f}" for v in per) + f" {mac:8.4f}", flush=True)

    ec = out["seeds"][ks[0]]["errcorr"]
    print(f"\n  КОРРЕЛЯЦИЯ ОШИБКИ (пункт 289, ворота 0.93), сид {ks[0]}, Пирсон", flush=True)
    print(f"  {'':>28s} " + " ".join(f"{c:>8s}" for c in CYPS) + f" {'сред.':>8s}", flush=True)
    for k in MEMBERS:
        v = np.array(ec["per_member"][k])
        print(f"  {'зонд против ' + k:>28s} " + " ".join(f"{x:8.3f}" for x in v)
              + f" {v.mean():8.3f}", flush=True)
    for key, lab in [("vs_five_mean", "зонд против средн. пятёрки"),
                     ("vs_shipped", "зонд против подаваемого")]:
        v = np.array(ec[key])
        print(f"  {lab:>28s} " + " ".join(f"{x:8.3f}" for x in v) + f" {v.mean():8.3f}",
              flush=True)
    print(f"  --- калибровка: наши собственные члены против средн. пятёрки ---", flush=True)
    for k in MEMBERS:
        v = np.array(ec["calibration_members_vs_five"][k])
        print(f"  {k:>28s} " + " ".join(f"{x:8.3f}" for x in v) + f" {v.mean():8.3f}",
              flush=True)

    print(f"\n  АНСАМБЛЬ (решающее)", flush=True)
    print(f"  {'арм':>34s} {'ранг ср.':>9s} {'sd':>8s} {'ST-RAE ср.':>11s}", flush=True)
    for nm, lab in [("A", "A подаваемый состав"),
                    ("D", "D + зонд шестым, равный вес"),
                    ("Dz", "Dz + зонд (мёрт. зона), равный вес"),
                    ("W", "W + зонд, подобранный вес")]:
        r = np.array([out["seeds"][k][nm][0] for k in ks])
        p = np.array([out["seeds"][k][nm][1] for k in ks])
        print(f"  {lab:>34s} {r.mean():9.4f} {r.std(ddof=1):8.4f} {p.mean():11.4f}", flush=True)

    print(f"\n  {'разность':>34s} {'Δранг ср.':>10s} {'sd':>8s} {'знак>0':>8s} "
          f"{'Δпара ср.':>10s}  ворота {GATE:+.4f} / пол {FLOOR:.4f}", flush=True)
    verdicts = {}
    for lab, arm in [("D - A  (РЕШАЮЩЕЕ, равный вес)", "D"),
                     ("Dz - A (чувствительность)", "Dz"),
                     ("W - A  (рука пункта 316)", "W")]:
        d = np.array([out["seeds"][k][arm][0] - out["seeds"][k]["A"][0] for k in ks])
        q = np.array([out["seeds"][k][arm][1] - out["seeds"][k]["A"][1] for k in ks])
        sd = d.std(ddof=1) if len(d) > 1 else 0.0
        sgn = int((d > 0).sum())
        if d.mean() >= FLOOR and sgn >= 3 and len(d) >= 4:
            v = "ПРОХОДИТ"
        elif d.mean() >= GATE and sgn >= 3 and len(d) >= 4:
            v = "ВНУТРИ ПОЛА"
        else:
            v = "не проходит"
        verdicts[arm] = (float(d.mean()), float(sd), sgn, v)
        print(f"  {lab:>34s} {d.mean():+10.4f} {sd:8.4f} {sgn:5d}/{len(d):<2d} "
              f"{q.mean():+10.4f}  {v}", flush=True)

    # per-enzyme anatomy of the deciding arm
    print(f"\n  D - A по ферментам (Δранг, среднее по сидам)", flush=True)
    per = np.array([[out["seeds"][k]["D"][2][e] - out["seeds"][k]["A"][2][e]
                     for e in range(len(CYPS))] for k in ks]).mean(0)
    print("  " + " ".join(f"{c} {v:+.4f}" for c, v in zip(CYPS, per)), flush=True)

    npos = sum(out["seeds"][k]["controls"]["weight_cells_positive"] for k in ks)
    ncell = sum(out["seeds"][k]["controls"]["weight_cells"] for k in ks)
    amax = max(out["seeds"][k]["controls"]["max_a"] for k in ks)
    print(f"\n  контроль руки 316: сетка весов {WGRID[0]}..{WGRID[-1]}, выбрано a>0 в "
          f"{npos}/{ncell} ячейках (фермент, фолд, сид), максимум a = {amax}", flush=True)

    print(f"\n  ПЕРЕВОД В МЕСТА НА ДОСКЕ. Одно место = {PLACE_STRAE:.4f} макро ST-RAE, "
          f"прочитано с доски 2026-09-25 (мы 0.5535, выше нас 0.5485).", flush=True)
    print(f"  {'арм':>34s} " + " ".join(f"{s:>16s}" for s in
          ("-1.62 причинный", "-0.97 ранги 43-73", "-0.58 ранги 48-68")), flush=True)
    for arm, lab in [("D", "D - A равный вес"), ("Dz", "Dz - A"), ("W", "W - A")]:
        dr = verdicts[arm][0]
        cells = []
        for slope in (-1.62, -0.97, -0.58):
            dstrae = slope * dr                      # ST-RAE moves DOWN when rank goes up
            cells.append(f"{-dstrae / PLACE_STRAE:+16.2f}")
        print(f"  {lab:>34s} " + " ".join(cells), flush=True)
    print(f"  (положительное = мест выиграно. Пункт 323: +0.0031 --- НИЖНЯЯ ГРАНИЦА цены "
          f"места, а не цена; окна у нашего ранга самые пологие на доске.)", flush=True)


if __name__ == "__main__":
    main()
