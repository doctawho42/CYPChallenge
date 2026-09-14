"""Docking contrast over the ensemble. Item 298 (pre-registered, blind).

This is the ablation item 298 promised, and it is written while `data/dock.npz` DOES NOT YET
EXIST -- the campaign is still running. That is the point: every choice below, including the
missing-value policy, is fixed before any affinity can be looked at. Item 298 fixed the four
acceptance conditions; this file fixes the arithmetic that feeds them.

Structure mirrors `verify/k93_overlay.py` exactly, because the comparison against the cheap
overlay proxy (item 296: +0.0017 on CYP1A2, a quarter of the floor) is the whole point of having
spent ~204 CPU-hours. Same base matrix, same cached members, same folds, same dead-zone pass, same
two cost economies. Only the block changes: an overlay onto one reference ligand becomes a docked
pose in the cavity itself.

THE FEATURE IS THE CONTRAST, not the four affinities. Row-centring across the four cavities cancels
size and lipophilicity -- which the ligand block already carries in 2295 columns -- and leaves
complementarity to a particular pocket. The contrast is therefore also the wrong-isoform control
built into the feature itself: item 296's CYP2C9 and CYP3A4 columns collapsed from +0.077 and +0.082
to +0.002 under centring, because there the raw score was mostly bulk.

THE ONE PLACE THIS CANNOT BE A COPY, and it is pre-registered here rather than decided later.
`data/overlay.npz` is dense: measured 0 NaN over 4905x4. `data/dock.npz` cannot be. Three of the
4905 rows have no 3D structure at all (`data/lig3d_train.sdf` holds 4902 records), and any docking
run may fail. Centring is a row operation, so a single NaN in a row makes all four of its columns
NaN and would poison `_oof_one`.

    Policy, fixed blind: centre first, then set every remaining NaN to 0.0.

Zero after centring means "no preference among the four pockets" -- the only filling that does not
invent an affinity for a molecule we never docked. Two gates guard it, and both can actually fire:
the imputed-row count is always printed, and a run where more than MAX_IMPUTED_FRAC of rows are
imputed REFUSES to print a verdict, because past that point the column is mostly fabrication and a
pass would be meaningless.

LIVE CELLS ARE TWO. After items 282-285 the per-enzyme member is kept only on CYP1A2 and CYP2D6.
CYP2C9 and CYP3A4 ship GP+ствол, so this arm cannot move them by construction -- they are a harness
control that must read exactly 0.000000 (item 298, condition 4).

Arms, identical in shape to k93 so the two campaigns are paired:
    A  ЦЕЛЕВАЯ             контраст СВОЕЙ полости, одна колонка
    B  НЕВЕРНАЯ ИЗОФОРМА   контраст ЧУЖОЙ полости (сдвиг на 1) -- what separates real pocket
                           complementarity from a volume descriptor any cavity would supply
    C  НЕНАПРАВЛЕННАЯ      все четыре колонки контраста каждому ферменту

The third falsifiable statement of item 298 is scored here too, because a prediction that no script
computes is a prediction that quietly disappears: CYP3A4's RAW affinities should correlate with
heavy-atom count more strongly than the other three cavities' do, reproducing item 199's +0.283
volume confound from the docking side. It needs no model fit, so it is reported even when the
ablation is skipped for want of member caches.

Reads data/dock.npz (built by `uv run python src/dock.py --assemble-only`), data/feats.npz,
results/preds/members_seed{4..7}.json. Writes results/logs/k97_dock.json.
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

from cypsplit import butina_folds
import submit as S
from submit import CYPS, _oof_one, _dz_oof, _keep

MEM = ["поферментно", "пул", "GP", "гребневая", "ствол"]
FLOOR = [0.0061, 0.0071, 0.0049, 0.0033]
# Refuse a verdict past this share of imputed rows: beyond it the column is mostly fabrication.
MAX_IMPUTED_FRAC = 0.01


def contrast(S_raw):
    """Per-molecule centring across the four cavities: removes size/greasiness, keeps pocket fit."""
    return S_raw - S_raw.mean(axis=1, keepdims=True)


def load_contrast(path):
    """Centre, then impute. Returns (C, n_imputed, raw). Policy pre-registered in the docstring."""
    z = np.load(path, allow_pickle=True)
    raw = z["train"].astype(np.float64)
    C = contrast(raw)
    bad = ~np.isfinite(C).all(axis=1)
    C[bad] = 0.0
    C[~np.isfinite(C)] = 0.0
    return C, int(bad.sum()), raw


def heavy_atoms(smiles):
    """Heavy-atom count per row, for item 298's volume-confound prediction."""
    from rdkit import Chem, RDLogger
    RDLogger.DisableLog("rdApp.*")
    out = np.full(len(smiles), np.nan)
    for i, s in enumerate(smiles):
        m = Chem.MolFromSmiles(s)
        if m is not None:
            out[i] = m.GetNumHeavyAtoms()
    return out


def volume_confound(raw, smiles):
    """Spearman of RAW affinity against heavy-atom count, per cavity. No model fit involved."""
    ha = heavy_atoms(smiles)
    out = {}
    for e, c in enumerate(CYPS):
        m = np.isfinite(raw[:, e]) & np.isfinite(ha)
        out[c] = float(spearmanr(raw[m, e], ha[m]).statistic) if m.sum() > 2 else float("nan")
    return out


def per_enzyme_member(Xs, y, LO, HI, mask, fold):
    """Dead-zone-passed per-enzyme member on an augmented matrix. Mirrors k93.per_enzyme_member."""
    P = _oof_one(Xs, y, mask, fold, False)
    out = []
    for e in range(len(CYPS)):
        m = mask[:, e]
        t = np.clip(P[e], LO[m, e], HI[m, e])
        out.append(_dz_oof("поферментно", Xs[m], t, fold[m], e))
    return out


def compose(members, e):
    keep = [np.asarray(members[k][e], float) for k in MEM if _keep(e, k)]
    return np.mean(keep, axis=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="4,5,6,7",
                    help="сиды с готовым кэшем членов; свежие для ЭТОЙ гипотезы (пункт 298)")
    a = ap.parse_args()
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    src = D + "dock.npz"
    if not _pl.Path(src).exists():
        print(f"нет {src}. Кампания ещё идёт либо блок не собран:\n"
              f"  uv run python src/dock.py --assemble-only", flush=True)
        return

    rows = pd.read_csv(D + "rows.csv")
    C, n_imp, raw = load_contrast(src)
    frac = n_imp / len(C)
    print(f"контраст докинга: {C.shape}, импутировано строк {n_imp} ({frac:.3%})", flush=True)
    if C.shape[0] != len(rows):
        print(f"  СТОП: {C.shape[0]} строк против {len(rows)} в rows.csv --- блок не выровнен",
              flush=True)
        return

    # Shape equality cannot tell an aligned block from a permuted or shifted one, and the positional
    # anchor that can is already on disk: data/lig3d_index.npz records exactly which rows of
    # rows.csv received a 3D structure, so the rows carrying a finite affinity must be precisely
    # that set. This check is here because its absence cost a campaign: src/dock.py was harvesting
    # 2270 of 4902 poses per cavity without a word, and this line reports it as a mismatch of 2632
    # rows instead of letting it arrive disguised as a 54 per cent imputation rate.
    ix = D + "lig3d_index.npz"
    if _pl.Path(ix).exists():
        zi = np.load(ix, allow_pickle=True)
        expect = set(int(v) for v in np.asarray(zi["ok_train"]).ravel().tolist())
        got = set(int(v) for v in np.flatnonzero(np.isfinite(raw).all(axis=1)).tolist())
        if got != expect:
            print(f"  СТОП: конечные строки блока не совпадают с lig3d_index.ok_train --- "
                  f"в блоке {len(got)}, ожидалось {len(expect)}, лишних {len(got - expect)}, "
                  f"недостающих {len(expect - got)}", flush=True)
            return
        print(f"  выравнивание сверено с lig3d_index: {len(got)} строк со структурой", flush=True)
    else:
        print(f"  ВНИМАНИЕ: нет {ix}, выравнивание проверено только по форме", flush=True)

    # Item 298's third prediction, scored whether or not the ablation runs below.
    vol = volume_confound(raw, list(rows.SMILES))
    print("\n  сырое сродство против числа тяжёлых атомов (пункт 199 дал +0.283):", flush=True)
    for c in CYPS:
        print(f"    {c:>8s} {vol[c]:+7.3f}", flush=True)
    worst = max(vol, key=lambda k: abs(vol[k]) if np.isfinite(vol[k]) else -1)
    print(f"    сильнее всех: {worst} -> предсказание пункта 298 "
          f"{'подтверждается' if worst == 'CYP3A4' else 'НЕ подтверждается'}", flush=True)

    if frac > MAX_IMPUTED_FRAC:
        print(f"\n  СТОП: импутировано {frac:.3%} строк при пороге {MAX_IMPUTED_FRAC:.3%}. "
              f"Вердикт не печатаю --- колонка в основном выдумана.", flush=True)
        return

    z = np.load(D + "feats.npz")
    base = np.hstack([z["FP"], z["DESC"], z["MECH"]])

    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    mask = ~np.isnan(y)
    S.LO, S.HI = LO, HI

    # Arm A gives enzyme e its OWN cavity column; arm B gives it the NEXT cavity's (wrong isoform);
    # arm C gives every enzyme all four. A and B differ in exactly one column, which is the point.
    ARMS = {"A целевая": [C[:, [e]] for e in range(4)],
            "B неверная изоформа": [C[:, [(e + 1) % 4]] for e in range(4)],
            "C ненаправленная": [C for _ in range(4)]}

    dst = RES + "logs/k97_dock.json"
    out = json.load(open(dst)) if _pl.Path(dst).exists() else {"seeds": {}}
    out["импутировано"] = {"строк": n_imp, "доля": frac}
    out["объёмный конфаунд"] = vol
    print(f"\n  {'сид':>3s} {'арм':>20s} {'фермент':>8s} {'подаётся':>9s} {'+докинг':>9s} "
          f"{'Δранг':>8s}", flush=True)
    for s in seeds:
        cache = RES + f"preds/members_seed{s}.json"
        if not _pl.Path(cache).exists():
            print(f"  сид {s}: нет {cache}, пропуск", flush=True)
            continue
        mem = json.load(open(cache))
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        res = out["seeds"].get(str(s), {})

        ref = {c: float(spearmanr(y[mask[:, e], e], compose(mem, e)).statistic)
               for e, c in enumerate(CYPS)}
        res["эталон"] = ref

        for arm, blocks in ARMS.items():
            t0 = time.time()
            # Both economies are k93's and are kept verbatim, because condition 4 depends on the
            # second one: where the shipped SOLO does not keep the per-enzyme member, compose()
            # ignores the rebuilt member entirely, so the delta is 0.000000 BY CONSTRUCTION rather
            # than by measurement. CYP2C9 and CYP3A4 ship GP+ствол, so this also halves arms A and B.
            aug = [np.asarray(mem["поферментно"][e], float) for e in range(len(CYPS))]
            shared = all(b is blocks[0] for b in blocks)
            if shared:
                Xs = np.hstack([base, blocks[0]])
                P = _oof_one(Xs, y, mask, fold, False)
                for e in range(len(CYPS)):
                    if not _keep(e, "поферментно"):
                        continue
                    m = mask[:, e]
                    t = np.clip(P[e], LO[m, e], HI[m, e])
                    aug[e] = _dz_oof("поферментно", Xs[m], t, fold[m], e)
            else:
                for e in range(len(CYPS)):
                    if not _keep(e, "поферментно"):
                        continue
                    Xs = np.hstack([base, blocks[e]])
                    one = np.zeros_like(mask); one[:, e] = mask[:, e]
                    P = _oof_one(Xs, y, one, fold, False)
                    m = mask[:, e]
                    t = np.clip(P[e], LO[m, e], HI[m, e])
                    aug[e] = _dz_oof("поферментно", Xs[m], t, fold[m], e)
            memA = dict(mem); memA["поферментно"] = [np.asarray(v).tolist() for v in aug]
            res[arm] = {}
            for e, c in enumerate(CYPS):
                r1 = float(spearmanr(y[mask[:, e], e], compose(memA, e)).statistic)
                res[arm][c] = r1
                print(f"  {s:3d} {arm:>20s} {c:>8s} {ref[c]:9.4f} {r1:9.4f} {r1-ref[c]:+8.4f}",
                      flush=True)
            print(f"      ({time.time()-t0:.0f} c)", flush=True)

        out["seeds"][str(s)] = res
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(dst, "w"), ensure_ascii=False, indent=1)

    # ---- verdict against item 298 -----------------------------------------------------------
    ks = [k for k in out["seeds"] if "эталон" in out["seeds"][k]]
    if not ks:
        return
    print(f"\n  {'арм':>20s} {'фермент':>8s} {'Δранг ср.':>10s} {'sd':>8s} {'знак':>6s} "
          f"{'пол':>8s} {'вердикт':>10s}", flush=True)
    summ = {}
    for arm in ARMS:
        summ[arm] = {}
        for e, c in enumerate(CYPS):
            d = np.array([out["seeds"][k][arm][c] - out["seeds"][k]["эталон"][c] for k in ks
                          if arm in out["seeds"][k]])
            if not len(d):
                continue
            # Item 298 fixed "sign holds on at least 3 of the 4 seeds". k93's 0.75*len(d) form is
            # identical at len(d)==4 (0.75*4 = 3.0, needs 3) but silently becomes 2-of-2 at
            # len(d)==2 and 1-of-1 at len(d)==1 -- so a short --seeds list or one lost member cache
            # would adopt a 204-CPU-hour result on fewer seeds than were pre-registered. Require
            # the counts outright. At four seeds this changes no number.
            ok = len(d) >= 3 and d.mean() > FLOOR[e] and (d > 0).sum() >= 3
            summ[arm][c] = {"mean": float(d.mean()),
                            "sd": float(d.std(ddof=1)) if len(d) > 1 else 0.0,
                            "sign": int((d > 0).sum()), "n": int(len(d)), "passes": bool(ok)}
            print(f"  {arm:>20s} {c:>8s} {d.mean():+10.4f} "
                  f"{d.std(ddof=1) if len(d) > 1 else 0:8.4f} {int((d>0).sum())}/{len(d):<2d} "
                  f"{FLOOR[e]:8.4f} {'ПРОХОДИТ' if ok else 'нет':>10s}", flush=True)

    print("\n  === условия пункта 298 ===", flush=True)
    for c in ("CYP1A2", "CYP2D6"):
        A = summ.get("A целевая", {}).get(c)
        B = summ.get("B неверная изоформа", {}).get(c)
        if not A or not B:
            continue
        beats = A["mean"] > B["mean"]
        print(f"  {c}: A {A['mean']:+.4f} (знак {A['sign']}/{A['n']}) против "
              f"B {B['mean']:+.4f} -> целевая {'бьёт' if beats else 'НЕ бьёт'} неверную изоформу; "
              f"порог {'взят' if A['passes'] else 'не взят'}   "
              f"=> {'ПРИНЯТО' if (A['passes'] and beats) else 'ОТКЛОНЕНО'}", flush=True)
    for c in ("CYP2C9", "CYP3A4"):
        A = summ.get("A целевая", {}).get(c)
        if A:
            print(f"  контроль стенда {c}: Δ {A['mean']:+.6f} "
                  f"(обязан быть 0.000000 --- поферментный член туда не входит)", flush=True)

    # Against the cheap proxy this campaign was supposed to beat (item 296: +0.0017 on CYP1A2).
    prev = RES + "logs/k93_overlay.json"
    A1 = summ.get("A целевая", {}).get("CYP1A2")
    if _pl.Path(prev).exists() and A1:
        p = json.load(open(prev)).get("итог", {}).get("A целевая", {}).get("CYP1A2")
        if p:
            print(f"\n  против дешёвого прокси (пункт 296): перекрытие {p['mean']:+.4f}, "
                  f"докинг {A1['mean']:+.4f} -> докинг "
                  f"{'бьёт' if A1['mean'] > p['mean'] else 'НЕ бьёт'} перекрытие", flush=True)

    out["итог"] = summ
    json.dump(out, open(dst, "w"), ensure_ascii=False, indent=1)
    print(f"\n  сохранено: {dst}", flush=True)


if __name__ == "__main__":
    main()
