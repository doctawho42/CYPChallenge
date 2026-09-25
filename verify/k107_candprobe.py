"""The candidate built with the sixth member (item 324): its controls, its placement, its gate.

Four questions, in the order that makes the later ones worth asking at all.

1. IS THE DEFAULT PATH STILL THE DEFAULT PATH. `src/submit.py --probe` writes
   activity_submission_noprobe.csv beside its output: the SAME run, the same member
   predictions, the composition WITHOUT the probe, and its own affine pair. That file must
   reproduce results/submission/activity_submission.csv -- the file on disk from before the
   probe existed, built at commit fde6b221 with the bundle ON -- element for element. This is
   a stronger statement than "the flag is off by default": it says the run that HAS the flag
   on still computes the old vector unchanged, and incidentally that --no-bundle does not
   touch the activity columns. A tolerance is not offered; the answer is 0.0 or the default
   path moved.

2. IS THE NEW MEMBER ACTUALLY A NEW MEMBER. A probe that had collapsed onto one of our five
   would leave the out-of-fold gain behind on the test set. So its test predictions are
   correlated against each of our members' test predictions, read from the run's own
   probe_parts.json. Near 1.0 refutes; the out-of-fold decorrelation was 0.868-0.900 against
   members that sit at 0.966-0.982 against each other. The range check is part of the same
   control: a probe predicting outside any plausible pIC50 range would be a broken fit that
   correlation alone could still call healthy.

3. WHAT SCALE DOES THE PLACED FILE NEED. `mu` is an absolute target mean and transfers from
   item 320's placement; `b` multiplies the spread of whatever vector it is given, so
   sd(q) = b*sd(p) and it does NOT. The target is the sd the board file HAS, read off that
   file rather than quoted, and b is re-derived against the new vector's spread. The
   placement is otherwise held fixed: this run is about the probe, and a placement change in
   the same file would confound the two.

4. DOES THE ORGANISERS' VALIDATOR ACCEPT IT. With expected_ids SET, and with five refuters
   that must each be REJECTED with a DIFFERENT message. One of them (R5) is the reason the
   others are not enough: it renames one molecule, so the file still has 750 rows, no
   duplicates and no non-finite values, and the validator ACCEPTS it when expected_ids is
   None. It is rejected only if the identity check is really running. A suite of refuters
   that fail either way would prove the validator rejects things, not that this call does.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import argparse, hashlib, json
import numpy as np
import pandas as pd
from scipy.stats import spearmanr
from validation.activity_validation import validate_activity_submission

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
COL = "{c}_pIC50_direct_inhibition"
MU_Y = [4.390, 4.851, 3.400, 4.815]     # item 320's placement, held fixed
SHIPPED = RES + "submission/activity_submission.csv"
BOARD = RES + "submission/activity_submission_cand2d6xc.csv"
LOG = RES + "logs/k107_candprobe.json"


def sha256(p):
    return hashlib.sha256(open(p, "rb").read()).hexdigest()


def col(df, c):
    return df[COL.format(c=c)].to_numpy(float)


def twin_control(twin_path):
    """Question 1. The probe-free twin against the file built before the probe existed."""
    a, b = pd.read_csv(twin_path), pd.read_csv(SHIPPED)
    if not (a.Molecule_Name.equals(b.Molecule_Name) and a.SMILES.equals(b.SMILES)):
        raise SystemExit("близнец и подаваемый файл различаются идентификаторами")
    out = {}
    for c in CYPS:
        d = np.abs(col(a, c) - col(b, c))
        out[c] = {"max_abs": float(d.max()), "n_exact": int((d == 0).sum()), "n": int(len(d))}
    # Non-vacuity: the same comparison against a vector that is KNOWN to differ must report a
    # difference. A check that reads 0.0 on everything proves nothing about what it compared.
    sham = np.abs(col(a, CYPS[0]) - col(b, CYPS[1]))
    out["контроль_несовпадения"] = {"max_abs": float(sham.max()),
                                    "замечание": "1A2 против 2C9, обязан быть > 0"}
    if sham.max() <= 0:
        raise SystemExit("контроль не сработал: разные столбцы сравнялись")
    return out


def decorrelation(parts_path):
    """Question 2. The probe against our own members, on the TEST predictions."""
    P = json.load(open(parts_path))
    per = {}
    for c in CYPS:
        blk = P["члены_по_ферментам"][c]
        pr = np.asarray(blk["зонд"], float)
        ours = {k: np.asarray(v, float) for k, v in blk["члены"].items()}
        mean_ours = np.mean(list(ours.values()), axis=0)
        per[c] = {
            "зонд": {"mean": float(pr.mean()), "sd": float(pr.std(ddof=1)),
                     "min": float(pr.min()), "max": float(pr.max())},
            "pearson_против_члена": {k: float(np.corrcoef(pr, v)[0, 1]) for k, v in ours.items()},
            "spearman_против_члена": {k: float(spearmanr(pr, v).statistic)
                                      for k, v in ours.items()},
            "pearson_против_среднего_наших": float(np.corrcoef(pr, mean_ours)[0, 1]),
            # Our own members against their mean, as the calibration for "what does a member
            # that is NOT decorrelated look like on this data".
            "калибровка_наши_против_среднего": {
                k: float(np.corrcoef(v, mean_ours)[0, 1]) for k, v in ours.items()},
            "членов": list(ours),
        }
        rng = per[c]["зонд"]
        if not (0.0 < rng["min"] and rng["max"] < 12.0):
            raise SystemExit(f"{c}: зонд вне правдоподобного диапазона pIC50 "
                             f"({rng['min']:.2f}..{rng['max']:.2f})")
        worst = max(per[c]["pearson_против_члена"].values())
        if worst > 0.99:
            raise SystemExit(f"{c}: зонд схлопнулся на существующий член (r = {worst:.4f})")
    return per


def derive_b(base_path):
    """Question 3. b so that the placed spread matches the board file's, mu held fixed."""
    new, board = pd.read_csv(base_path), pd.read_csv(BOARD)
    old_base = pd.read_csv(SHIPPED)
    rec = {}
    for e, c in enumerate(CYPS):
        tgt = float(col(board, c).std(ddof=1))
        sd_new = float(col(new, c).std(ddof=1))
        sd_old = float(col(old_base, c).std(ddof=1))
        b_old = tgt / sd_old
        b_new = tgt / sd_new
        # The whole transfer argument rests on the board file BEING recalib(shipped, b, mu).
        # That is checked rather than read off the meta: q rebuilt here from the shipped file
        # and item 320's own constants must reproduce the board file. If it does not, the
        # target sd below is the sd of some other file and none of this is about the board.
        p_old = col(old_base, c)
        q_old = MU_Y[e] + 1.0 * (p_old - p_old.mean()) * (tgt / sd_old)
        rec[c] = {"цель_sd": tgt, "sd_старого_вектора": sd_old, "sd_нового_вектора": sd_new,
                  "b_старое": b_old, "b_новое": b_new,
                  "b_новое_округл": round(b_new, 6),
                  "отношение": b_new / b_old, "mu_y": MU_Y[e],
                  "восстановлен_файл_доски_max_abs": float(np.abs(q_old - col(board, c)).max())}
        if rec[c]["восстановлен_файл_доски_max_abs"] > 1e-9:
            raise SystemExit(
                f"{c}: файл на доске НЕ воспроизводится как аффинное преобразование "
                f"подаваемого (расхождение "
                f"{rec[c]['восстановлен_файл_доски_max_abs']:.3e}) --- цель sd относится "
                f"не к тому файлу")
    return rec


REFUTERS = {
    "R1 нет столбца": lambda d: d.drop(columns=[COL.format(c="CYP2D6")]),
    "R2 дубликат имени": lambda d: pd.concat([d, d.iloc[[0]]], ignore_index=True),
    "R3 пропуск в значении": lambda d: d.assign(
        **{COL.format(c="CYP1A2"): lambda x: x[COL.format(c="CYP1A2")].mask(
            x.index == 3, np.nan)}),
    "R4 бесконечность": lambda d: d.assign(
        **{COL.format(c="CYP3A4"): lambda x: x[COL.format(c="CYP3A4")].mask(
            x.index == 7, np.inf)}),
    # The discriminating one: 750 rows, no duplicates, everything finite. Only the identity
    # check can catch it, so it is ACCEPTED when expected_ids is None and rejected when set.
    "R5 переименована молекула": lambda d: d.assign(
        Molecule_Name=d.Molecule_Name.mask(d.index == 11, "NOT-A-REAL-MOLECULE-324")),
}


def gate(cand_path, tmpdir):
    """Question 4. The organisers' validator on the candidate, then five refuters."""
    ids = set(pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv").Molecule_Name)
    ok, errs = validate_activity_submission(cand_path, expected_ids=ids)
    res = {"кандидат": {"принято": bool(ok), "ошибки": list(errs)},
           "expected_ids_из": "data/cyp-challenge-TEST-BLINDED.csv", "n_ids": len(ids),
           "опровергатели": {}}
    if not ok:
        raise SystemExit(f"валидатор отверг кандидата: {errs}")
    d0 = pd.read_csv(cand_path)
    _pl.Path(tmpdir).mkdir(parents=True, exist_ok=True)
    seen = {}
    for i, (name, f) in enumerate(REFUTERS.items()):
        p = str(_pl.Path(tmpdir) / f"refuter{i}.csv")
        f(d0.copy()).to_csv(p, index=False)
        rok, rerr = validate_activity_submission(p, expected_ids=ids)
        nok, _nerr = validate_activity_submission(p)      # the same file, clause disarmed
        res["опровергатели"][name] = {"принято_с_expected_ids": bool(rok),
                                      "принято_без_expected_ids": bool(nok),
                                      "первая_ошибка": (list(rerr) or [""])[0],
                                      "все_ошибки": list(rerr)}
        if rok:
            raise SystemExit(f"опровергатель {name} НЕ был отвергнут -- гейт не работает")
        msg = res["опровергатели"][name]["первая_ошибка"]
        if msg in seen:
            raise SystemExit(f"{name} и {seen[msg]} отвергнуты ОДНИМ сообщением: {msg}")
        seen[msg] = name
    r5 = res["опровергатели"]["R5 переименована молекула"]
    if not r5["принято_без_expected_ids"]:
        raise SystemExit("R5 отвергнут и без expected_ids -- он не проверяет этот аргумент, "
                         "и про сам expected_ids набор по-прежнему ничего не доказывает")
    return res


def compare(cand_path):
    """The candidate against the file currently on the board. This WILL move: the vector
    changed, unlike every placement-only candidate before it, and by how much is the point."""
    a, b = pd.read_csv(cand_path), pd.read_csv(BOARD)
    out = {"sha256_кандидата": sha256(cand_path), "sha256_доски": sha256(BOARD),
           "по_ферментам": {}}
    for c in CYPS:
        x, yb = col(a, c), col(b, c)
        out["по_ферментам"][c] = {
            "spearman_кандидат_против_доски": float(spearmanr(x, yb).statistic),
            "pearson": float(np.corrcoef(x, yb)[0, 1]),
            "кандидат": {"mean": float(x.mean()), "sd": float(x.std(ddof=1)),
                         "min": float(x.min()), "max": float(x.max())},
            "доска": {"mean": float(yb.mean()), "sd": float(yb.std(ddof=1)),
                      "min": float(yb.min()), "max": float(yb.max())},
            "сдвинутых_пар_ранга": int((np.argsort(np.argsort(x))
                                        != np.argsort(np.argsort(yb))).sum()),
        }
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--rundir", required=True, help="выход прогона src/submit.py --probe")
    ap.add_argument("--cand", default="", help="размещённый кандидат; без него --- только b")
    ap.add_argument("--tmp", default="/tmp/k107_refuters")
    a = ap.parse_args()
    rd = a.rundir if a.rundir.endswith("/") else a.rundir + "/"

    out = {}
    print("1. КОНТРОЛЬ-БЛИЗНЕЦ: состав без зонда против подаваемого файла")
    out["близнец"] = twin_control(rd + "activity_submission_noprobe.csv")
    for c in CYPS:
        v = out["близнец"][c]
        print(f"   {c:>8s} max|разность| {v['max_abs']:.3e}, побитово равных "
              f"{v['n_exact']}/{v['n']}")
    print(f"   контроль несовпадения (1A2 против 2C9): "
          f"{out['близнец']['контроль_несовпадения']['max_abs']:.3f} > 0")

    print("\n2. ДЕКОРРЕЛЯЦИЯ ЗОНДА НА ТЕСТЕ")
    out["декорреляция"] = decorrelation(rd + "probe_parts.json")
    for c in CYPS:
        v = out["декорреляция"][c]
        pm = v["pearson_против_члена"]
        cal = v["калибровка_наши_против_среднего"]
        print(f"   {c:>8s} зонд {v['зонд']['min']:.2f}..{v['зонд']['max']:.2f} "
              f"(ср. {v['зонд']['mean']:.2f}, sd {v['зонд']['sd']:.2f}); "
              f"r против членов " + ", ".join(f"{k} {x:.3f}" for k, x in pm.items()))
        print(f"   {'':>8s} наши члены против своего среднего: "
              + ", ".join(f"{k} {x:.3f}" for k, x in cal.items()))

    print("\n3. ПЕРЕВЫВОД b (mu удержано)")
    out["b"] = derive_b(rd + "activity_submission.csv")
    print(f"   {'фермент':>8s} {'цель sd':>8s} {'sd стар':>8s} {'sd нов':>8s} "
          f"{'b стар':>7s} {'b нов':>7s} {'отнош.':>7s}")
    for c in CYPS:
        v = out["b"][c]
        print(f"   {c:>8s} {v['цель_sd']:8.4f} {v['sd_старого_вектора']:8.4f} "
              f"{v['sd_нового_вектора']:8.4f} {v['b_старое']:7.3f} {v['b_новое']:7.4f} "
              f"{v['отношение']:7.3f}   (доска восстановлена с "
              f"{v['восстановлен_файл_доски_max_abs']:.1e})")
    print("   команда: uv run python src/recalib.py --b "
          + ",".join(f"{out['b'][c]['b_новое_округл']:.6f}" for c in CYPS)
          + " --mu-y " + ",".join(f"{m:.3f}" for m in MU_Y))

    if a.cand:
        print("\n4. ГЕЙТ ОРГАНИЗАТОРОВ И ОПРОВЕРГАТЕЛИ")
        out["гейт"] = gate(a.cand, a.tmp)
        print(f"   кандидат: принято (expected_ids из {out['гейт']['expected_ids_из']}, "
              f"{out['гейт']['n_ids']} имён)")
        for k, v in out["гейт"]["опровергатели"].items():
            print(f"   {k:>26s}: отвергнут; без expected_ids "
                  f"{'ПРИНЯТ' if v['принято_без_expected_ids'] else 'тоже отвергнут'}"
                  f"  <- {v['первая_ошибка'][:78]}")

        print("\n5. КАНДИДАТ ПРОТИВ ФАЙЛА НА ДОСКЕ")
        out["против_доски"] = compare(a.cand)
        print(f"   sha256 кандидата {out['против_доски']['sha256_кандидата']}")
        print(f"   {'фермент':>8s} {'spearman':>9s} {'ср. канд':>9s} {'sd канд':>8s} "
              f"{'ср. доска':>10s} {'sd доска':>9s}")
        for c in CYPS:
            v = out["против_доски"]["по_ферментам"][c]
            print(f"   {c:>8s} {v['spearman_кандидат_против_доски']:9.4f} "
                  f"{v['кандидат']['mean']:9.3f} {v['кандидат']['sd']:8.4f} "
                  f"{v['доска']['mean']:10.3f} {v['доска']['sd']:9.4f}")

    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)
    print(f"\nзаписано: {LOG}")


if __name__ == "__main__":
    main()
