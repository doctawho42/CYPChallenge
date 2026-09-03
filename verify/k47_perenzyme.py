"""Re-read every saved ablation per enzyme, because the macro floor buried the per-enzyme effects.

Why this is not a new experiment. Item 165 measured what had never been written down: the noise
floor is **per enzyme** 0.0061, 0.0071, 0.0049 and 0.0033, against 0.0036 on the macro. The macro
averages four enzymes and is therefore quieter than three of them -- and this file has spent a
hundred and sixty items judging per-enzyme claims against the macro figure of 0.007.

The arithmetic that follows is the same one that hid the NCGC panel's CYP1A2 effect until item 157
split it out: **an intervention that helps one enzyme by 0.04 and does nothing to the other three
shows up on the macro as 0.010**, which sits at the floor and reads as a null. Item 81 is the proof
that this happens here -- the mechanistic block is +0.0454 on CYP2D6 and -0.0019 on CYP1A2, and its
macro is +0.0148.

So before building anything new, every ablation whose predictions are on disk gets read again,
per enzyme, against that enzyme's own floor. The blocks were built for chemistry that is
enzyme-specific by construction:

    PBF, планарность         узкая плоская щель CYP1A2
    объём, SASA              большая полость CYP3A4
    кислоты, анион при 7.4   Arg108 в CYP2C9
    геометрия основного N    Glu216 в CYP2D6 --- единственная, что строилась намеренно

None of these has any reason to help all four enzymes, and every one of them was scored on the
macro.

This computes nothing new. It reads saved out-of-fold predictions, pairs each arm against its
file's baseline, and reports the per-enzyme difference with the spread across seeds. Runs in
under a minute.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES

import glob
import json
import os

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
# Пол по ферментам из пункта 165: sd ранга по сидам при неизменной руке.
FLOOR = {"CYP1A2": 0.0061, "CYP2C9": 0.0071, "CYP2D6": 0.0049, "CYP3A4": 0.0033}
BASE_HINT = ("без ", "база", "независимо", "L2 по метке", "эталон", "FP+DESC+MECH")


def main():
    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    truth = {}
    for c in CYPS:
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        truth[c] = (tr.loc[m, col].to_numpy(), int(m.sum()))

    print(f"пол по ферментам (пункт 165): " + "  ".join(f"{c[3:]} {FLOOR[c]:.4f}" for c in CYPS))
    print(f"пол по макро: 0.0036 --- он МЕНЬШЕ трёх из четырёх\n")

    out = []
    for path in sorted(glob.glob(RES + "preds/oof_*.json")):
        fn = os.path.basename(path)
        try:
            blob = json.load(open(path))
        except Exception:
            continue
        P = blob.get("preds", blob)
        if not isinstance(P, dict):
            continue
        # ключи бывают «рука|фермент» и «сид|рука|фермент»
        recs = {}
        for k, v in P.items():
            parts = k.split("|")
            if parts[-1] not in CYPS:
                continue
            c = parts[-1]
            seed = parts[0] if parts[0].isdigit() else "0"
            arm = "|".join(parts[1:-1]) if parts[0].isdigit() else "|".join(parts[:-1])
            if len(v) != truth[c][1]:
                continue
            recs.setdefault(arm, {}).setdefault(seed, {})[c] = np.asarray(v, float)
        arms = list(recs)
        if len(arms) < 2:
            continue
        base = next((a for a in arms if any(h in a for h in BASE_HINT)), None)
        if base is None:
            continue
        for arm in arms:
            if arm == base:
                continue
            seeds = sorted(set(recs[arm]) & set(recs[base]))
            if not seeds:
                continue
            d = {c: [] for c in CYPS}
            for s in seeds:
                for c in CYPS:
                    if c in recs[arm][s] and c in recs[base][s]:
                        y = truth[c][0]
                        d[c].append(spearmanr(y, recs[arm][s][c]).statistic
                                    - spearmanr(y, recs[base][s][c]).statistic)
            if not all(d[c] for c in CYPS):
                continue
            r = {"файл": fn.replace("oof_", "").replace(".json", ""), "рука": arm[:34],
                 "сидов": len(seeds)}
            for c in CYPS:
                r[c[3:]] = float(np.mean(d[c]))
            r["МАКРО"] = float(np.mean([r[c[3:]] for c in CYPS]))
            # сколько ферментов вышли за СВОЙ пол, и куда
            r["выше пола"] = sum(1 for c in CYPS if r[c[3:]] > FLOOR[c])
            r["ниже пола"] = sum(1 for c in CYPS if r[c[3:]] < -FLOOR[c])
            out.append(r)

    df = pd.DataFrame(out)
    if df.empty:
        print("нечего читать")
        return
    # Интересны те, у кого МАКРО у пола, а хоть один фермент --- заметно выше своего.
    df["скрыто макро"] = (df["МАКРО"].abs() < 0.010) & (df["выше пола"] >= 1)
    df = df.sort_values("скрыто макро", ascending=False)

    show = ["файл", "рука", "сидов", "МАКРО"] + [c[3:] for c in CYPS] + ["выше пола", "ниже пола"]
    hid = df[df["скрыто макро"]]
    print("=" * 100)
    print("РУКИ, У КОТОРЫХ МАКРО У ПОЛА, А ХОТЯ БЫ ОДИН ФЕРМЕНТ ВЫШЕ СВОЕГО")
    print("=" * 100)
    print(hid[show].round(4).to_string(index=False) if len(hid) else "  таких нет")
    print("\n" + "=" * 100)
    print("ВСЁ ОСТАЛЬНОЕ")
    print("=" * 100)
    print(df[~df["скрыто макро"]][show].round(4).to_string(index=False))
    df.to_csv(RES + "perenzyme.csv", index=False)
    print(f"\nсохранено: {RES}perenzyme.csv")
    print("""
Как читать. Верхняя таблица --- руки, которые файл закрыл по макро, а по своему ферменту они
за полом. Это НЕ автоматически находки: разность по одному ферменту шумнее макро, и часть
строк там окажется случайной. Но каждая такая строка --- это вопрос, который был закрыт
измерением не той величины, и его надо задать заново.

«выше пола» и «ниже пола» считают ферменты, у которых |разность| больше СОБСТВЕННОГО пола из
пункта 165, а не общего 0.007.

Рука, у которой один фермент выше своего пола, а другой ниже своего, --- не провал, а размен.
Файл до сих пор мог видеть такие только как ноль.""")


if __name__ == "__main__":
    main()
