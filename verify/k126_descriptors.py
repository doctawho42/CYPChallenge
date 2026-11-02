"""Do DESCRIPTOR-pretrained checkpoints fail because the target is a descriptor, or because the
descriptors overlap our own features? Three checkpoints chosen to have no overlap at all.

WHAT IS UNDER TEST. Item 349 refined item 345's criterion: across six measured cases the separating
variable is whether the pretraining target is an ADME/activity ENDPOINT (three successes, one of
them a pure distillation of another model's predictions) or a DESCRIPTOR (two failures). The clause
item 345 believed in -- "not already present in our features" -- was shown to do none of the work,
because chemberta-MTR satisfies it and fails at -0.0070.

But both descriptor failures are weak evidence for the strong claim. `rdkit2d` is pretrained on
exactly our DESC block, so it is the maximally overlapping case by construction, and
chemberta-MTR's targets are computed descriptors of the same 2D family. **Neither tests a descriptor
target that our features do not carry at all.** If overlap is what matters after all, a
non-overlapping descriptor checkpoint should work, and item 349's split is wrong.

THE THREE, chosen for minimal overlap with FP (2048 Morgan counts r=2), DESC (217 RDKit 2D) and
MECH (30 mechanistic):

    e3fp      1024 целей, 3D-отпечаток       3D у нас в отгруженной матрице НЕТ
    atompair  2048 целей, атомные пары       другое семейство, не круговые отпечатки Моргана
    usrcat      60 целей, 3D-форма           3D, и самая узкая цель из трёх

`e3fp` and `atompair` are deliberately WIDE (1024 and 2048 targets) so that a failure cannot be
explained by the narrow-bottleneck reasoning that item 349's prediction rested on -- `surrogate_adme`
has 25 targets and works, so width is not the variable, but a 60-target failure beside a
2048-target failure settles it inside this experiment rather than by cross-reference.

THE 3D PRECONDITION, checked before building rather than assumed. A 3D-descriptor checkpoint only
tests the overlap clause if 3D information is (a) absent from the shipped matrix and (b) not useless
for this task. Both hold: `data/shape3d.npz` (4905x16) appears zero times in the shipped feature
path, and items 119/165/279 measured those sixteen columns at **+0.0034 macro** -- just under the
0.0036 floor -- with **+0.0087 on CYP2D6 at sign 4/4**, its per-enzyme profile being
+0.0060 / -0.0003 / +0.0087 / -0.0008 over CYPS. So 3D carries real signal here, and little of it.

PREDICTION, fixed before any number for these three exists, as items 347 and 349 did.

    ГЛАВНОЕ: все три проваливаются -- одиночный ранг ниже порога 0.60 пункта 348 и вклад
    седьмым членом не выше нуля при знаке 0/4 -- потому что класс цели дескрипторный.

    ЧТО ЭТО ОПРОВЕРГНЕТ: любой из трёх, положительный при 4/4 выше пола, означает, что
    разделяет всё-таки ПЕРЕСЕЧЕНИЕ с нашими признаками, а не класс цели, и пункт 349 неверен.
    Из трёх я ставлю на e3fp: 1024 цели, чистое 3D, нулевое пересечение. usrcat при 60 целях
    наименее вероятен.

    ВТОРОЕ, по ферментам: если 3D вообще переносится, оно должно проявиться на CYP2D6, где
    шестнадцать столбцов формы стоят +0.0087 при 4/4. Положительная клетка 2D6 при
    отрицательном макро повторила бы блок формы и была бы поферментной находкой, а не
    опровержением -- и именно поэтому поферментные клетки печатаются здесь, а не только макро.

NO NEW DEVICE. Drives `verify/k121_chemeleon.py --emb/--tag/--out`, the same code path items 345,
348 and 349 used, including the alignment control these nameless embeddings need.

    uv run python verify/k126_descriptors.py --seeds 0,1,2,3

Reads data/emb_chemprop_{e3fp,atompair,usrcat}.npz from `src/embed.py`. Writes
results/logs/k126_{tag}.json. Exit: 0 ran, 2 a control failed, 3 input missing.
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

from submit import CYPS
import k121_chemeleon as K121

FLOOR = 0.0036
THRESHOLD = 0.60
CANDS = [
    ("e3fp", 1024, "3D-отпечаток, 3D в матрице НЕТ -> ставка на него, если пункт 349 неверен"),
    ("atompair", 2048, "атомные пары, не Морган -> широкая неперекрывающаяся 2D-цель"),
    ("usrcat", 60, "3D-форма, самая узкая цель -> наименее вероятен"),
]
# Item 279's per-enzyme profile for the sixteen shape columns, over CYPS, on the bare member.
SHAPE3D = [+0.0060, -0.0003, +0.0087, -0.0008]
# Item 349's classes, four-seed means on this device.
REF = [("эндпойнт (измеренный)", "chemprop_medium", 0.6261, +0.0138),
       ("эндпойнт (измеренный)", "chemprop_chemeleon", 0.6153, +0.0033),
       ("эндпойнт (дистиллированный)", "surrogate_adme", 0.6263, +0.0038),
       ("представление чужой модели", "minimol", 0.5990, -0.0013),
       ("дескриптор, ПЕРЕСЕКАЕТ наши", "chemprop_rdkit2d", 0.5686, -0.0089),
       ("дескриптор, не пересекает", "chemberta_77m_mtr", 0.5075, -0.0070),
       ("случайные веса той же формы", "rand2048 (1 сид)", 0.4833, float("nan"))]


def run(tag, seeds):
    emb = D + "emb_chemprop_%s.npz" % tag
    out = RES + "logs/k126_%s.json" % tag
    if not _pl.Path(emb).exists():
        print("нет %s --- сначала прогоните src/embed.py в окружении энкодера" % emb)
        return None
    for s in seeds:
        if K121.main(["--seed", str(s), "--emb", emb, "--tag", tag, "--out", out]) != 0:
            print("ОТКАЗ: k121 на сиде %d для %s" % (s, tag))
            return None
    return json.load(open(out))["seeds"]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    got = {}
    for tag, _n, _why in CANDS:
        print("\n" + "#" * 96)
        print("# %s" % tag)
        print("#" * 96, flush=True)
        r = run(tag, seeds)
        if r is None:
            return 3
        got[tag] = r

    print("\n" + "=" * 104)
    print("ДЕСКРИПТОРНЫЕ ЦЕЛИ БЕЗ ПЕРЕСЕЧЕНИЯ С НАШИМИ ПРИЗНАКАМИ: класс цели или пересечение?")
    print("=" * 104)
    print("  %-30s %-22s %9s %10s" % ("класс цели", "чекпойнт", "одиночн.", "Δ членом"))
    last = None
    for cls, nm, so, g in REF:
        print("  %-30s %-22s %9.4f %10s"
              % (cls if cls != last else "", nm, so, "%+.4f" % g if g == g else "-"))
        last = cls
    print("  " + "-" * 100)
    res = {}
    for tag, ntgt, _why in CANDS:
        r = got[tag]
        ks = sorted(r)
        so = float(np.mean([r[k]["solo_new"][0] for k in ks]))
        g = [r[k]["gain_seventh"] for k in ks]
        p = [r[k]["gain_pair_seventh"] for k in ks]
        cr = float(np.mean([float(np.mean(r[k]["errcorr_new_vs_five"])) for k in ks]))
        # per-enzyme seventh-member gain, the arm the shape block was positive on
        per = [float(np.mean([r[k]["seven"][2][e] - r[k]["six"][2][e] for k in ks]))
               for e in range(len(CYPS))]
        per_sign = [sum(1 for k in ks if r[k]["seven"][2][e] - r[k]["six"][2][e] > 0)
                    for e in range(len(CYPS))]
        print("  %-30s %-22s %9.4f %10s   знак %d/%d | Δпара %+.4f | corr %.3f | целей %d"
              % ("дескриптор, НЕ пересекает", tag, so, "%+.4f" % float(np.mean(g)),
                 sum(1 for x in g if x > 0), len(g), float(np.mean(p)), cr, ntgt))
        res[tag] = {"solo": so, "gain": float(np.mean(g)),
                    "gain_sd": float(np.std(g, ddof=1)) if len(g) > 1 else float("nan"),
                    "gain_sign": "%d/%d" % (sum(1 for x in g if x > 0), len(g)),
                    "gain_pair": float(np.mean(p)), "corr_vs_five": cr,
                    "per_enzyme": per, "per_enzyme_sign": per_sign, "n_targets": ntgt}

    print("\n  ПО ФЕРМЕНТАМ, вклад седьмым членом (блок формы пункта 279 для сравнения):")
    print("    %-14s %s" % ("", "  ".join("%9s" % c for c in CYPS)))
    print("    %-14s %s" % ("shape3d (16 ст.)", "  ".join("%+9.4f" % v for v in SHAPE3D)))
    for tag, _n, _w in CANDS:
        v, sg = res[tag]["per_enzyme"], res[tag]["per_enzyme_sign"]
        print("    %-14s %s" % (tag, "  ".join("%+9.4f" % x for x in v)))
        print("    %-14s %s" % ("", "  ".join("%9s" % ("%d/%d" % (s, len(seeds))) for s in sg)))

    print("\n  ПРОВЕРКА ПРЕДСКАЗАНИЯ (записано в докстроке до прогона: все три проваливаются):")
    refuted = []
    for tag, _n, _w in CANDS:
        v = res[tag]
        failed = v["solo"] < THRESHOLD and v["gain"] <= 0
        if not failed:
            refuted.append(tag)
        print("    %-10s одиночн. %.4f против %.2f | Δ %+.4f знак %s -> %s"
              % (tag, v["solo"], THRESHOLD, v["gain"], v["gain_sign"],
                 "провалился, как предсказано" if failed else "НЕ ПРОВАЛИЛСЯ"))
    print("\n  ЧТО ЭТО ЗНАЧИТ:")
    if refuted:
        print("    ОПРОВЕРГНУТО пунктом %s: разделяет ПЕРЕСЕЧЕНИЕ с признаками, а не класс цели."
              % ", ".join(refuted))
    else:
        print("    Пункт 349 держится: дескрипторная цель проваливается и БЕЗ пересечения с")
        print("    нашими признаками, на целях шириной 60, 1024 и 2048. Разделяет класс цели.")
    d2d6 = [(t, res[t]["per_enzyme"][2], res[t]["per_enzyme_sign"][2]) for t, _n, _w in CANDS]
    pos = [t for t, v, s in d2d6 if v > 0 and s == len(seeds)]
    print("    CYP2D6: %s" % ("3D переносится на 2D6 у %s --- поферментная находка, как блок формы"
                              % ", ".join(pos) if pos else
                              "ни у одного кандидата знак не держится 4/4 --- 3D не переносится"))
    print("\n  пол ранга %.4f, пол ST-RAE 0.007, порог пункта 348 около %.2f" % (FLOOR, THRESHOLD))

    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump({"threshold": THRESHOLD, "shape3d_per_enzyme": SHAPE3D, "reference": REF,
               "result": res, "refuted_349": refuted},
              open(RES + "logs/k126_verdict.json", "w"), ensure_ascii=False, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
