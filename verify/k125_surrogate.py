"""Does pretraining on PREDICTED ADME endpoints behave like measured activity, or like computed
descriptors? Two checkpoints from the organisers' own repository, through item 345's device.

WHERE THIS COMES FROM. Item 345 bought a criterion with a negative: a frozen checkpoint earns a
member slot when its pretraining task is ACTIVITY-like and is not already present in our features.
Two checkpoints work (`chemprop_medium` +0.0138, `chemprop_chemeleon` +0.0033, both pretrained on
CYP/ADME activity) and four fail (`rdkit2d` -0.0089 and chemberta-MTR -0.0070, both pretrained on
computed descriptors; two masked language models). Item 348 then measured that a candidate's own
standalone macro rank separates the two groups with no exception, at a threshold near 0.60.

A CORRECTION FIRST, because it is what made this experiment invisible. Item 345 and `src/embed.py`
both say the organisers published "sixteen pretraining checkpoints, one per descriptor family", and
item 345 concluded "the sixteen are not a stock of candidates". The repository
(`openadmet/chemprop-foundation-pretraining-weights`) holds **fifteen** checkpoint files, and two of
them are NOT descriptor families:

    surrogate_adme   25 целей: псевдодескрипторы из surrogate-данных Novartis
    minimol         512 целей: выученные вложения из MPNN-основы MiniMol

The organisers' own write-up ("throwing everything and the kitchen sink at chemeleon") says the
surrogate values were "predicted by a model trained on their highly consistent, internal corpus" --
so this is not measured assay data either. It is a DISTILLATION: a D-MPNN trained to reproduce
another model's ADME predictions. `minimol` is the same shape one level further out, a D-MPNN
trained to reproduce a representation that was itself trained on quantum and biological tasks.

WHY THAT IS THE INFORMATIVE EXPERIMENT. The two working checkpoints differ from the four failures in
two ways at once -- their task is ADME/CYP-shaped AND their labels are measured -- and nothing so
far separates those. These two candidates hold the shape and drop the measurement.

    если surrogate_adme РАБОТАЕТ -> действующее начало это ФОРМА задачи, и остальные
                                    тринадцать чекпойнтов снова становятся запасом кандидатов
    если surrogate_adme ПРОВАЛЕН  -> действующее начало это ИЗМЕРЕННЫЕ метки, и дистилляция
                                    чужой модели ничего не добавляет, как и дистилляция дескрипторов

PREDICTION, fixed here before the measurement exists, as item 347 did. **Both fail**: standalone
macro rank below 0.60 and a seventh-member contribution that is negative at sign 0/4. The reasoning
is that every distilled target measured in this project has failed -- `rdkit2d` reproduces computed
descriptors and chemberta-MTR reproduces computed descriptors, at -0.0089 and -0.0070 -- and 25
targets is a far narrower bottleneck than the 200 of `rdkit2d`. The prediction is held weakly for
`surrogate_adme`, because the ADME shape is the one property that distinguishes the two checkpoints
that work, and not at all for `minimol`, where the criterion is silent and 512 targets is wide.

The architecture is the same as `rdkit2d`'s -- d_h 2048, depth 6, LeakyReLU, 34.86 MB -- so item
343's matched random-weight control transfers directly and costs nothing to re-run: `rand2048`
scores 0.4833 standalone. Anything near that is the architecture, not the weights.

NO NEW DEVICE. This file only drives `verify/k121_chemeleon.py --emb/--tag/--out`, which is the same
code path items 345 and 348 used, including the alignment control these nameless embeddings need.

    uv run python verify/k125_surrogate.py --seeds 0,1,2,3

Reads data/emb_chemprop_{surrogate_adme,minimol}.npz, written by `src/embed.py` in the encoder
environment. Writes results/logs/k125_{tag}.json. Exit: 0 ran, 2 a control failed, 3 input missing.
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

import k121_chemeleon as K121

FLOOR = 0.0036
THRESHOLD = 0.60          # item 348, read off thirteen candidates rather than pre-registered
CANDS = [
    ("surrogate_adme", "25 целей, предсказанные ADME-эндпойнты -> ПРЕДСКАЗАН ПРОВАЛ (слабо)"),
    ("minimol", "512 целей, вложения био-обученной модели -> критерий МОЛЧИТ"),
]
# Every figure here is a four-seed mean read out of a log by verify/README.md item 348's script,
# except rand2048, which exists at seed 0 only (item 343 ran one seed).
REF = [
    ("измеренная активность", "chemprop_medium", 0.6261, +0.0138),
    ("измеренная активность", "chemprop_chemeleon", 0.6153, +0.0033),
    ("вычисленные дескрипторы", "chemprop_rdkit2d", 0.5686, -0.0089),
    ("вычисленные дескрипторы", "chemberta_77m_mtr", 0.5075, -0.0070),
    ("случайные веса, та же архитектура", "rand2048 (1 сид)", 0.4833, float("nan")),
]


def run(tag, seeds):
    """Drive k121 for one candidate, one seed at a time, and return its log."""
    emb = D + "emb_chemprop_%s.npz" % tag
    out = RES + "logs/k125_%s.json" % tag
    if not _pl.Path(emb).exists():
        print("нет %s --- сначала прогоните src/embed.py в окружении энкодера" % emb)
        return None
    for s in seeds:
        rc = K121.main(["--seed", str(s), "--emb", emb, "--tag", tag, "--out", out])
        if rc != 0:
            print("ОТКАЗ: k121 на сиде %d вернул %d для %s" % (s, rc, tag))
            return None
    return json.load(open(out))["seeds"]


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    got = {}
    for tag, _pred in CANDS:
        print("\n" + "#" * 96)
        print("# %s" % tag)
        print("#" * 96, flush=True)
        r = run(tag, seeds)
        if r is None:
            return 3
        got[tag] = r

    print("\n" + "=" * 100)
    print("ПРЕДОБУЧЕНИЕ НА ПРЕДСКАЗАННЫХ ЭНДПОЙНТАХ: форма задачи или измеренные метки?")
    print("=" * 100)
    print("  %-34s %-26s %9s %10s" % ("класс предобучения", "чекпойнт", "одиночн.", "Δ членом"))
    last = None
    for cls, nm, so, g in REF:
        print("  %-34s %-26s %9.4f %10s"
              % (cls if cls != last else "", nm, so, "%+.4f" % g if g == g else "-"))
        last = cls
    print("  " + "-" * 96)
    verdict = {}
    for tag, pred in CANDS:
        r = got[tag]
        ks = sorted(r)
        so = float(np.mean([r[k]["solo_new"][0] for k in ks]))
        sd = float(np.std([r[k]["solo_new"][0] for k in ks], ddof=1)) if len(ks) > 1 else float("nan")
        g = [r[k]["gain_seventh"] for k in ks]
        p = [r[k]["gain_pair_seventh"] for k in ks]
        cr = float(np.mean([float(np.mean(r[k]["errcorr_new_vs_five"])) for k in ks]))
        pos = sum(1 for x in g if x > 0)
        print("  %-34s %-26s %9.4f %10s   sd %.4f | знак %d/%d | Δпара %+.4f | corr %.3f"
              % ("дистилляция чужой модели", tag, so, "%+.4f" % float(np.mean(g)),
                 sd, pos, len(g), float(np.mean(p)), cr))
        verdict[tag] = {"solo": so, "solo_sd": sd, "gain": float(np.mean(g)),
                        "gain_sign": "%d/%d" % (pos, len(g)),
                        "gain_pair": float(np.mean(p)), "corr_vs_five": cr,
                        "prediction": pred}

    print("\n  ПРОВЕРКА ПРЕДСКАЗАНИЯ (оно записано в докстроке до прогона):")
    for tag, pred in CANDS:
        v = verdict[tag]
        failed = v["solo"] < THRESHOLD and v["gain"] <= 0
        print("    %-16s предсказан провал: %s | одиночн. %.4f против порога %.2f, Δ %+.4f -> %s"
              % (tag, "да" if "ПРОВАЛ" in pred else "критерий молчал",
                 v["solo"], THRESHOLD, v["gain"],
                 "ПРОВАЛИЛСЯ, как предсказано" if failed else "НЕ ПРОВАЛИЛСЯ"))
    s_ok = verdict["surrogate_adme"]
    print("\n  ЧТО ЭТО ЗНАЧИТ для механизма критерия:")
    if s_ok["solo"] >= THRESHOLD and s_ok["gain"] > 0:
        print("    surrogate_adme РАБОТАЕТ -> действующее начало это ФОРМА задачи (ADME-эндпойнты),")
        print("    а не измеренность меток. Остальные тринадцать чекпойнтов снова запас кандидатов.")
    else:
        print("    surrogate_adme ПРОВАЛЕН -> действующее начало это ИЗМЕРЕННЫЕ метки. Дистилляция")
        print("    чужой модели ведёт себя как дистилляция дескрипторов, и форма задачи не спасает.")
    print("\n  пол ранга %.4f, пол ST-RAE 0.007; порог пункта 348 около %.2f" % (FLOOR, THRESHOLD))

    _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
    json.dump({"threshold": THRESHOLD, "reference": REF, "verdict": verdict},
              open(RES + "logs/k125_verdict.json", "w"), ensure_ascii=False, indent=1)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
