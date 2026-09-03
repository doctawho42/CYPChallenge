"""Oracle C: how much of the ensemble's pairwise error lives on pairs its members argue about.

Over 200 items this file has measured ensemble spread per COMPOUND and never per PAIR: how
often the five members disagree about which of two compounds is the stronger inhibitor. That is
a different object, and it is the one the metric pays for after the affine pair, since the pair
is strictly increasing and preserves every pairwise comparison exactly.

**Why Kendall and not Spearman.** The proposal is "impose the true order on contested pairs and
keep the model's elsewhere". As a construction on Spearman that is ill-posed: fixing an
arbitrary subset of pairwise comparisons need not yield a consistent total order, and cycles
are not a hypothetical when members disagree. Kendall's tau has no such problem because it IS a
sum over pairs,

    tau = (concordant - discordant) / C(n, 2),

so "repair every contested discordant pair" is exactly delta_tau = 2 * D_contested / C(n,2),
with no ordering to construct. Kendall and Spearman are monotonically related and the affine
pair preserves both, so nothing about the criterion is given up by measuring the one that
decomposes.

What decides whether there is a target here is not the ceiling but the LIFT: contested pairs are
some fraction of all pairs and carry some fraction of the discordance. If those two fractions
are equal, disagreement among members carries no information about where the ensemble is wrong,
a second cascade has nothing to aim at, and the direction closes for an hour's work. If
discordance concentrates on contested pairs, the target exists and is addressable without
touching the pairs that are already settled -- which is the rare kind of intervention that
cannot lose much.

Three further questions the same pass answers, because they cost nothing once the sign matrices
exist and they decide what a cascade could even do:

    on contested pairs, is the ensemble MEAN better or worse than a coin;
    is any single member systematically better there than the mean;
    does majority vote among the five beat the mean there.

If the mean is already the best available combination on contested pairs, then the disagreement
is noise rather than signal and there is nothing to route.

Reads the member cache written by verify/k58_dzsubmit.py --cache. Minutes.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()

import argparse
import json

import numpy as np
import pandas as pd

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--cache", required=True, help="кэш членов из k58_dzsubmit.py")
    ap.add_argument("--which", default="dzp", choices=["parts", "dzp"],
                    help="parts: члены без прохода; dzp: с мёртвой зоной (подаётся)")
    a = ap.parse_args()

    rows = pd.read_csv(D + "rows.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    C = json.load(open(a.cache))
    members = [(k, [np.asarray(v, float) for v in P]) for k, P in C[a.which]]
    print(f"кэш {a.cache}, набор {a.which}, сид {C.get('seed')}, "
          f"членов {len(members)}: " + ", ".join(k for k, _ in members) + "\n", flush=True)

    T = []
    for e, c in enumerate(CYPS):
        col = f"{c}_pIC50_direct_inhibition"
        m = tr[col].notna().to_numpy()
        y = tr.loc[m, col].to_numpy()
        Pm = np.stack([P[e] for _, P in members])          # (члены, n)
        ens = Pm.mean(0)
        n = len(y)
        iu = np.triu_indices(n, 1)

        ty = np.sign(y[iu[0]] - y[iu[1]])
        keep = ty != 0                                     # связи в истине не считаем
        ty = ty[keep]
        te_ = np.sign(ens[iu[0]] - ens[iu[1]])[keep]
        sm = np.stack([np.sign(P[e][iu[0]] - P[e][iu[1]])[keep] for _, P in members])

        contested = ~np.all(sm == sm[0], axis=0)
        disc = te_ != ty
        npair = len(ty)

        f_cont = contested.mean()
        f_disc = disc.mean()
        f_disc_given_cont = disc[contested].mean() if contested.any() else np.nan
        share = (disc & contested).sum() / max(disc.sum(), 1)
        lift = share / max(f_cont, 1e-12)

        tau = (~disc).mean() - disc.mean()
        tau_fixed = tau + 2.0 * (disc & contested).sum() / npair

        # что можно было бы сделать НА спорных парах
        acc_mean = (~disc)[contested].mean() if contested.any() else np.nan
        maj = np.sign(sm[:, contested].sum(0))
        acc_maj = (maj == ty[contested]).mean() if contested.any() else np.nan
        acc_best = max(((sm[k, contested] == ty[contested]).mean()
                        for k in range(len(members))), default=np.nan)

        T.append({"фермент": c, "n": n, "пар": npair,
                  "спорных": f_cont, "несогл.": f_disc,
                  "несогл|спорн": f_disc_given_cont,
                  "доля несогл. на спорных": share, "подъём": lift,
                  "tau": tau, "tau почин.": tau_fixed, "потолок": tau_fixed - tau,
                  "точн. среднего": acc_mean, "точн. больш-ва": acc_maj,
                  "лучший член": acc_best})

    df = pd.DataFrame(T).set_index("фермент")
    print("СКОЛЬКО СПОРНОГО И ГДЕ ЖИВЁТ ОШИБКА")
    print(df[["n", "пар", "спорных", "несогл.", "несогл|спорн",
              "доля несогл. на спорных", "подъём"]].round(4).to_string())
    print("\nПОТОЛОК ВТОРОГО КАСКАДА (Кендалл)")
    print(df[["tau", "tau почин.", "потолок"]].round(4).to_string())
    print("\nЕСТЬ ЛИ ЧТО МАРШРУТИЗИРОВАТЬ НА СПОРНЫХ ПАРАХ")
    print(df[["точн. среднего", "точн. больш-ва", "лучший член"]].round(4).to_string())
    print(f"\nмакро: спорных {df['спорных'].mean():.3f}, подъём {df['подъём'].mean():.3f}, "
          f"потолок Кендалла {df['потолок'].mean():+.4f}")
    print("""
Как читать. ПОДЪЁМ --- главное число. Это доля несогласованных пар, лежащих на спорных,
делённая на долю спорных среди всех. Единица означает, что расхождение членов не знает
ничего о том, где ансамбль ошибается: спорные пары ошибочны ровно с той же частотой, что и
любые, мишени нет, направление закрыто. Заметно больше единицы --- мишень есть.

ПОТОЛОК --- прирост Кендалла, если бы все спорные несогласованные пары стали согласованными.
Он недостижим по построению (в него подставлена истина), и опыт пунктов 126 и 128 в том, что
реальное составляет от потолка малую долю.

Последняя таблица говорит, чем вообще можно было бы починить. Если точность среднего на
спорных парах не ниже точности большинства и лучшего одиночного члена, то на этих парах
среднее уже лучшее из доступного, расхождение --- шум, а не сигнал, и маршрутизировать
нечего, каким бы ни был потолок.""")


if __name__ == "__main__":
    main()
