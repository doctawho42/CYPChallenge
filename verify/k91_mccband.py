"""Полоса выборки для MCC на ПОДАВАЕМОМ плече классификационного трека. Апаратура, не выигрыш.

Зачем. Классификационный трек --- треть лидерборда, и у него НЕТ фальсифицируемого предсказания
к раскрытию 25 сентября: все четыре предрегистрированные полосы (пункты 246, 253, 255, 256)
относятся к ST-RAE, а по MCC полосу не считал никто (METHOD.md прямо это говорит; ни одного
бутстрэпа ни в одном из скриптов, считающих MCC). Значит любое число MCC с лидерборда нечем будет
прочитать: попало оно в ожидаемое или нет --- вопрос без заготовленного ответа.

Почему нельзя было взять готовое. Поштучных вневыборочных вероятностей ПОДАВАЕМОЙ связки нет
нигде: `bundle.json` и `tdicalib.json` --- только сводные записи (seed/cyp/arm/mcc), `oof_tdif.json`
держит поштучные, но для других плеч (предсказания как признаки), а `tdi_probs.json` от 28 августа
--- прямой классификатор, то есть плечо, которое больше не подаётся. Поэтому связка пересчитывается
здесь теми же функциями, что в `verify/k76_bundle.py`, --- плечо гарантированно то же, что мерил
пункт 250, а не похожее на него.

Подаваемое плечо (пункты 244, 250): P(ворота) из четырёхчленного ансамбля с мёртвой зоной на
pi_TDI, калиброванного вне фолда, умножается на P(сдвиг) из классификатора на (Delta > log10 2);
произведение калибруется и режется plug-in-порогом, подобранным НА ОБУЧАЮЩИХ фолдах.

ДВЕ РАЗНЫЕ ВЕЛИЧИНЫ, и путать их нельзя --- пункт 276 стоил табло ровно этой ошибки (разность
сравнивали с шумом одиночного счёта):

  1. CI ОРГАНИЗАТОРОВ --- то, что они напечатают рядом с нашим MCC. Бутстрэп С ВОЗВРАЩЕНИЕМ на
     своём же n, 1000 повторов, сид 0 --- ровно как `evaluation/config.py:BOOTSTRAP_SAMPLES = 1000`
     и `evaluation/utils.py:BOOTSTRAP_SEED = 0`. Отвечает: насколько неопределён наш MCC на ЭТОМ
     наборе соединений.
  2. ПОЛОСА РАСКРЫТИЯ --- подвыборка БЕЗ возвращения до тестового размера. Отвечает на другой
     вопрос: насколько MCC гуляет от того, КАКИЕ соединения попали в раскрываемый набор. Это
     аналог полосы ST-RAE из `k79_bandfix` (n=750 и n=375) и именно она делает предсказание
     фальсифицируемым.

Вырожденные ресэмплы (нет ни одного положительного) не выбрасываются: `matthews_corrcoef`
возвращает на них 0.0, и организаторы это специально отметили в `config.py`, так что наша полоса
обязана вести себя так же. Их доля печатается --- при доле метки 0.07-0.22 она не ноль.

Сид по умолчанию ОДИН. Полоса --- утверждение о выборке тестовых соединений при подаваемом
разбиении, а не о разбросе по сидам; точечная оценка здесь поэтому сидовая, а не четырёхсидовое
среднее пункта 250 (3A4 0.3510, 2D6 0.1235, макро 0.2373), и расхождение с ним ожидаемо и
записано. Четыре сида стоят вчетверо дороже и полосу не уточняют.

Пишет results/preds/bundle_oof.json (поштучные вневыборочные y/score/cal/dec --- отсутствующий
артефакт) и results/logs/k91_mccband.json.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parent))
from cyppaths import D, RES

import argparse, json, time
import numpy as np
import pandas as pd
from sklearn.metrics import matthews_corrcoef

from cypsplit import butina_folds
from k75_armens import member_oof, dz_refit, TDI, GATE, L
from k76_bundle import platt_oof, oof_clf, decide      # same functions as item 250

KINDS = ["поферментно", "пул", "GP", "гребневая"]
N_BOOT = 1000          # evaluation/config.py:BOOTSTRAP_SAMPLES
BOOT_SEED = 0          # evaluation/utils.py:BOOTSTRAP_SEED
SIZES = [750, 375]     # same convention as k79_bandfix's ST-RAE band


def mcc(y, d):
    return float(matthews_corrcoef(y, d))


def organiser_ci(y, d, rng):
    """Bootstrap WITH replacement at the set's own n -- the CI the organisers print."""
    n = len(y)
    vals, degen = np.empty(N_BOOT), 0
    for b in range(N_BOOT):
        i = rng.integers(0, n, n)
        yb, db = y[i], d[i]
        if yb.sum() == 0 or yb.sum() == n:
            degen += 1
        vals[b] = mcc(yb, db)
    return vals, degen


def reveal_band(y, d, n, rng):
    """Subsample WITHOUT replacement to n -- the spread over WHICH compounds are revealed."""
    N = len(y)
    if n >= N:
        return None, 0
    vals, degen = np.empty(N_BOOT), 0
    for b in range(N_BOOT):
        i = rng.choice(N, n, replace=False)
        yb, db = y[i], d[i]
        if yb.sum() == 0 or yb.sum() == n:
            degen += 1
        vals[b] = mcc(yb, db)
    return vals, degen


def summarise(v):
    return {"mean": float(v.mean()), "sd": float(v.std(ddof=1)),
            "lo95": float(np.percentile(v, 2.5)), "hi95": float(np.percentile(v, 97.5))}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0", help="разбиения; полоса --- про выборку соединений, "
                                                "не про сиды, поэтому по умолчанию один")
    ap.add_argument("--out", default=RES + "logs/k91_mccband.json")
    a = ap.parse_args()

    z = np.load(D + "feats.npz")
    X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    rows = pd.read_csv(D + "rows.csv")
    tdi = (pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv")
             .set_index("Molecule_Name").reindex(rows.Molecule_Name))
    inh = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
             .set_index("Molecule_Name").reindex(rows.Molecule_Name))

    out, per_compound = {}, {}
    for seed in [int(s) for s in a.seeds.split(",") if s.strip()]:
        fold, _ = butina_folds(list(rows.SMILES), seed=seed)
        print(f"\n=== сид {seed}: воспроизводим подаваемую связку вне фолда ===", flush=True)
        Xs, ys, folds, LOs, HIs, labs, shifts, gates = [], [], [], [], [], [], [], []
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
            labs.append(lab[m].astype(bool).to_numpy())
            shifts.append((ta - dr)[m]); gates.append((ta > GATE)[m])
            print(f"  {c}: n {int(m.sum())}, доля метки {lab[m].astype(bool).mean():.3f}",
                  flush=True)

        dz = {}
        for k in KINDS:
            t0 = time.time()
            P = member_oof(k, Xs, ys, folds)
            tg = [np.clip(P[e], LOs[e], HIs[e]) for e in range(len(TDI))]
            dz[k] = dz_refit(k, Xs, tg, folds)
            print(f"    плечо, член {k}: {time.time() - t0:.0f} с", flush=True)
        arm_ens = [np.mean([dz[k][e] for k in KINDS], 0) for e in range(len(TDI))]

        res = {}
        for e, c in enumerate(TDI):
            y, fi = labs[e], folds[e]
            dp = oof_clf(Xs[e], shifts[e] > L, fi, f"{c} сдвиг")
            gp = platt_oof(arm_ens[e] - GATE, gates[e], fi, raw=True)
            score = gp * dp                          # СВЯЗКА, exactly k76's shipped arm
            dec, cal = decide(score, y, fi)
            point = mcc(y, dec)
            print(f"  {c}: MCC(точечно, сид {seed}) {point:+.4f}   "
                  f"доля положительных решений {dec.mean():.3f}", flush=True)

            rng = np.random.default_rng(BOOT_SEED)
            ci, dg = organiser_ci(y, dec, rng)
            rec = {"n": int(len(y)), "prevalence": float(y.mean()), "mcc": point,
                   "ci_organisers": summarise(ci), "ci_degenerate": int(dg), "reveal": {}}
            print(f"      CI организаторов (с возвр., n={len(y)}, 1000): "
                  f"[{rec['ci_organisers']['lo95']:+.4f}, {rec['ci_organisers']['hi95']:+.4f}]"
                  f"  sd {rec['ci_organisers']['sd']:.4f}  вырожденных {dg}", flush=True)
            for n in SIZES:
                v, dgn = reveal_band(y, dec, n, np.random.default_rng(BOOT_SEED + n))
                if v is None:
                    print(f"      полоса раскрытия n={n}: пропущена (n >= {len(y)})", flush=True)
                    continue
                s = summarise(v); s["degenerate"] = int(dgn)
                rec["reveal"][str(n)] = s
                print(f"      полоса раскрытия n={n} (без возвр.): "
                      f"[{s['lo95']:+.4f}, {s['hi95']:+.4f}]  sd {s['sd']:.4f}  "
                      f"полуширина {(s['hi95'] - s['lo95']) / 2:.4f}  вырожденных {dgn}",
                      flush=True)
            res[c] = rec
            per_compound[f"{seed}|{c}"] = {"y": y.astype(int).tolist(),
                                           "score": score.tolist(),
                                           "cal": cal.tolist(),
                                           "dec": dec.astype(int).tolist()}

        # Макро: MCC двух ферментов усредняется ВНУТРИ каждого ресэмпла, а не после.
        # Наборы соединений у ферментов разные, поэтому каждый ресэмплится своим потоком.
        for tag, sizes in (("ci_organisers", [None]), ("reveal", SIZES)):
            for n in sizes:
                vs = []
                for e, c in enumerate(TDI):
                    y = labs[e]
                    d = np.asarray(per_compound[f"{seed}|{c}"]["dec"], bool)
                    if tag == "ci_organisers":
                        v, _ = organiser_ci(y, d, np.random.default_rng(BOOT_SEED))
                    else:
                        v, _ = reveal_band(y, d, n, np.random.default_rng(BOOT_SEED + n))
                    if v is None:
                        vs = []
                        break
                    vs.append(v)
                if not vs:
                    continue
                m = np.mean(vs, axis=0)
                key = "макро CI организаторов" if tag == "ci_organisers" else f"макро полоса n={n}"
                res[key] = summarise(m)
                print(f"  {key}: [{res[key]['lo95']:+.4f}, {res[key]['hi95']:+.4f}]  "
                      f"sd {res[key]['sd']:.4f}  среднее {res[key]['mean']:+.4f}", flush=True)

        out[str(seed)] = res
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(a.out, "w"), ensure_ascii=False, indent=1)
        json.dump(per_compound, open(RES + "preds/bundle_oof.json", "w"))
        print(f"\n  сохранено: {a.out}\n  сохранено: {RES}preds/bundle_oof.json", flush=True)

    print("""
Как читать. Точечный MCC здесь --- ОДНОСИДОВЫЙ, поэтому он не обязан совпадать с
четырёхсидовым средним пункта 250 (3A4 0.3510, 2D6 0.1235, макро 0.2373); расхождение между
ними --- сидовый разброс, а не дефект.

CI организаторов --- то, что напечатают рядом с нашим числом. Полоса раскрытия --- то, во что
наше число имеет право попасть, и только она делает предсказание к 25 сентября проверяемым.
Складывать их нельзя: первая про неопределённость на данном наборе, вторая про то, какой набор
достанется. Если полоса окажется шире любого выигрыша, когда-либо померенного на этом треке
(калибровка +0.0133, связка +0.0127 при поле 0.0076) --- это и есть результат, а не разочарование:
он означает, что одно раскрытие не способно различить наши плечи по MCC.""")


if __name__ == "__main__":
    main()
