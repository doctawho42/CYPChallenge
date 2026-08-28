"""Чем на самом деле является «сдвиг центра усадки» и откуда берётся его выигрыш.

Три вещи, которые k1 оставил неразобранными.

1. Параметризация вводит в заблуждение. c + l*(p - c) = [mu + l*(p - mu)] + (1-l)*off,
   то есть предсказания смещаются на (1-l)*off, а не на off. При l около 0.8 «оптимум
   +0.40» из k1 — это +0.08 в шкале pIC50. Семейство (c, l) при этом вырождается в
   двухпараметрическое a + b*p с a = c*(1-l), b = l, так что сам c определён только при
   l != 1, и вблизи единицы его оптимум неустойчив.

2. Откуда выигрыш. Гипотеза «бесплатный обед на цензурированных полосах» — предсказание
   уезжает вверх внутри широкой полосы неактивного и ничего не стоит — НЕВЕРНА: доля
   предсказаний внутри полосы при сдвиге падает. Выигрыш приходит из средней зоны и с
   активных, а платят за него неактивные, и платят меньше.

3. Переносится ли (c, l), подобранные на обычных фолдах, на оценочный набор, обогащённый
   активными. Колонка «оракул» — те же параметры, подобранные на самом обогащённом наборе;
   это недостижимая сверху граница, и разрыв до неё показывает цену незнания."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd, json
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from cypsplit import butina_folds

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
rows = pd.read_csv(D + "rows.csv")
tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
        .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
oof = json.load(open(RES + "preds/oof.json"))
fold, _ = butina_folds(list(rows.SMILES))
OFF = np.round(np.arange(-1.5, 2.51, 0.05), 2)
LAM = np.round(np.arange(0.30, 1.31, 0.02), 2)


def pack(c):
    col = f"{c}_pIC50_direct_inhibition"; m = tr[col].notna().to_numpy()
    return (tr.loc[m, col].to_numpy(), tr.loc[m, col + "_conf_low"].to_numpy(),
            tr.loc[m, col + "_conf_high"].to_numpy(),
            np.asarray(oof[f"FP+DESC+MECH|{c}"]), fold[m])


def S(y, p, lo, hi):
    return strae(y, p, y_true_upper=hi, y_true_lower=lo)


def num(p, lo, hi):
    return np.maximum(p - hi, 0) + np.maximum(lo - p, 0)


print("=" * 80)
print("1. c + l*(p - c) = [mu + l*(p - mu)] + (1-l)*off")
print("   смещение предсказаний равно (1-l)*off, а не off")
print("=" * 80)
print(f"{'фермент':8s} {'off*':>6s} {'l*':>5s} {'сдвиг предск.':>13s} {'ST-RAE p':>9s} "
      f"{'c=mu':>8s} {'c=mu+off':>9s}")
FIT = {}
for c in CYPS:
    y, lo, hi, p, f = pack(c)
    pA = np.zeros_like(p); pB = np.zeros_like(p); offs = []; lamA = []; lamB = []
    for k in range(5):
        a, b = f != k, f == k
        mu = y[a].mean()
        lA = LAM[np.argmin([S(y[a], mu + L * (p[a] - mu), lo[a], hi[a]) for L in LAM])]
        G = np.array([[S(y[a], mu + o + L * (p[a] - mu - o), lo[a], hi[a]) for L in LAM]
                      for o in OFF])
        i, j = np.unravel_index(G.argmin(), G.shape)
        offs.append(OFF[i]); lamA.append(lA); lamB.append(LAM[j])
        pA[b] = mu + lA * (p[b] - mu)
        pB[b] = mu + OFF[i] + LAM[j] * (p[b] - mu - OFF[i])
    o, lB = np.mean(offs), np.mean(lamB)
    FIT[c] = (o, lB, np.mean(lamA), pA, pB)
    print(f"{c:8s} {o:+6.2f} {lB:5.2f} {(1-lB)*o:+13.3f} {S(y,p,lo,hi):9.4f} "
          f"{S(y,pA,lo,hi):8.4f} {S(y,pB,lo,hi):9.4f}")
MAC = [0.0, 0.0, 0.0]
for c in CYPS:
    y, lo, hi, p, f = pack(c)
    for i, q in enumerate([p, FIT[c][3], FIT[c][4]]):
        MAC[i] += S(y, q, lo, hi) / len(CYPS)
print(f"{'МАКРО':8s} {'':6s} {'':5s} {'':13s} {MAC[0]:9.4f} {MAC[1]:8.4f} {MAC[2]:9.4f}")

print()
print("=" * 80)
print("2. Откуда добавка от сдвига: числитель B(со сдвигом) - A(центр=mu)")
print("=" * 80)
for c in CYPS:
    y, lo, hi, p, f = pack(c); o, lB, lA, pA, pB = FIT[c]
    z = np.digitize(y, [4.5, 5.5])
    nA, nB = num(pA, lo, hi), num(pB, lo, hi)
    tot = nB.sum() - nA.sum()
    print(f"\n{c}  сдвиг предсказаний {(1-lB)*o:+.3f}   d(числителя) {tot:+.1f}")
    print(f"  {'зона':>12s} {'n':>5s} {'A':>8s} {'B':>8s} {'d':>8s} {'доля':>6s}")
    for k, nm in enumerate(["pIC50<4.5", "4.5-5.5", ">5.5"]):
        m = z == k; dd = nB[m].sum() - nA[m].sum()
        print(f"  {nm:>12s} {m.sum():5d} {nA[m].sum():8.1f} {nB[m].sum():8.1f} {dd:+8.1f} "
              f"{dd/tot*100 if tot else 0:5.0f}%")
    print(f"  внутри полосы A {((pA>=lo)&(pA<=hi)).mean()*100:.1f}% -> "
          f"B {((pB>=lo)&(pB<=hi)).mean()*100:.1f}%")

print()
print("=" * 80)
print("3. Перенос. Карта подобрана на обычных фолдах; оценка на обогащённом наборе.")
print("   «оракул» — (off, l), подобранные на самом обогащённом наборе, недостижимо")
print("=" * 80)
print(f"{'фермент':8s} {'набор':>8s} {'n':>5s} {'ybar':>6s} | {'сырое':>7s} {'c=mu':>7s} "
      f"{'c=mu+off':>8s} | {'ор. off':>7s} {'l':>5s} {'оракул':>7s}")
for c in CYPS:
    y, lo, hi, p, f = pack(c); o, lB, lA, pA, pB = FIT[c]
    mu = y.mean()
    for nm, q in [("случ.", 0.0), ("top50%", 0.50), ("top25%", 0.75), ("top10%", 0.90)]:
        m = y >= (np.quantile(y, q) if q > 0 else -np.inf)
        yy, ll, hh = y[m], lo[m], hi[m]
        G = np.array([[S(yy, mu + oo + L * (p[m] - mu - oo), ll, hh) for L in LAM]
                      for oo in OFF])
        i, j = np.unravel_index(G.argmin(), G.shape)
        print(f"{c if q == 0 else '':8s} {nm:>8s} {m.sum():5d} {yy.mean():6.3f} | "
              f"{S(yy,p[m],ll,hh):7.4f} {S(yy,pA[m],ll,hh):7.4f} {S(yy,pB[m],ll,hh):8.4f} | "
              f"{OFF[i]:+7.2f} {LAM[j]:5.2f} {G[i,j]:7.4f}")
