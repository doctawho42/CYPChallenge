"""Метка TDI — это «достаточно потентен» И «сдвинулся». Алерты про второе.

h2 проверял алерты против самой метки и получил ноль. Здесь то же самое проверяется
против сдвига Delta = pIC50(TDI) - pIC50(прямое), обусловленного на потентность (pi > 4),
то есть против той половины определения, о которой алерты вообще что-то говорят.

Там сигнал есть и он сильный. Заодно видно, почему метка сама по себе шумная: медиана
Delta на активных 3A4 равна +0.294 против порога log10(2) = 0.301 — отсечка проходит
ровно по середине распределения."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
from scipy.stats import mannwhitneyu
RDLogger.DisableLog("rdApp.*")

ALERTS = {"циклопропиламин": "[NX3][CH1]1[CH2][CH2]1", "метилендиоксифенил": "c1cc2OCOc2cc1",
          "фуран": "c1ccoc1", "тиофен": "c1ccsc1", "катехол/гидрохинон": "c1cc(O)c(O)cc1",
          "первичный ароматич. амин": "[NX3;H2][c]", "терминальный алкин": "[CX2;H1]#[CX2]",
          "терминальный алкен": "[CX3;H2]=[CX3;H1]", "тиазолидиндион": "O=C1NC(=O)CS1"}
PAT = {k: Chem.MolFromSmarts(v) for k, v in ALERTS.items()}
t = pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv")
mols = [Chem.MolFromSmiles(s) for s in t.SMILES]
H = {k: np.array([bool(m and m.HasSubstructMatch(PAT[k])) for m in mols]) for k in PAT}
ANY = np.zeros(len(mols), bool)
for v in H.values():
    ANY |= v
L2 = np.log10(2)

for enz in ["CYP3A4", "CYP2D6"]:
    a = t[f"{enz}_pIC50_direct_inhibition"].to_numpy()
    b = t[f"{enz}_pIC50_TDI_condition"].to_numpy()
    m = ~np.isnan(a) & ~np.isnan(b)
    d = b - a
    act = m & (a > 4)
    print(f"\n=== {enz}: обе величины у {m.sum()}, из них активных (pi>4) {act.sum()} ===")
    print(f"    медиана Delta: все {np.median(d[m]):+.3f} | активные {np.median(d[act]):+.3f}"
          f" | порог log10(2) = {L2:.3f}")
    print(f"{'алерт':28s} {'n(акт)':>7s} {'медD|есть':>10s} {'медD|нет':>9s} "
          f"{'доля D>log2':>12s} {'p':>8s}")
    for k in PAT:
        h = H[k] & act
        if h.sum() < 12:
            continue
        o = act & ~H[k]
        p = mannwhitneyu(d[h], d[o], alternative="greater").pvalue
        print(f"{k:28s} {h.sum():7d} {np.median(d[h]):10.3f} {np.median(d[o]):9.3f} "
              f"{(d[h]>L2).mean():12.3f} {p:8.4f}")
    h = ANY & act; o = act & ~ANY
    p = mannwhitneyu(d[h], d[o], alternative="greater").pvalue
    print(f"{'ЛЮБОЙ':28s} {h.sum():7d} {np.median(d[h]):10.3f} {np.median(d[o]):9.3f} "
          f"{(d[h]>L2).mean():12.3f} {p:8.4f}   (без алерта доля {(d[o]>L2).mean():.3f})")
