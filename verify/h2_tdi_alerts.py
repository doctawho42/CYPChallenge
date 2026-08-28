"""Несёт ли реакционная способность сигнал по TDI, до всякой квантовой химии.

Классические структурные алерты механизм-зависимой инактивации — мотивы, которые фермент
окисляет в реакционноспособный продукт. Если реакционная ось вообще работает, она должна
быть видна уже на SMARTS, без BDE и без Срнеца.

Результат отрицательный: против самой МЕТКИ TDI алерты не работают (на 3A4 максимальный
MCC около -0.013). Продолжение в h3_alerts_delta.py, где то же самое проверяется против
сдвига Delta, а не против метки, — и там сигнал есть."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D
import numpy as np, pandas as pd
from rdkit import Chem, RDLogger
RDLogger.DisableLog("rdApp.*")

ALERTS = {
    "терминальный алкин":        "[CX2;H1]#[CX2]",
    "метилендиоксифенил":        "c1cc2OCOc2cc1",
    "фуран":                     "c1ccoc1",
    "тиофен":                    "c1ccsc1",
    "первичный ароматич. амин":  "[NX3;H2][c]",
    "гидразин/гидразид":         "[NX3][NX3]",
    "тиомочевина/тиоамид":       "[NX3][CX3]=[SX1]",
    "циклопропиламин":           "[NX3][CH1]1[CH2][CH2]1",
    "терминальный алкен":        "[CX3;H2]=[CX3;H1]",
    "катехол/гидрохинон":        "c1cc(O)c(O)cc1",
    "аминофенол":                "c1cc(O)c([NX3;H2,H1])cc1",
    "тиазолидиндион":            "O=C1NC(=O)CS1",
    "N-метилпиперазин":          "[CH3][NX3]1[CH2][CH2][NX3][CH2][CH2]1",
    "бензил. C-H у гетероцикла": "[cX3][CH2][NX3]",
}
PAT = {k: Chem.MolFromSmarts(v) for k, v in ALERTS.items()}
bad = [k for k, v in PAT.items() if v is None]
if bad:
    print("не скомпилировались:", bad)

tdi = pd.read_csv(D + "cyp-challenge-TRAIN_TDI.csv")
mols = [Chem.MolFromSmiles(s) for s in tdi.SMILES]
ok = np.array([m is not None for m in mols])
H = {k: np.array([bool(m and PAT[k] and m.HasSubstructMatch(PAT[k])) for m in mols])
     for k in PAT if PAT[k]}
ANY = np.zeros(len(mols), bool)
for v in H.values():
    ANY |= v


def mcc(y, p):
    tp = int((p & y).sum()); tn = int((~p & ~y).sum())
    fp = int((p & ~y).sum()); fn = int((~p & y).sum())
    d = np.sqrt(float(tp + fp) * (tp + fn) * (tn + fp) * (tn + fn))
    return 0.0 if d == 0 else (tp * tn - fp * fn) / d


for enz in ["CYP3A4", "CYP2D6"]:
    col = f"{enz}_is_TDI"; m = tdi[col].notna().to_numpy() & ok
    y = tdi.loc[m, col].astype(bool).to_numpy()
    print(f"\n=== {enz}: размечено {m.sum()}, положительных {y.mean():.3f} ===")
    print(f"{'алерт':30s} {'встреч.':>8s} {'P(TDI|есть)':>12s} {'P(TDI|нет)':>11s} "
          f"{'лифт':>6s} {'MCC':>7s}")
    out = []
    for k in H:
        h = H[k][m]
        if h.sum() < 15:
            continue
        p1 = y[h].mean(); p0 = y[~h].mean()
        out.append((k, int(h.sum()), p1, p0, p1 / max(p0, 1e-9), mcc(y, h)))
    for k, n, p1, p0, lift, mc in sorted(out, key=lambda r: -abs(r[5])):
        print(f"{k:30s} {n:8d} {p1:12.3f} {p0:11.3f} {lift:6.2f} {mc:+7.3f}")
    a = ANY[m]
    print(f"{'ЛЮБОЙ из алертов':30s} {a.sum():8d} {y[a].mean():12.3f} {y[~a].mean():11.3f} "
          f"{y[a].mean()/max(y[~a].mean(),1e-9):6.2f} {mcc(y,a):+7.3f}")
