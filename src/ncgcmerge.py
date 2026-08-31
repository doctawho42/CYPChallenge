"""Fold the six NCGC assays into one panel table, and measure the overlap before using it.

What this is for. Items 60 to 63 measured why the external ChEMBL set failed: it is a sample
selected by *publishability*, its labels sit 0.22 to 0.60 above ours, and no shift correction
repairs that because the distortion is selection rather than bias. The NCGC panel is the opposite
construction -- one laboratory, one protocol, one campaign, the whole library run at sixteen
concentrations -- so there is no publication filter to correct for. `src/fetch1851.py` downloaded
it; nothing has been merged.

What is here, as fetched:

    AID   изоформа      строк    с AC50    схема
    410   CYP1A2         9174      5169    'Log of AC50', 16 колонок с QC и моделью кривой
    883   CYP2C9        10296      3068    'Fit_LogAC50', 11 колонок
    891   CYP2D6        10296      3637    'Fit_LogAC50'
    899   CYP2C19       10296      3698    'Fit_LogAC50'  --- пятая изоформа, которой у нас нет
    884   CYP3A4        14115      6837    'Fit_LogAC50'
    411   люцифераза    72335      2262    контрскрин

The two schemas are not a nuisance to paper over: AID 410 alone carries `Curve Fit Model` and the
QC columns, so a quality filter that exists for CYP1A2 does not exist for the other four. That is
stated in the output rather than silently equalised.

Units. Both spellings are log10(AC50 in molar), so pAC50 = -x. The medians land at 5.30 (CYP1A2)
and 5.10 (CYP3A4), against our own labels' 5.13 -- the same scale, which is the first thing that
had to be true and was not true of ChEMBL.

Inactives are kept as a censored row, not dropped. A compound with no fitted AC50 in a qHTS run
is a measurement -- "not active up to the top concentration" -- and throwing it away would rebuild
exactly the selection bias that killed the ChEMBL merge. They are marked, and it is the ablation's
choice whether to use them.

**The overlap with the blinded test set is measured here and reported, and nothing is decided.**
If NCGC contains compounds that are in `cyp-challenge-TEST-BLINDED.csv` with a fitted AC50 on the
same enzyme, then training on the panel hands the model measured answers for test rows. That is
public prior data rather than a leak in the ordinary sense, but whether it may be used is a
question about the competition's rules and not one this script can settle. The number is printed
first, before anything else, so it cannot be skipped.

Matching is on the RDKit canonical SMILES of the parent structure. NCGC ships
`PUBCHEM_EXT_DATASOURCE_SMILES` with the assay records, so no compound lookup is needed.

Writes data/ncgc/panel.csv. ~4 minutes, most of it canonicalisation.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D

import argparse

import numpy as np
import pandas as pd
from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")

ASSAYS = {410: "CYP1A2", 883: "CYP2C9", 891: "CYP2D6", 899: "CYP2C19", 884: "CYP3A4"}
LUC = 411
CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]


def canon(smiles):
    """Каноническая форма без соли-противоиона: берётся крупнейший фрагмент."""
    m = Chem.MolFromSmiles(smiles) if isinstance(smiles, str) else None
    if m is None:
        return None
    frags = Chem.GetMolFrags(m, asMols=True, sanitizeFrags=False)
    if len(frags) > 1:
        m = max(frags, key=lambda f: f.GetNumHeavyAtoms())
    try:
        return Chem.MolToSmiles(m)
    except Exception:
        return None


def load(aid):
    d = pd.read_csv(D + f"ncgc/aid_{aid}.csv", low_memory=False)
    col = "Log of AC50" if "Log of AC50" in d.columns else "Fit_LogAC50"
    out = pd.DataFrame({
        "cid": pd.to_numeric(d.get("PUBCHEM_CID"), errors="coerce"),
        "smiles_raw": d.get("PUBCHEM_EXT_DATASOURCE_SMILES"),
        "outcome": d.get("PUBCHEM_ACTIVITY_OUTCOME"),
        "score": pd.to_numeric(d.get("PUBCHEM_ACTIVITY_SCORE"), errors="coerce"),
        "logac50": pd.to_numeric(d[col], errors="coerce"),
    })
    for src, dst in (("Hill Coefficient", "hill"), ("Fit_HillSlope", "hill"),
                     ("Curve R2", "r2"), ("Fit_R2", "r2")):
        if src in d.columns and dst not in out.columns:
            out[dst] = pd.to_numeric(d[src], errors="coerce")
    for src in ("Curve Fit Model", "Fit_CurveClass", "Curve_Description"):
        if src in d.columns:
            out["curve"] = d[src].astype(str)
            break
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=D + "ncgc/panel.csv")
    a = ap.parse_args()

    raw = {aid: load(aid) for aid in list(ASSAYS) + [LUC]}

    # Канонизируется каждый уникальный SMILES один раз, а не каждая строка.
    uniq = pd.unique(pd.concat([r.smiles_raw for r in raw.values()]).dropna())
    print(f"уникальных SMILES в панели: {len(uniq)}, канонизирую...", flush=True)
    cmap = {s: canon(s) for s in uniq}
    bad = sum(v is None for v in cmap.values())
    print(f"  RDKit не разобрал: {bad}\n")

    luc = raw[LUC].assign(smiles=raw[LUC].smiles_raw.map(cmap))
    luc_active = set(luc.loc[luc.outcome == "Active", "smiles"].dropna())

    parts = []
    for aid, enz in ASSAYS.items():
        d = raw[aid].assign(smiles=raw[aid].smiles_raw.map(cmap), enzyme=enz)
        d = d[d.smiles.notna()].copy()
        d["pAC50"] = -d.logac50
        d["censored"] = d.pAC50.isna() & (d.outcome == "Inactive")
        d["luc"] = d.smiles.isin(luc_active)
        parts.append(d.drop(columns=["smiles_raw", "logac50"]))
    panel = pd.concat(parts, ignore_index=True)

    # --- ПЕРЕСЕЧЕНИЯ. Печатаются первыми и до всего остального. ---
    rows = pd.read_csv(D + "rows.csv")
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    tr_c = {canon(s) for s in rows.SMILES}
    te_c = {canon(s) for s in te.SMILES}
    tr_c.discard(None); te_c.discard(None)

    print("=" * 74)
    print("ПЕРЕСЕЧЕНИЕ ПАНЕЛИ С НАШИМИ НАБОРАМИ (канонический SMILES, крупнейший фрагмент)")
    print("=" * 74)
    print(f"{'фермент':9s} {'строк':>7s} {'с AC50':>7s} {'∩ обучение':>11s} "
          f"{'∩ ТЕСТ':>7s} {'∩ТЕСТ с AC50':>13s}")
    tot_te = 0
    for enz in ASSAYS.values():
        d = panel[panel.enzyme == enz]
        it = d[d.smiles.isin(te_c)]
        tot_te += int(it.pAC50.notna().sum())
        print(f"{enz:9s} {len(d):7d} {int(d.pAC50.notna().sum()):7d} "
              f"{int(d.smiles.isin(tr_c).sum()):11d} {len(it):7d} "
              f"{int(it.pAC50.notna().sum()):13d}")
    print(f"\nмолекул теста в панели: {len({s for s in panel.smiles if s in te_c})} из {len(te_c)}")
    print(f"измеренных AC50 на тестовых молекулах по нашим четырём ферментам: {tot_te}")
    if tot_te:
        print("""
!!! ЭТО РЕШЕНИЕ НЕ СКРИПТА. Панель содержит измеренные AC50 для молекул ослеплённого
    теста. Это открытые данные, а не утечка в обычном смысле, но обучение на них даёт
    модели ответы на тестовые строки, и допустимо ли это --- вопрос правил соревнования.
    Число напечатано; решать человеку. Пока оно не решено, ablncgc.py считает руку
    «без тестовых молекул» контролем и сравнивает с ней.""")
    print("=" * 74 + "\n")

    print(f"{'фермент':9s} {'активных':>9s} {'цензурир.':>10s} {'люцифераза+':>12s} "
          f"{'медиана pAC50':>14s}")
    for enz in ASSAYS.values():
        d = panel[panel.enzyme == enz]
        print(f"{enz:9s} {int(d.pAC50.notna().sum()):9d} {int(d.censored.sum()):10d} "
              f"{int(d.luc.sum()):12d} {d.pAC50.median():14.2f}")
    print(f"\nQC-колонки есть только у AID 410: {'curve' in raw[410].columns}. "
          f"У остальных фильтра по модели кривой нет, и это не выравнивается.")

    panel.to_csv(a.out, index=False)
    print(f"\nсохранено: {a.out}  ({len(panel)} строк, "
          f"{panel.smiles.nunique()} уникальных структур)")


if __name__ == "__main__":
    main()
