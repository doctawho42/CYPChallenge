"""Fetch the NCGC cytochrome panel from PubChem, with its own luciferase counter-screen.

Why this and not more ChEMBL. Items 60 to 63 measured what went wrong with the external set the
organisers supplied: ChEMBL is a sample selected by *publishability*, its labels sit 0.22 to 0.60
above ours, and no shift correction repairs that because the distortion is selection rather than
bias. The NCGC panel is the opposite construction -- one laboratory, one protocol, one campaign,
the whole library screened, so there is no publication filter to correct for.

What it supplies that our training set structurally lacks:

    AID   изоформа                    соединений      у нас
    410   CYP1A2                            8354       1412
    883   CYP2C9                            9385       1285
    891   CYP2D6                            9385       1493
    899   CYP2C19                           9385          -
    884   CYP3A4                           13076       2335
    411   контрскрин на люциферазу         71298          -

Four to six times the compounds on the same enzymes, a full dose-response with sixteen
concentrations and a fitted AC50 rather than a single point, an unselected population, and CYP2C19
as a fifth isoform -- without which a proteochemometric enzyme coordinate has nothing to be fitted
on. Our own matrix carries all four enzymes for exactly 41 compounds; this one carries the panel
for the whole library.

The luciferase caveat, and why it is a measurement here rather than a caution. This is a
luciferin-coupled readout, and inhibitors of firefly luciferase masquerade as inhibitors of CYP in
exactly this format. The usual answer is care. The better answer is that NCGC ran the counter-
screen on the same library and published it as AID 411, covering 8422 of the 13076 compounds in
the CYP3A4 assay. So interference stops being a hypothesis to reason around and becomes a column:
excluded, or carried as a covariate, and either choice is checkable.

What is fetched per assay: the fitted `Fit_LogAC50`, the curve class and description, the activity
outcome and score. `Fit_CurveClass` matters as much as the potency -- NCGC's classes distinguish a
complete sigmoid from a partial or single-point response, and pooling those without the
distinction is how a qHTS set poisons a model.

Nothing is merged here. This script only downloads and caches; the join, the source indicator and
the interference handling belong to a separate ablation, because item 63's lesson is that the
merge is where external data goes wrong.

Writes data/ncgc/aid_<id>.csv and data/ncgc/cid_smiles.csv.
"""
import argparse
import pathlib
import time
import urllib.request

import pandas as pd

ROOT = pathlib.Path(__file__).resolve().parents[1]
D = ROOT / "data" / "ncgc"
BASE = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
ASSAYS = {410: "CYP1A2", 883: "CYP2C9", 891: "CYP2D6", 899: "CYP2C19",
          884: "CYP3A4", 411: "luciferase"}
KEEP = ["PUBCHEM_CID", "PUBCHEM_ACTIVITY_OUTCOME", "PUBCHEM_ACTIVITY_SCORE",
        "Fit_LogAC50", "Fit_HillSlope", "Fit_R2", "Fit_CurveClass",
        "Curve_Description", "Fit_InfiniteActivity", "Potency"]


def get(url, tries=4, pause=3.0):
    """PubChem отдаёт 503 при нагрузке, поэтому повтор с паузой обязателен."""
    for k in range(tries):
        try:
            with urllib.request.urlopen(url, timeout=180) as r:
                return r.read().decode("utf-8", "replace")
        except Exception as e:
            if k == tries - 1:
                raise
            print(f"    повтор {k+1}: {type(e).__name__}", flush=True)
            time.sleep(pause * (k + 1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-smiles", action="store_true")
    a = ap.parse_args()
    D.mkdir(parents=True, exist_ok=True)

    cids = set()
    for aid, name in ASSAYS.items():
        out = D / f"aid_{aid}.csv"
        if out.exists():
            df = pd.read_csv(out)
            print(f"AID {aid} ({name}): уже есть, {len(df)} строк")
        else:
            t0 = time.time()
            txt = get(f"{BASE}/assay/aid/{aid}/CSV")
            import io
            df = pd.read_csv(io.StringIO(txt), low_memory=False)
            # Первые строки CSV PubChem --- служебные описания типов, у них нет CID.
            df = df[pd.to_numeric(df.get("PUBCHEM_CID"), errors="coerce").notna()]
            cols = [c for c in KEEP if c in df.columns]
            df = df[cols]
            df.to_csv(out, index=False)
            print(f"AID {aid} ({name}): {len(df)} строк, {len(cols)} колонок, "
                  f"{time.time()-t0:.0f} с", flush=True)
        if aid != 411:
            cids |= set(pd.to_numeric(df.PUBCHEM_CID, errors="coerce").dropna().astype(int))

    print(f"\nуникальных соединений по пяти изоформам: {len(cids)}")
    sm = D / "cid_smiles.csv"
    if a.skip_smiles or sm.exists():
        print("SMILES пропущены" if a.skip_smiles else f"SMILES уже есть: {sm}")
        return
    ids = sorted(cids)
    rows, step = [], 200
    for i in range(0, len(ids), step):
        chunk = ",".join(str(x) for x in ids[i:i + step])
        txt = get(f"{BASE}/compound/cid/{chunk}/property/CanonicalSMILES/CSV")
        import io
        rows.append(pd.read_csv(io.StringIO(txt)))
        if (i // step) % 10 == 0:
            print(f"  SMILES {i+len(rows[-1])}/{len(ids)}", flush=True)
        time.sleep(0.25)          # PubChem просит не больше пяти запросов в секунду
    pd.concat(rows, ignore_index=True).to_csv(sm, index=False)
    print(f"сохранено: {sm}")


if __name__ == "__main__":
    main()
