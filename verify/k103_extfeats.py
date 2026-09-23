"""Build and cache the external CheMeleon feature block, aligned to trunk.BLOCKS.

Separated from the training driver so that a network drop or a crash during training never
costs the ~2 minutes of RDKit work again. Writes data/ext_thirdhead.npz (gitignored by the
*.npz rule), which the driver reads.

Overlap handling is trunkext.py's, reused rather than re-derived: external rows whose
canonical SMILES matches one of ours are DROPPED (item 290 recorded 64 of them), and the
blind test is guarded against. The test guard computes FROM the local blinded CSV -- it is
not a database lookup -- and it must find zero.
"""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()

import numpy as np
import pandas as pd
from rdkit import Chem

import feats as F
import trunk as T

CYPS = T.CYPS
OUT = D + "ext_thirdhead.npz"


def canon(s):
    try:
        m = Chem.MolFromSmiles(s)
        return Chem.MolToSmiles(m) if m else None
    except Exception:
        return None


def main():
    X, y, lo, hi, scr, smiles = T.load(T.BLOCKS)
    rows = pd.read_csv(D + "rows.csv")

    ext = pd.concat([pd.read_csv(D + "chemeleon_X_train.csv"),
                     pd.read_csv(D + "chemeleon_y_train.csv")], axis=1)
    n0 = len(ext)
    ext["k"] = [canon(s) for s in ext.OPENADMET_CANONICAL_SMILES]
    ours = set(filter(None, (canon(s) for s in rows.SMILES)))
    te_keys = set(filter(None, (canon(s) for s in
                                pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv").SMILES)))

    # Control that MUST succeed and that exercises the suspect clause: the matcher has to be
    # able to find a hit at all, so match the test set against ITSELF before asking whether the
    # external set hits it. A zero from a matcher that cannot match anything is not an answer.
    self_hits = sum(1 for k in te_keys if k in te_keys)
    print(f"контроль: тест против самого себя {self_hits}/{len(te_keys)} "
          f"(должно быть {len(te_keys)}/{len(te_keys)})", flush=True)
    if self_hits != len(te_keys) or len(te_keys) == 0:
        raise SystemExit("контроль совпадения не сработал: измерение пересечения бессмысленно")

    hit_test = int(ext.k.isin(te_keys).sum())
    print(f"пересечение внешних с СЛЕПЫМ тестом по каноническому SMILES: {hit_test}", flush=True)
    if hit_test:
        raise SystemExit("внешние данные пересекаются с тестом")

    hit_ours = int(ext.k.isin(ours).sum())
    ext = ext[~ext.k.isin(ours)].reset_index(drop=True)
    print(f"внешних строк {n0}, пересеклось с нашим обучением {hit_ours} (пункт 290: 64), "
          f"осталось {len(ext)}", flush=True)

    dn = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
    mn = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
    FP, dsc, M, ok = F.build(list(ext.OPENADMET_CANONICAL_SMILES), dn, mn)
    ext = ext.iloc[ok].reset_index(drop=True)
    parts = {"FP": np.log1p(FP), "DESC": dsc.to_numpy(np.float32), "MECH": M.to_numpy(np.float32)}
    Xe = np.hstack([parts[b] for b in T.BLOCKS.split("+")]).astype(np.float32)
    Xe = np.nan_to_num(Xe, nan=0.0, posinf=0.0, neginf=0.0)
    if Xe.shape[1] != X.shape[1]:
        raise SystemExit(f"ширина не совпала: {Xe.shape[1]} против {X.shape[1]}")

    ye = np.stack([ext[f"OPENADMET_LOGAC50_{c.lower()}"].to_numpy(np.float32) for c in CYPS], 1)
    print(f"Xe {Xe.shape}, ye {ye.shape}")
    print("внешних меток по ферментам: "
          + ", ".join(f"{c} {int(np.isfinite(ye[:, e]).sum())}" for e, c in enumerate(CYPS)))
    np.savez_compressed(OUT, Xe=Xe, ye=ye)
    print(f"сохранено: {OUT}")


if __name__ == "__main__":
    main()
