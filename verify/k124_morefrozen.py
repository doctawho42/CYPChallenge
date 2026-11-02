"""Three more frozen representations, already on disk, through the probe device --- and a criterion
with a prediction written down before the run.

WHERE THIS COMES FROM. Item 346 found that the sub-floor members ADD: `chemeleon` alone is +0.0033
of macro rank and the gbm head on the 600 columns is +0.0023, while together they are +0.0047 with
the pair at -0.0085, both above their floors at sign 4/4. Item 345 then bought a criterion with a
negative: the organisers' `rdkit2d` encoder, run through the same device, is **-0.0089 at sign 0/4**,
because it was pretrained to predict the very RDKit descriptors our DESC block already contains. So
the rule is not "pretrained encoders help" but:

    a frozen checkpoint earns a member slot when its pretraining task is ACTIVITY-like and is not
    already present in our features.

THREE CANDIDATES THAT COST NOTHING. `data/` already holds `emb_deepchem_chemberta_77m_mtr.npz`
(384 dims), `emb_ibm_molformer_xl_both_10pct.npz` (768) and
`emb_seyonec_pubchem10m_smiles_bpe_450k.npz` (768), each with `train` and `test`. They were closed
by items 61/154 --- but as the SUBSTITUTION arm of `src/ablpre.py` ("эмбеддинг один", a
HistGradientBoosting head replacing our feature matrix), never as a RidgeCV probe whose prediction
enters the mean as its own member. That is the same distinction item 324 drew when it reopened
`chemprop_medium` after items 61/117/154 had closed embeddings-as-columns. `seyonec` has zero
mentions in this file at all.

THE PREDICTION, fixed before the measurement. `chemberta_77m_mtr` is pretrained by multi-task
regression on computed molecular descriptors, so by the criterion above it should fail the way
`rdkit2d` failed. `molformer_xl` and `pubchem10m_smiles_bpe` are masked language models over SMILES:
neither descriptors nor activity, so the criterion is SILENT on them and whatever they do extends it
rather than tests it. Writing this down is the point --- a criterion that only ever explains results
after the fact explains nothing.

THE ORDER CONTROL, because these files carry no `molecule_name`. For `molformer` the committed
`results/preds/oof_pre_ibm_molformer_xl_both_10pct.json` holds the "эмбеддинг один" predictions, and
re-running `ablpre`'s own head on the stored array reproduces them at **0.000e+00** on CYP2C9 against
a shifted control at 2.158 --- so that array is in `rows.csv` order and the pipeline's contract
holds. The same reproduction is run here for each candidate rather than assumed from the family.

    uv run python verify/k124_morefrozen.py --seeds 0,1,2,3

Writes results/logs/k124_morefrozen.json. Exit: 0 ran, 2 a control failed, 3 an input is missing.
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
import pandas as pd
from sklearn.ensemble import HistGradientBoostingRegressor

from cypsplit import butina_folds, fold_digest
from submit import CYPS

import k106_probe as K
from k106_probe import score, members, probe_oof, combine_add, load_emb

LOG = RES + "logs/k124_morefrozen.json"
FLOOR = 0.0036
GOLDEN = "2d93c19815e14261"
# (tag, dims, what the criterion predicts)
CANDS = [
    ("deepchem_chemberta_77m_mtr", "дескрипторы (MTR) -> ПРЕДСКАЗАН ПРОВАЛ"),
    ("ibm_molformer_xl_both_10pct", "языковая модель -> критерий МОЛЧИТ"),
    ("seyonec_pubchem10m_smiles_bpe_450k", "языковая модель -> критерий МОЛЧИТ"),
]
ABLPRE_KW = dict(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
                 l2_regularization=1.0, random_state=0)


def order_control(tag, E, rows, tr, fold):
    """Reproduce `src/ablpre.py`'s own committed arm from the stored array, on one enzyme.

    These files carry no `molecule_name`, so the order cannot be checked against a stored name
    list. It can be checked against a committed NUMBER: if the rows were permuted, the head would
    not reproduce. CYP2C9 is used because it has the fewest labels and so is the cheapest.
    """
    p = _pl.Path(RES + "preds/oof_pre_%s.json" % tag)
    if not p.exists():
        print("      %-36s committed-предсказаний нет, порядок НЕ проверен" % tag, flush=True)
        return None
    ref = json.load(open(p))["preds"].get("0|эмбеддинг один|CYP2C9")
    if ref is None:
        print("      %-36s в committed-файле нет нужного арма" % tag, flush=True)
        return None
    col = "CYP2C9_pIC50_direct_inhibition"
    m = tr[col].notna().to_numpy()
    y = tr.loc[m, col].to_numpy()
    Ei, fi = E[m], fold[m]
    q = np.zeros_like(y)
    for f in range(5):
        trn, te = fi != f, fi == f
        if te.sum() == 0:
            continue
        q[te] = HistGradientBoostingRegressor(**ABLPRE_KW).fit(Ei[trn], y[trn]).predict(Ei[te])
    r = np.asarray(ref, float)
    d = float(np.abs(q - r).max())
    off = float(np.abs(q[:-1] - r[1:]).max())
    print("      %-36s порядок: max|Δ| %.2e, сдвинутый контроль %.3f -> %s"
          % (tag, d, off, "ПОДТВЕРЖДЁН" if d < 1e-9 and off > 1e-3 else "НЕ ПОДТВЕРЖДЁН"),
          flush=True)
    return d < 1e-9 and off > 1e-3


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default="0")
    a = ap.parse_args(argv)
    seeds = [int(x) for x in a.seeds.split(",") if x.strip()]

    rows = pd.read_csv(D + "rows.csv")
    z = np.load(D + "feats.npz")
    XFULL = np.hstack([z["FP"], z["DESC"], z["MECH"]])
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr[f"{c}_pIC50_direct_inhibition"].to_numpy(float) for c in CYPS], 1)
    LO = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_low"].to_numpy(float) for c in CYPS], 1)
    HI = np.stack([tr[f"{c}_pIC50_direct_inhibition_conf_high"].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)
    K.S.LO, K.S.HI = LO, HI

    emb = {}
    for tag, _ in CANDS:
        p = _pl.Path(D + "emb_%s.npz" % tag)
        if not p.exists():
            print("нет %s" % p)
            return 3
        E = np.load(p)["train"].astype(np.float64)
        if E.shape[0] != len(rows):
            print("ОТКАЗ: %s имеет %d строк против %d" % (tag, E.shape[0], len(rows)))
            return 2
        emb[tag] = E

    out = json.load(open(LOG)) if _pl.Path(LOG).exists() else {}
    out.setdefault("meta", {}).update(floor=FLOOR,
                                      predictions={t: w for t, w in CANDS})
    out.setdefault("seeds", {})

    for s in seeds:
        fold, _ = butina_folds(list(rows.SMILES), seed=s)
        dg = fold_digest(fold)
        print("\n  сид %d: дайджест %s" % (s, dg), flush=True)
        if s == 0 and dg != GOLDEN:
            print("ОТКАЗ: разбиение сдвинулось, сид 0 дал %s" % dg)
            return 2

        if s == 0:
            for tag, _ in CANDS:
                order_control(tag, emb[tag], rows, tr, fold)

        pdz = members(s, XFULL, Y, M, fold, LO, HI)
        probe = [probe_oof(load_emb(rows), Y, M, fold, e)[0] for e in range(len(CYPS))]
        chem = json.load(open(RES + "preds/k122_chemeleon_probe_seed%d.json" % s))
        if chem.get("digest") != dg:
            print("ОТКАЗ: кэш chemeleon на другом разбиении")
            return 2
        chem = [np.asarray(v, float) for v in chem["P"]]
        g = json.load(open(RES + "preds/k114_heads_seed%d.json" % s))
        if g.get("digest") != dg:
            print("ОТКАЗ: кэш k114 на другом разбиении")
            return 2
        gb = [np.asarray(v, float) for v in g["heads"]["gbm"]]

        kept = [[Q[e] for k, Q in pdz if K._keep(e, k)] for e in range(len(CYPS))]
        base = [np.mean(kept[e] + [probe[e], chem[e], gb[e]], 0) for e in range(len(CYPS))]
        s_base = score(base, Y, LO, HI, M, fold)
        six = combine_add(pdz, probe)
        s_six = score(six, Y, LO, HI, M, fold)
        print("      КОНТРОЛЬ: шесть %.4f, лучший состав (восемь) %.4f, +%.4f (пункт 346: +0.0047)"
              % (s_six[0], s_base[0], s_base[0] - s_six[0]), flush=True)

        res = {}
        for tag, pred in CANDS:
            new = [probe_oof(emb[tag], Y, M, fold, e)[0] for e in range(len(CYPS))]
            sn = score(new, Y, LO, HI, M, fold)
            with_it = [np.mean(kept[e] + [probe[e], chem[e], gb[e], new[e]], 0)
                       for e in range(len(CYPS))]
            sw = score(with_it, Y, LO, HI, M, fold)
            cr = {nm: float(np.mean([np.corrcoef(new[e] - Y[M[:, e], e],
                                                 ref[e] - Y[M[:, e], e])[0, 1]
                                     for e in range(len(CYPS))]))
                  for nm, ref in (("зонд", probe), ("chemeleon", chem), ("gbm", gb),
                                  ("состав", base))}
            res[tag] = {"solo": sn, "gain_rank": sw[0] - s_base[0],
                        "gain_pair": sw[1] - s_base[1], "corr": cr, "prediction": pred}
            print("      %-36s одиночн. %.4f | Δранг к восьми %+.4f | Δпара %+.4f | "
                  "corr с составом %.3f"
                  % (tag, sn[0], sw[0] - s_base[0], sw[1] - s_base[1], cr["состав"]), flush=True)

        out["seeds"][str(s)] = {"six": s_six, "base8": s_base, "arms": res}
        _pl.Path(RES + "logs").mkdir(parents=True, exist_ok=True)
        json.dump(out, open(LOG, "w"), ensure_ascii=False, indent=1)

    print("\n" + "=" * 100)
    print("ТРИ ВЛОЖЕНИЯ С ДИСКА, УСТРОЙСТВОМ ЗОНДА, сверх состава из восьми")
    print("=" * 100)
    print("  %-36s %9s %9s %7s %9s  %s"
          % ("кандидат", "одиночн.", "Δранг", "знак", "Δпара", "предсказание критерия"))
    for tag, pred in CANDS:
        g = [out["seeds"][str(s)]["arms"][tag]["gain_rank"] for s in seeds
             if str(s) in out["seeds"]]
        p = [out["seeds"][str(s)]["arms"][tag]["gain_pair"] for s in seeds
             if str(s) in out["seeds"]]
        so = [out["seeds"][str(s)]["arms"][tag]["solo"][0] for s in seeds
              if str(s) in out["seeds"]]
        if not g:
            continue
        print("  %-36s %9.4f %+9.4f %5d/%d %+9.4f  %s"
              % (tag, float(np.mean(so)), float(np.mean(g)),
                 int(sum(1 for x in g if x > 0)), len(g), float(np.mean(p)), pred))
    print("\n  пол ранга %.4f, пол ST-RAE 0.007" % FLOOR)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
