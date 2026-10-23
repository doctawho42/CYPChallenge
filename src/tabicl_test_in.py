"""Build the TabICL TEST bundle: train and test halves of both representations, with the control
that the 750 unnamed rows of `surrogate_adme` are in the test table's order.

Item 374 measured the two in-context members out of fold and said in as many words that nothing
could ship until test predictions existed. This assembles their input.

THE PROBLEM THIS FILE EXISTS TO SOLVE. `emb_chemprop_medium_test.npz` carries `molecule_name`, so
its order is checked element by element against the test table exactly as `submit.py:_probe_load`
checks it, rolled control included. `emb_chemprop_surrogate_adme.npz` carries no names at all --- its
`train` and `test` arrays are positional. The train half has a control already: item 349's and
`k121.align_control`'s, a probe that collapses from 0.5533 to -0.0192 when the matrix is rolled by
one row. **The test half has no labels, so that control cannot be run on it**, and writing 750
predictions against a row order nobody checked is how item 302 overwrote four working receptor files.

THE CONTROL USED INSTEAD, and why it must succeed. Both test arrays describe the same 750 molecules
through two chemprop encoders. If their rows correspond, a molecule's nearest NEIGHBOUR inside the
test set should largely agree between the two embeddings; if one array is shifted, the agreement
collapses to chance. The statistic is the median rank, in `medium`'s neighbour ordering, of the
neighbour `surrogate` picks --- 1 is perfect agreement, 375 is chance at n=750. Three readings are
taken and all three are required:

    обучение   the same statistic on the TRAIN halves, where BOTH orders are independently
               verified (medium by name, surrogate by align_control).  This is the calibration:
               it says what "correct" looks like for this pair of encoders, so the test figure is
               read against a measurement rather than against a guess.
    тест       the figure that decides.  It must be close to the calibration.
    сдвиг      the test arrays with surrogate rolled by one row.  It must collapse.  A control
               that cannot fail is not a control (CLAUDE.md), and rolling leaves the matrix, its
               scale and its neighbour structure untouched while breaking only correspondence.

Writes data/tabicl_test_in.npz. Nothing is fitted here.
"""

import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "verify"))
from cyppaths import D, RES, tutorial
tutorial()

import numpy as np
import pandas as pd

from cypsplit import butina_folds, fold_digest
from submit import CYPS
from k121_chemeleon import align_control

E_TR = D + "emb_chemprop_medium.npz"
E_TE = D + "emb_chemprop_medium_test.npz"
S_BOTH = D + "emb_chemprop_surrogate_adme.npz"
OUT = D + "tabicl_test_in.npz"
GOLDEN = "2d93c19815e14261"


def named_order(path, names, what):
    """`submit.py:_probe_load`'s check, reused rather than reimplemented: names must match element
    for element AND the same comparison run on the names rolled by one must FAIL."""
    z = np.load(path, allow_pickle=True)
    E = z["train"].astype(np.float64)        # the test file stores its array under `train` too
    nm = np.asarray(z["molecule_name"], object)
    r = np.asarray(names, object)
    if E.shape[0] != len(r):
        raise SystemExit("%s: %d строк против %d в таблице" % (path, E.shape[0], len(r)))
    ok, rolled = bool((nm == r).all()), bool((np.roll(nm, 1) == r).all())
    if not ok or rolled:
        raise SystemExit("%s: порядок не совпал (совпал=%s, сдвинутый=%s)" % (path, ok, rolled))
    if not np.isfinite(E).all():
        raise SystemExit("%s: во вложении есть не-конечные значения" % path)
    print("    %-9s %s, имена совпали: %s; сдвинутый контроль (ожидается False): %s"
          % (what, E.shape, ok, rolled), flush=True)
    return E


def _nn(A):
    """Index of each row's nearest OTHER row, and the full ordering, by cosine distance."""
    Z = A - A.mean(0)
    Z /= np.maximum(np.linalg.norm(Z, axis=1, keepdims=True), 1e-12)
    G = Z @ Z.T
    np.fill_diagonal(G, -np.inf)
    return G


def nn_agree(A, B):
    """Median rank, in A's neighbour ordering, of the nearest neighbour B picks. 1 is perfect."""
    GA, GB = _nn(A), _nn(B)
    nnB = GB.argmax(1)
    order = np.argsort(-GA, axis=1)
    rank = np.empty(len(A), int)
    for i in range(len(A)):
        rank[i] = int(np.where(order[i] == nnB[i])[0][0]) + 1
    return float(np.median(rank)), float((rank == 1).mean())


def main():
    rows = pd.read_csv(D + "rows.csv")
    te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")
    tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
            .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
    Y = np.stack([tr["%s_pIC50_direct_inhibition" % c].to_numpy(float) for c in CYPS], 1)
    M = ~np.isnan(Y)

    print("  ВЛОЖЕНИЯ С ИМЕНАМИ")
    Etr = named_order(E_TR, rows.Molecule_Name, "medium")
    Ete = named_order(E_TE, te.Molecule_Name, "medium тест")

    print("\n  surrogate_adme: имён в файле нет")
    z = np.load(S_BOTH, allow_pickle=True)
    Str, Ste = z["train"].astype(np.float64), z["test"].astype(np.float64)
    if Str.shape[0] != len(rows) or Ste.shape[0] != len(te):
        raise SystemExit("surrogate: %s / %s против %d / %d строк"
                         % (Str.shape, Ste.shape, len(rows), len(te)))
    fold0, _ = butina_folds(list(rows.SMILES), seed=0)
    if fold_digest(fold0) != GOLDEN:
        raise SystemExit("разбиение сдвинулось: %s" % fold_digest(fold0))
    if not align_control(Str, Y, M, fold0, e=0):
        raise SystemExit("ОТКАЗ: порядок обучающей половины surrogate не подтверждён")
    print("    обучающая половина: контроль метками пройден (пункт 349, k121.align_control)")

    print("\n  СОГЛАСИЕ СОСЕДЕЙ между medium и surrogate --- контроль порядка ТЕСТА")
    cal_med, cal_top1 = nn_agree(Etr, Str)
    tst_med, tst_top1 = nn_agree(Ete, Ste)
    rol_med, rol_top1 = nn_agree(Ete, np.roll(Ste, 1, axis=0))
    print("    %-10s медианный ранг %7.1f   совпал первым %5.1f%%" % ("обучение", cal_med,
                                                                     100 * cal_top1))
    print("    %-10s медианный ранг %7.1f   совпал первым %5.1f%%" % ("тест", tst_med,
                                                                      100 * tst_top1))
    print("    %-10s медианный ранг %7.1f   совпал первым %5.1f%%" % ("сдвиг", rol_med,
                                                                      100 * rol_top1))
    print("    случайно при n=%d ожидалось бы около %.0f" % (len(te), len(te) / 2))
    if not rol_med > 10 * tst_med:
        raise SystemExit("ОТКАЗ: сдвинутый контроль не развалился (%.1f против %.1f) --- "
                         "статистика не проверяет соответствие строк" % (rol_med, tst_med))
    if not tst_med < 5 * cal_med:
        raise SystemExit("ОТКАЗ: согласие на тесте (%.1f) далеко от калибровки (%.1f) --- "
                         "порядок 750 строк surrogate не подтверждён" % (tst_med, cal_med))
    print("    ПОДТВЕРЖДЕНО: тест близок к калибровке, сдвиг развалился")

    np.savez_compressed(OUT, E_tr=Etr.astype(np.float32), E_te=Ete.astype(np.float32),
                        S_tr=Str.astype(np.float32), S_te=Ste.astype(np.float32),
                        Y=Y, cyps=np.array(CYPS),
                        mol_test=np.asarray(te.Molecule_Name, object))
    import os
    print("\n  записано %s, %.1f МБ" % (OUT, os.path.getsize(OUT) / 1e6))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
