"""Ковариантный сдвиг обучение -> тест в пространстве самих предсказаний.

Одна и та же модель применяется к обоим наборам, и сравниваются распределения её
собственных выходов. Это честнее прокси по ближайшему соседу (k4): отображение
признаки -> предсказание одно и то же, поэтому разница в выходах — это разница входов,
а не разница качества прокси. Признаки теста строятся feats.build() с переиндексацией
по именам колонок, как в src/submit.py.

БЛОК 2 ЗДЕСЬ ПРОВАЛИЛСЯ, и оставлен как провалившийся. Веса важности из логистической
регрессии «обучение против теста» на всех 2295 признаках вырождаются: AUC разделения
0.817, медианный вес 0.001, эффективный размер выборки 116 из 4905, оптимумы упираются
в края сетки. Любая метрика под такими весами описывает сотню молекул, а не набор.
Числа блока 3 поэтому недействительны и в выводы не идут. Рабочая версия — k6_shift1d.py,
где отношение плотностей оценивается по одной оси, вдоль которой аффинная перекалибровка
и работает."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1] / "src"))
from cyppaths import D, RES, tutorial
tutorial()
import time, numpy as np, pandas as pd, json
from evaluation.custom_scoring_functions import rae_soft_threshold_absolute_error as strae
from sklearn.ensemble import HistGradientBoostingRegressor
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_predict
from sklearn.preprocessing import StandardScaler
from scipy.stats import ks_2samp
import feats as F

CYPS = ["CYP1A2", "CYP2C9", "CYP2D6", "CYP3A4"]
rows = pd.read_csv(D + "rows.csv")
tr = (pd.read_csv(D + "cyp-challenge-TRAIN_inhibition.csv")
        .set_index("Molecule_Name").loc[rows.Molecule_Name].reset_index())
z = np.load(D + "feats.npz")
X = np.hstack([z["FP"], z["DESC"], z["MECH"]])
dnames = pd.read_csv(D + "desc_names.csv", header=None)[0].tolist()
mnames = pd.read_csv(D + "mech_names.csv", header=None)[0].tolist()
te = pd.read_csv(D + "cyp-challenge-TEST-BLINDED.csv")

cache = _pl.Path(D) / "test_feats.npz"
if cache.exists():
    zz = np.load(cache)
    XT = np.hstack([zz["FP"], zz["DESC"], zz["MECH"]])
else:
    t0 = time.time()
    tFP, tD, tM, ok = F.build(list(te.SMILES), dnames, mnames)
    assert len(ok) == len(te), "часть тестовых SMILES не разобралась — строки бы поехали"
    np.savez_compressed(cache, FP=tFP, DESC=tD.to_numpy(np.float32),
                        MECH=tM.to_numpy(np.float32))
    XT = np.hstack([tFP, tD.to_numpy(np.float32), tM.to_numpy(np.float32)])
    print(f"признаки теста построены за {time.time()-t0:.0f} с")
assert XT.shape[1] == X.shape[1], (XT.shape, X.shape)
oof = json.load(open(RES + "preds/oof.json"))

print()
print("=" * 84)
print("1. Распределение предсказаний: вне фолда на обучении против теста")
print("=" * 84)
print(f"{'фермент':8s} {'ybar обуч':>9s} {'сред. OOF':>10s} {'сред. тест':>11s} {'сдвиг':>7s} "
      f"{'sd OOF':>7s} {'sd тест':>8s} {'KS p':>9s} {'доля>5.5':>16s}")
PT = {}
for c in CYPS:
    col = f"{c}_pIC50_direct_inhibition"; m = tr[col].notna().to_numpy()
    y = tr.loc[m, col].to_numpy(); p = np.asarray(oof[f"FP+DESC+MECH|{c}"])
    mdl = HistGradientBoostingRegressor(max_iter=300, learning_rate=0.06, max_leaf_nodes=31,
                                        l2_regularization=1.0, random_state=0).fit(X[m], y)
    pt = mdl.predict(XT); PT[c] = pt
    print(f"{c:8s} {y.mean():9.3f} {p.mean():10.3f} {pt.mean():11.3f} {pt.mean()-p.mean():+7.3f} "
          f"{p.std():7.3f} {pt.std():8.3f} {ks_2samp(pt,p).pvalue:9.2e} "
          f"{np.mean(p>5.5)*100:6.1f}%->{np.mean(pt>5.5)*100:5.1f}%")
np.savez_compressed(D + "test_pred.npz", **PT)
print("предсказания на тесте сохранены в data/test_pred.npz — их читает k6_shift1d.py")

print()
print("=" * 84)
print("2. ПРОВАЛИВШИЙСЯ БЛОК: веса важности на всех признаках")
print("=" * 84)
sc = StandardScaler().fit(np.vstack([X, XT]))
Xs, XTs = sc.transform(X), sc.transform(XT)
Z = np.vstack([Xs, XTs]); lab = np.r_[np.zeros(len(Xs)), np.ones(len(XTs))]
pr = cross_val_predict(LogisticRegression(C=0.05, max_iter=3000), Z, lab, cv=5,
                       method="predict_proba")[:, 1]
auc = float(np.mean(pr[lab == 1][:, None] > pr[lab == 0][None, :]))
w = pr[:len(Xs)] / (1 - pr[:len(Xs)] + 1e-9)
w = np.clip(w, 0, np.quantile(w, 0.99)); w /= w.mean()
neff = w.sum() ** 2 / (w ** 2).sum()
print(f"AUC разделения обучение/тест = {auc:.3f}   (0.5 = наборы неразличимы)")
print(f"веса: медиана {np.median(w):.3f}, 90-й процентиль {np.quantile(w,0.90):.2f}, "
      f"эффективный размер {neff:.0f} из {len(w)}")
print(f"\nВЫВОД: {neff:.0f} из {len(w)} — это не перевзвешивание, а подвыборка из сотни")
print("молекул. Взвешенная метрика под такими весами не измеряет ничего переносимого;")
print("результат отброшен, продолжение в k6_shift1d.py.")
