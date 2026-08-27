"""Насколько часто в обучающей выборке встречается тот класс азота, на котором
правило pKa провалилось на кофеине (амид-фланкированный имидазол/пиридин)."""
import sys as _sys, pathlib as _pl
_sys.path.insert(0, str(_pl.Path(__file__).resolve().parents[1]))
from cyppaths import D, RES, tutorial
tutorial()
import numpy as np, pandas as pd, sys
sys.path.insert(0,'" + D + "')
from rdkit import Chem, RDLogger
RDLogger.DisableLog('rdApp.*')
import pka as PK
tr=pd.read_csv('" + D + "cyp-challenge-TRAIN_inhibition.csv')
mols=[Chem.MolFromSmiles(s) for s in tr.SMILES]
ok=[m for m in mols if m]
print(f"молекул: {len(ok)}")
cls=[]; val=[]
for m in ok:
    r=PK.most_basic_pka(m)
    if isinstance(r,tuple): v,c=r[0],r[1]
    else: v,c=r,"?"
    cls.append(c); val.append(v)
cls=pd.Series(cls); val=np.array(val,float)
print("\nраспределение назначенных классов:")
print(cls.value_counts().to_string())
print(f"\nоценка pKa: медиана {np.nanmedian(val):.2f}, доля протонированных при 7.4: {np.mean(1/(1+10**(7.4-val))>0.5):.3f}")
# сколько молекул имеют имидазол/пиридиновый азот рядом с карбонилом (случай кофеина)
patt=Chem.MolFromSmarts("[n;$(n1cncc1),$(n1ccccc1)]~[#6]~[CX3]=[OX1]")
amide_ring=Chem.MolFromSmarts("[nX3]([#6])[CX3]=[OX1]")
c1=sum(1 for m in ok if m.HasSubstructMatch(patt))
c2=sum(1 for m in ok if m.HasSubstructMatch(amide_ring))
print(f"\nмолекул с ароматическим N через атом от карбонила: {c1} ({100*c1/len(ok):.1f}%)")
print(f"молекул с N-ациламидом в кольце (мотив кофеина):     {c2} ({100*c2/len(ok):.1f}%)")
risky=cls.isin(["imidazole","pyridine"]).to_numpy()
print(f"молекул, у которых самый основной центр отнесён к имидазолу/пиридину: {risky.sum()} ({100*risky.mean():.1f}%)")
both=np.array([ok[i].HasSubstructMatch(patt) for i in range(len(ok))]) & risky
print(f"пересечение (класс имидазол/пиридин И карбонил рядом) -- зона риска: {both.sum()} ({100*both.mean():.1f}%)")
