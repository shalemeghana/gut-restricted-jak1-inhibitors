import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import MolStandardize
from rdkit.Chem.Draw import IPythonConsole
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Draw
df = pd.read_csv('Clean.csv')
smis = df['Smiles']
ms2 = []
for smi in smis:
    m = Chem.MolFromSmiles(smi,sanitize=False)
    ma = rdMolStandardize.Cleanup(m)
    mb = rdMolStandardize.DisconnectOrganometallics(ma)
    mc = Chem.MolStandardize.rdMolStandardize.Normalize(mb)
    md = Chem.MolStandardize.rdMolStandardize.Reionize(mc)
    me = MolStandardize.rdMolStandardize.RemoveFragments(md)
    mf = rdMolStandardize.FragmentParent(me)
    uncharger = rdMolStandardize.Uncharger()
    mg = uncharger.uncharge(mf)
    te = rdMolStandardize.TautomerEnumerator()
    mh = te.Canonicalize(mg)
    mh.UpdatePropertyCache(strict=False)
    mi = rdMolStandardize.Normalize(mh)
    ms2.extend([mi])
smiles = [Chem.MolToSmiles(mol) for mol in ms2]
