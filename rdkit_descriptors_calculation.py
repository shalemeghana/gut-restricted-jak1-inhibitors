from rdkit import Chem
from rdkit.Chem import Descriptors
import pandas as pd

df = pd.read_csv('7908_NTox.csv')
df['Label'].replace({'Non-toxic':0,'Toxic':1}, inplace=True)
Smiles = df['Smiles']
Label = df['Label']
mols = [Chem.MolFromSmiles(smi) for smi in Smiles]
def getMolDescriptors(mol, missingVal=None):
    res={}
    for nm, fn in Descriptors._descList:
        try:
            val = fn(mol)
        except:
            import traceback
            traceback.print_exc()
            val = missingVal
        res[nm] = val
    return res        
all_desc = [getMolDescriptors(m) for m in mols]
desc = pd.DataFrame(all_desc)
desc.head()
