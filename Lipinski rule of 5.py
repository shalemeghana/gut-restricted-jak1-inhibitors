#import rdkit packages
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import Lipinski
from rdkit.Chem import Crippen

#import required packages for dataframe and numericals
import pandas as pd
import numpy as np

#import datasets
data = pd.read_csv('output_24.csv')
Smiles = data['Smiles']
#import sdf file
#sdf_file = Chem.SDMolSupplier('chembl_33.sdf')

# Get the molecules in the file
"""molecules = []
for mol in sdf_file:
    molecules.append(mol)
"""
#create blank values to store values
emw , hbad, hbdd, clogp, emw1, hbad1, hbdd1, clogp1, error_indices = ([] for i in range(9))

#apply loop for calculating descriptors
for i, smiles in enumerate(Smiles):
    try:
        mol = Chem.MolFromSmiles(smiles)
        R1 = Lipinski.NumHAcceptors(mol)
        R2 = Lipinski.NumHDonors(mol)
        R3 = Crippen.MolLogP(mol)
        R4 = Descriptors.ExactMolWt(mol)
        hbad.append(R1 > 10) 
        hbad1.append(R2)
        hbdd.append(R2 > 5)
        hbdd1.append(R2)
        clogp.append(R3 >= 5)
        clogp1.append(R3)
        emw.append(R4>= 500)
        emw1.append(R4)
    except Exception as e:
        error_indices.append(i)
        print(f"Error occurred for SMILES {smiles} at index {i}: {e}")
print("Indices of problematic SMILES strings:", error_indices)
#create a dataframe with the calculated datas    
results = pd.DataFrame({
    'Smiles': Smiles,
    'Exact_Mol_Wt_Values': emw1,
    'Exact_Mol_Wt_Dev' : emw,
    'H_bond_Acc_Values':hbad1,
    'H_bond_Acc' : hbad,
    'H_bond_Don_Values' : hbdd1,
    'H_bond_Don' : hbdd,
    'CLogP_Values': clogp1,
    'CLogP': clogp,
})

#out the dataframe frame containing the calculated values with either deviating or not
results.to_csv('Lipinski Rule Deviations results.csv', index=False)

#load the Dataframe containing the the calculated values for analyzing the frequency either deviating or not
Data = pd.read_csv('Lipinski Rule Deviations results.csv')
MW = Data.Exact_Mol_Wt_Dev.value_counts().rename_axis('MW').reset_index()
HBA = Data.H_bond_Acc.value_counts().rename_axis('HBA').reset_index()
HBD = Data.H_bond_Don.value_counts().rename_axis('HBD').reset_index()
CLP = Data.CLogP.value_counts().rename_axis('cLogP').reset_index()
df_concat = pd.concat([MW, HBA, HBD, CLP], axis=1)

#write ot the csv file of the calculated frequency of each rule
df_concat.to_csv('Frequency.csv', index=False)

#create a dataframe to detect the molecules deviating two rules at a time
value = True
results1 = results[results['Exact_Mol_Wt_Dev'] & results['CLogP']== value]
results1.to_csv('Exact_Mol_Wt_Dev_&_CLogP.csv')
results2 = results[results['Exact_Mol_Wt_Dev'] & results['H_bond_Don']== value]
results2.to_csv('Exact_Mol_Wt_Dev_&_H_bond_Don.csv')
results3 = results[results['Exact_Mol_Wt_Dev'] & results['H_bond_Acc']== value]
results3.to_csv('Exact_Mol_Wt_Dev_&_H_bond_Acc.csv')
results4 = results[results['Exact_Mol_Wt_Dev']==value]
results4.to_csv('Exact_Mol_Wt_Dev.csv')
results5 = results[results['CLogP']==value]
results5.to_csv('CLogP.csv')
results6 = results[results['H_bond_Don']==value]
results6.to_csv('H_bond_Don.csv')
results7 = results[results['H_bond_Acc']==value]
results7.to_csv('H_bond_Acc.csv')
