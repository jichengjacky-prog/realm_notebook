#script skeleton given by Summer to get stats like logp on ligands from SMILES strings
#this will ideally be used to predict toxicity correlation with logP in drugs ordered for screens
#feed the script a smiles string, and it will give you desired stats

from rdkit import Chem
from rdkit.Chem import Descriptors
import os,sys

smi = sys.argv[1]

# Function to compute key descriptors
def compute_descriptors(smiles):
	mol = Chem.MolFromSmiles(smiles)
	if mol is None:
		return None

	desc = {}
	desc['MW'] = Descriptors.MolWt(mol)
	desc['cLogP'] = Descriptors.MolLogP(mol)
	desc['TPSA'] = Descriptors.TPSA(mol)
	desc['HBD'] = Descriptors.NumHDonors(mol)
	desc['HBA'] = Descriptors.NumHAcceptors(mol)
	desc['RotatableBonds'] = Descriptors.NumRotatableBonds(mol)
	desc['AromaticRings'] = Descriptors.NumAromaticRings(mol)
	desc['Fsp3'] = Descriptors.FractionCSP3(mol)

	# Approximate logD at pH 7.4 using ionizable groups (weak estimate)
	# For neutral molecules logD ≈ logP
	desc['logD7.4'] = desc['cLogP'] # Simple approximation; more advanced methods require pKa

	return desc

# Example compounds
smiles_list = {
	"Suvorexant": "CC(C)C1=CC=C(C=C1)C2=CC=CC=C2C3=NC(=O)C4=CC=CC=C34", # example
	"Seltorexant": "CC1=CC=CC=C1C2=NC3=C(C=CC=C3)C(=O)N2C" # example
}

#for name, smi in smiles_list.items():
#	d = compute_descriptors(smi)
#	print(name, d)
d = compute_descriptors(smi)
print(d)
