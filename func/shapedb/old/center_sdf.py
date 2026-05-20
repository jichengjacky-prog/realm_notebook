import os,sys

my_file = sys.argv[1]
out_file = my_file.split(".sdf")[0] + "_centered.sdf"


from rdkit import Chem
import numpy as np

suppl = Chem.SDMolSupplier(my_file, removeHs=False)

writer = Chem.SDWriter(out_file)

for mol in suppl:
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    centroid = coords.mean(axis=0)
    for i in range(mol.GetNumAtoms()):
        pos = conf.GetAtomPosition(i)
        conf.SetAtomPosition(i, pos - centroid)
    writer.write(mol)
writer.close()
