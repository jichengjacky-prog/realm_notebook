import argparse
import os
from pathlib import Path

from rdkit import Chem
import numpy as np

parser = argparse.ArgumentParser(
    description="Convert all MOL2 files in a shapedb directory to centered SDF files."
)
parser.add_argument(
    "shapedb_dir",
    nargs="?",
    default="out/shapedb",
    help="Directory containing MOL2 files to convert (default: out/shapedb)",
)
args = parser.parse_args()

shapedb_dir = Path(args.shapedb_dir)
if not shapedb_dir.exists():
    raise FileNotFoundError(f"Shapedb directory not found: {shapedb_dir}")

for root, dirs, files in os.walk(shapedb_dir):
    for file in files:
        if file.endswith(".mol2"):
            mol2_path = Path(root) / file

            mol = Chem.MolFromMol2File(str(mol2_path), removeHs=False)
            if mol is None:
                print(f"Failed to read {mol2_path}")
                continue

            sdf_path = mol2_path.with_name(mol2_path.stem + "_centered.sdf")

            conf = mol.GetConformer()
            coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
            centroid = coords.mean(axis=0)
            for i in range(mol.GetNumAtoms()):
                pos = conf.GetAtomPosition(i)
                conf.SetAtomPosition(i, pos - centroid)

            writer = Chem.SDWriter(str(sdf_path))
            writer.write(mol)
            writer.close()

            print(f"Processed: {mol2_path} -> {sdf_path}")

print("Done!")
