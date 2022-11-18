from rdkit import Chem
from oddt.toolkits.extras.rdkit import *

def pdbqtConversion(file):
    mol = Chem.MolFromPDBFile(file)
    with open(file[:-4]+'.pdbqt', 'w') as out:
        pdbqt = MolToPDBQTBlock(mol, flexible=False, addHs=True, computeCharges=True)
        out.write(pdbqt)
        out.close()
