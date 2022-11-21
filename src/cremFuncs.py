from rdkit import Chem
from oddt.toolkits.extras.rdkit import *
from crem.crem import mutate_mol, grow_mol
import os

"""
Functions to generate novel ligands from one .pdbqt file.
"""


def smilesConversion(file):
    with open(file, "r") as pdbqt:
        mol = MolFromPDBQTBlock(pdbqt.read())
        return mol


def generateNovelCompoundsList(mol, grow_mode, radius):
    modes = ["Mutate", "Grow"]
    if grow_mode == "Mutate":
        new_mols = list(
            mutate_mol(mol, db_name="./data/crem_db/crem_db_sa2.db", radius=radius)
        )
        return new_mols
    elif grow_mode == "Grow":
        new_mols = list(
            grow_mol(mol, db_name="./data/crem_db/crem_db_sa2.db", radius=radius)
        )
        return new_mols


def compoundListToPDBQT(mol_list, outdir_name):
    for idx, i in enumerate(mol_list):
        with open(os.path.join(outdir_name, f"newmol_{idx}.pdbqt"), "w") as pdbqt:
            mol = Chem.MolFromSmiles(i)
            Chem.AllChem.EmbedMolecule(mol)
            pdbqt_block = MolToPDBQTBlock(mol, addHs=True)
            pdbqt.write(pdbqt_block)
            pdbqt.close()
