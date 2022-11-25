from rdkit import Chem
from oddt.toolkits.extras.rdkit import *
from crem.crem import mutate_mol, grow_mol
import os

"""
Functions to generate novel ligands from one .pdbqt file.
"""


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
