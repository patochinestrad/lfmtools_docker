from rdkit import Chem
from oddt.toolkits.extras.rdkit import *
from crem.crem import mutate_mol, grow_mol
import shutil
import os
import subprocess

"""
Functions to generate novel ligands from one .pdbqt file.
"""


def readPDBQT(file):
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


def compoundListToPDB(mol_list, outdir_name):
    for idx, i in enumerate(mol_list):
        with open(os.path.join(outdir_name, f"newmol_{idx}.pdb"), "w") as pdb:
            mol = Chem.MolFromSmiles(i)
            Chem.AllChem.EmbedMolecule(
                mol, useExpTorsionAnglePrefs=True, useBasicKnowledge=True
            )
            Chem.AllChem.UFFOptimizeMolecule(mol)
            pdb_block = Chem.MolToPDBBlock(mol)
            pdb.write(pdb_block)
            pdb.close()


def pdbsToPDBQT(indir_name, outdir_name):
    for i in [
        os.path.join(indir_name, i) for i in os.listdir(indir_name) if i.endswith("pdb")
    ]:
        subprocess.call(
            [
                "python",
                "src/prepare_ligand4.py",
                "-v",
                "-l",
                i,
                "-o",
                i[:-4] + ".pdbqt",
                "-A",
                "bonds_hydrogens",
            ]
        )
    for i in [os.path.join(i) for i in os.listdir(indir_name) if i.endswith("pdbqt")]:
        shutil.move(i, outdir_name)


# testmol = readPDBQT("/home/patricio/Documents/streamlit_test/ligands_test/C8.pdbqt")
# new_mols = generateNovelCompoundsList(testmol, "Mutate", 2)
# compoundListToPDB(
#     new_mols, "/home/patricio/Documents/streamlit_test/ligands_test/C8_pdbs/"
# )

pdbsToPDBQT(
    "/home/patricio/Documents/streamlit_test/ligands_test/C8_pdbs/",
    "/home/patricio/Documents/streamlit_test/ligands_test/C8_pdbs/C8_pdbqts/",
)
