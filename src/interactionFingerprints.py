from openbabel import openbabel as ob
import os
from rdkit.Chem import MolFromPDBFile, CombineMols, MolToPDBFILE


def convertDirToPDB(directory):
    obconversion = ob.OBConversion()
    obconversion.SetInAndOutFormats("pdbqt", "pdb")

    for file in [os.path.join(directory, file) for file in os.listdir(directory)]:
        obmol = ob.OBMol()
        obconversion.ReadFile(obmol, file)
        obconversion.WriteFile(obmol, file[:-2])


def receptorAndLigandsToComplex(receptor, ligands: list):
    receptor_mol = MolFromPDBFile(receptor, sanitize=False)
    complexes = []
    for ligand in ligands:
        ligand_mol = MolFromPDBFile(ligand, sanitize=False)
        complex = CombineMols(receptor_mol, ligand_mol)
        complexes.append(complex)

    for complex in complexes:


receptorAndLigandsToComplex("/home/patricio/Documents/streamlit_test/PLIP/6AZUB.pdb")
