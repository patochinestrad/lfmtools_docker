from rdkit import Chem
from crem.crem import mutate_mol, grow_mol
import shutil
import os
import subprocess

"""
Functions to generate novel ligands from one .pdbqt file.
"""


def readPDB(file):
    with open(file, "r") as pdbqt:
        mol = Chem.MolFromPDBBlock(pdbqt.read())
        return mol


def generateNovelCompoundsList(mol, grow_mode, radius):
    if grow_mode == "Mutate":
        new_mols = list(
            mutate_mol(
                mol,
                db_name="crem_db_sa2.db",
                radius=radius,
            )
        )
        return new_mols
    elif grow_mode == "Grow":
        new_mols = list(
            grow_mol(
                mol,
                db_name="crem_db_sa2.db",
                radius=radius,
            )
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


def ligandToPDBQT(indir_name, outdir_name):
    for i in [
        os.path.join(indir_name, i) for i in os.listdir(indir_name) if i.endswith("pdb")
    ]:
        print(i)
        subprocess.call(
            [
                "python",
                "src/prepare_ligand4.py",
                "-l",
                i,
                "-o",
                i[:-4] + ".pdbqt",
                "-A",
                "bonds_hydrogens",
            ]
        )
    for i in [
        os.path.join(indir_name, i)
        for i in os.listdir(indir_name)
        if i.endswith("pdbqt")
    ]:
        shutil.move(i, outdir_name)


def receptorToPDBQT(indir_name, outdir_name):
    for i in [
        os.path.join(indir_name, i) for i in os.listdir(indir_name) if i.endswith("pdb")
    ]:
        subprocess.call(
            [
                "python",
                "src/prepare_receptor",
                "-r",
                i,
                "-A",
                "bonds_hydrogens",
                "-o",
                i[:-4] + ".pdbqt",
            ]
        )
    for i in [
        os.path.join(indir_name, i)
        for i in os.listdir(indir_name)
        if i.endswith("pdbqt")
    ]:
        shutil.move(i, outdir_name)


def generate_pdb(center, size):
    # Extract center and size
    cx, cy, cz = center
    sx, sy, sz = size

    # Calculate the half sizes
    hx, hy, hz = sx / 2, sy / 2, sz / 2

    # Define the vertices based on center and half sizes
    vertices = [
        (cx - hx, cy - hy, cz - hz),
        (cx - hx, cy - hy, cz + hz),
        (cx - hx, cy + hy, cz - hz),
        (cx + hx, cy - hy, cz - hz),
        (cx + hx, cy - hy, cz + hz),
        (cx + hx, cy + hy, cz - hz),
        (cx - hx, cy + hy, cz + hz),
        (cx + hx, cy + hy, cz + hz),
    ]

    # Define the center point
    center_point = (cx, cy, cz)

    # Start writing the PDB content
    pdb_content = []

    # Write the vertices
    for i, vertex in enumerate(vertices, start=1):
        pdb_content.append(
            f"HETATM{str(i).rjust(5)} PS{i}  PSD P   1    {vertex[0]:8.3f}{vertex[1]:8.3f}{vertex[2]:8.3f}  0.00  0.00          PSDOPS"
        )

    # Write the center point
    pdb_content.append(
        f"HETATM{str(len(vertices) + 1).rjust(5)} PS{len(vertices) + 1}  PSD P   1    {center_point[0]:8.3f}{center_point[1]:8.3f}{center_point[2]:8.3f}  0.00  0.00          PSDOPS"
    )

    # Write the CONECT records
    conect_records = [
        "CONECT    1    2    4    3",
        "CONECT    2    1    5    7",
        "CONECT    3    1    6    7",
        "CONECT    4    1    5    6",
        "CONECT    5    2    4    8",
        "CONECT    6    3    4    8",
        "CONECT    7    2    3    8",
        "CONECT    8    5    6    7",
    ]

    pdb_content.extend(conect_records)
    pdb_content.append("END")

    # Join the PDB content into a single string
    pdb_string = "\n".join(pdb_content)

    return pdb_string
    # with open("temp/grid_box.pdb", "w") as f:
    # f.write(pdb_string)
