from rdkit import Chem, DataStructs
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem import MACCSkeys
import os
import itertools
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


def readPDB(file, check_file=True):
    with open(file, "r") as pdb:
        mol = Chem.MolFromPDBBlock(pdb.read(), sanitize=check_file)
        return mol


def molsToDict(mol_list, semisanitize=False):
    if semisanitize:
        new_mols = [readPDB(i, check_file=True) for i in mol_list]
        new_mols = [m.UpdatePropertyCache(strict=False) for m in new_mols]
        for m in new_mols:
            m = Chem.SanitizeMol(
                m, Chem.SanitizeFlags.SANITIZE_KEKULIZE, catchErrors=True
            )
    else:
        new_mols = [readPDB(i) for i in mol_list]

    none_mols = {
        mol_list[i]: new_mols[i] for i in range(len(mol_list)) if new_mols[i] is None
    }
    mol_dict = {
        mol_list[i]: new_mols[i]
        for i in range(len(mol_list))
        if new_mols[i] is not None
    }
    return mol_dict, none_mols


def morganFPDict(mol_dict):
    fp_dict = {
        key: Chem.AllChem.GetMorganFingerprintAsBitVect(
            value,
            2,
            useFeatures=True,
            nBits=1024,
        )
        for key, value in mol_dict.items()
    }

    return fp_dict


def murckoScaffoldFPDict(mol_dict):
    fp_dict = {
        key: Chem.AllChem.GetMorganFingerprintAsBitVect(
            MurckoScaffold.GetScaffoldForMol(value),
            2,
            useFeatures=True,
            nBits=1024,
        )
        for key, value in mol_dict.items()
    }

    return fp_dict


def maccsKeysFPDict(mol_dict):
    fp_dict = {key: MACCSkeys.GenMACCSKeys(
        value) for key, value in mol_dict.items()}

    return fp_dict


def tanimotoSimilarity(fp_dict):
    keys_sim = []
    for pair in itertools.product(fp_dict.keys(), repeat=2):
        key1 = pair[0]
        key_name1 = pair[0][:-4]
        key2 = pair[1]
        key_name2 = pair[1][:-4]
        fpsim = DataStructs.FingerprintSimilarity(fp_dict[key1], fp_dict[key2])
        keys_sim.append([key_name1, key_name2, fpsim])

    return keys_sim


def plotSimilarityHeatmap(similarity_list, fptype, annotate=False):
    simdf = pd.DataFrame(similarity_list)
    simdf[0] = [os.path.basename(i) for i in simdf[0]]
    simdf[1] = [os.path.basename(i) for i in simdf[1]]
    simdf = simdf.set_index([0, 1])[2].unstack()
    if annotate:
        ax = sns.heatmap(
            simdf,
            vmin=0,
            vmax=1,
            cbar_kws={"label": "Tanimoto Similarity Coefficient"},
            annot=True,
        )
    else:
        ax = sns.heatmap(
            simdf,
            vmin=0,
            vmax=1,
            cbar_kws={"label": "Tanimoto Similarity Coefficient"},
        )
    ax.set(
        title=f"Tanimoto Coefficient of Similarity for {fptype} fingerprint",
        xlabel="Molecule",
        ylabel="Molecule",
    )
    return ax
