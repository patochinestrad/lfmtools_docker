import urllib
import streamlit as st
import pandas as pd
from prody import parsePDB, writePDB, pathPDBFolder
from tempfile import TemporaryDirectory
import os


def getEntryTXT(upid):
    try:
        url = urllib.request.urlopen(f"https://rest.uniprot.org/uniprotkb/{upid}.txt")
        response = url.read().decode("utf-8")
        return response
    except urllib.error.HTTPError as e:
        st.warning(
            f"""There was an error while trying to get your UniProt ID:
             *{e}*. 
            Maybe it is wrong?"""
        )
        return e


def getPDBTable(response):
    pdbtable = [
        i.split(";") for i in response.splitlines() if i.startswith("DR   PDB;")
    ]
    table = []
    for i in pdbtable:
        pdbname = i[1].strip(" ")
        technique = i[2]
        resolution = i[3]
        chains_coverage = i[4]
        chains = chains_coverage.split("=")[0].split("/")
        chains = [i.strip(" ") for i in chains]
        coverage = chains_coverage.split("=")[1].split("-")
        coverage = [i.strip(".") for i in coverage]
        table.append([pdbname, technique, resolution, chains, coverage])
    pdbdataframe = pd.DataFrame(
        table, columns=["PDB", "Technique", "Resolution [A]", "Chains", "Coverage"]
    )
    return pdbdataframe


def listLigands(entry, chainid):
    with TemporaryDirectory() as tempDir:
        pathPDBFolder(tempDir)
        ppdb = parsePDB(entry, chain=chainid)
        hetatm = ppdb.select("hetatm and not water")
        return list(set([i for i in hetatm.getResnames()]))


def downloadPDB(entry, chainid, hetatoms, dir):
    with TemporaryDirectory() as tempDir:
        pathPDBFolder(tempDir)
        ppdb = parsePDB(entry, chain=chainid)
        if hetatoms:
            list_of_strings = ["or hetatm resname " + i for i in hetatoms]
            string = " ".join(list_of_strings)
            protein = ppdb.select("protein not water " + string)

        else:
            protein = ppdb.select("protein and not water")
        return writePDB(os.path.join(dir, f"{entry}{chainid}.pdb"), protein)
