import urllib
import streamlit as st
import pandas as pd
from prody import parsePDB, writePDB, pathPDBFolder
from tempfile import TemporaryDirectory
import os
import nglview as nv


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
        if hetatm:
            return list(set([i for i in hetatm.getResnames()]))
        else:
            return st.write("This entry has no ligand")


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


def displayPDB(pdbname, assembly_number):
    return st.components.v1.html(
        f"""
                                    <html lang="en">
                                        <head>
                                            <meta charset="utf-8" />
                                            <meta name="viewport" content="width=device-width, user-scalable=yes, minimum-scale=1.0, maximum-scale=10.0">
                                            
                                            
                                            
                                            <script src="https://cdn.jsdelivr.net/npm/babel-polyfill/dist/polyfill.min.js"></script>
                                            
                                            <script src="https://cdn.jsdelivr.net/npm/@webcomponents/webcomponentsjs/webcomponents-lite.js" charset="utf-8"></script>
                                            
                                            <script src="https://cdn.jsdelivr.net/npm/@webcomponents/webcomponentsjs/custom-elements-es5-adapter.js" charset="utf-8"></script>

                                            <link rel="stylesheet" type="text/css" href="https://www.ebi.ac.uk/pdbe/pdb-component-library/css/pdbe-molstar-3.1.0.css">
                                            <script type="text/javascript" src="https://www.ebi.ac.uk/pdbe/pdb-component-library/js/pdbe-molstar-component-3.1.0.js"></script>
                                            <style>
                                                #myViewer{{
                                                float:left;
                                                width:800px;
                                                height: 800px;
                                                position:relative;
                                                }}
                                            </style>
                                        </head>
                                        <body>

                                        
                                            
                                    
                                            
                                            <div id="myViewer">
                                                <pdbe-molstar id="pdbeMolstarComponent" molecule-id="{pdbname.lower()}" hide-controls="false" assembly-id="{assembly_number}"></pdbe-molstar>
                                            </div>
                                            
                                        </body>
                                    </html>
                                    """,
        width=900,
        height=1200,
    )
