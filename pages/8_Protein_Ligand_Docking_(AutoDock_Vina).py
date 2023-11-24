import streamlit as st
from src import autodockFuncs, helper_funcs
from src import uniprotFuncs

st.title("Ligand PDB to PDBQT conversion")

st.subheader(
    "Input one receptor file in PDBQT format and one or more ligands in PDBQT format to perform the docking"
)

receptor = st.file_uploader("Upload receptor file", type="pdbqt")
ligands = st.file_uploader(
    "Upload ligands files", type="pdbqt", accept_multiple_files=True
)

if receptor and ligands:
    if st.checkbox("Display receptor"):
        uniprotFuncs.speckmolView(receptor)
