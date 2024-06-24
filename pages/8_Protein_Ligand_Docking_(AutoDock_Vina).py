import streamlit as st
import st_speckmol
from src import autodockFuncs, helper_funcs
from src import uniprotFuncs

st.title("Protein Ligand Docking (AutoDock Vina)")

st.subheader(
    "Input one receptor file in PDBQT format and one or more ligands in PDBQT format to perform the docking"
)

receptor = st.file_uploader("Upload receptor file", type="pdbqt")
ligands = st.file_uploader(
    "Upload ligands files", type="pdbqt", accept_multiple_files=True
)

if receptor and ligands:
    uniprotFuncs.speckmolView(receptor)
    """
    Hacer una caja con el centro y los lados del grid
    """
    center = st.number_input
    st_speckmol(
        """C        0.000000        0.000000        0.000000
H        0.000000        0.000000        1.089000
H        1.026719        0.000000       -0.363000
H       -0.513360       -0.889165       -0.363000
H       -0.513360        0.889165       -0.363000"""
    )
