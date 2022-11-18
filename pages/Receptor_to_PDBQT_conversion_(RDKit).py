import streamlit as st
import os
from src import RDKitPDBQTConversion, helper_funcs
import shutil
import glob

st.title("Protein preparation for docking")

st.subheader("Upload one or several PDB files to be converted to PDBQT with RDKit")

path = st.text_input("Input folder with PDB files")
if os.path.exists(path) == False:
    st.warning("Please enter a valid path")


filenames = helper_funcs.pdb_selector(path)

if st.checkbox("Convert to PDBQT"):
    if st.button("Run"):
        res_pdbqt = os.path.join(path, "res_pdbqt_conversion")
        try:
            os.mkdir(res_pdbqt)
        except FileExistsError:
            st.warning("Output folder already exists")
        with st.spinner("Converting to pdbqt. Please wait."):
            for i in filenames:
                RDKitPDBQTConversion.pdbqtConversion(i)
            for i in glob.glob(path + "/*.pdbqt"):
                try:
                    shutil.move(i, res_pdbqt)
                except:
                    st.warning(
                        "PDBQT Conversion for %s was already done. Press Ok to continue with the next one."
                        % i
                    )
                    os.remove(i)
            st.success("Done")
