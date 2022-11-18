import streamlit as st
import os
import subprocess
import shutil
import pandas as pd
import pandas_profiling
from streamlit_pandas_profiling import st_profile_report
from src import fpocketAnalysis, helper_funcs
from stmol import showmol
import py3Dmol

st.title("Fpocket Pocket Analysis")

st.subheader("Upload one or several PDB files to be analyzed with Fpocket")

path = st.text_input("Input folder with PDB files")
if os.path.exists(path) == False:
    st.warning("Please enter a valid path")


filenames = helper_funcs.pdb_selector(path)

if st.checkbox("Run Fpocket"):
    if st.button("Run"):
        with st.spinner("Analyzing pockets. Please wait."):
            for i in filenames:
                subprocess.run(["fpocket", "-f", i])
        st.success("Done")

if st.checkbox("Analyse Results"):
    if st.button("Analyse"):
        res_fpocket = os.path.join(path, "res_fpocket")
        try:
            os.mkdir(res_fpocket)
        except FileExistsError:
            pass
        for i in os.listdir(path):
            if i.endswith("_out"):
                shutil.move(os.path.join(path, i), res_fpocket)
        fpocketAnalysis.fpocketAnalysis(res_fpocket)

if st.checkbox("Display results"):
    res_parsed = os.path.join(path, "res_fpocket/res_parsed")
    if os.path.exists(res_parsed):
        show = st.radio("Choose protein to show", os.listdir(res_parsed))
        df = pd.read_csv(os.path.join(res_parsed, show), sep=",")
        st.dataframe(df.drop(["Unnamed: 0", "ProteinName"], axis=1))
        pr = df.drop(
            ["Unnamed: 0", "ProteinName", "PocketResidues", "PocketName"], axis=1
        ).profile_report()

        if st.checkbox("Display profile report"):
            st_profile_report(pr)
    else:
        st.write(
            "Please enter a valid path to the PDB folder with Fpocket results already"
        )
