import os
import streamlit as st


def pdb_selector(folder_path):
    if not os.path.exists(folder_path):
        return st.write("Please enter a valid path")
    filenames = [i for i in os.listdir(folder_path) if i.endswith(".pdb")]
    if filenames:
        select_all = st.checkbox("Select all files")
        if select_all:
            return [os.path.join(folder_path, i) for i in filenames]
        selected_filenames = st.multiselect("Select a file", filenames)
        return [os.path.join(folder_path, i) for i in selected_filenames]
    else:
        return st.write("No PDB file in this directory")
