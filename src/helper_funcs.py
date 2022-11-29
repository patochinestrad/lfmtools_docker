import os
import zipfile
import streamlit as st


def pdb_selector(folder_path):
    filenames = [i for i in os.listdir(folder_path) if i.endswith(".pdb")]
    if filenames:
        select_all = st.checkbox("Select all files")
        if select_all:
            return [os.path.join(folder_path, i) for i in filenames]
        selected_filenames = st.multiselect("Select a file", filenames)
        return [os.path.join(folder_path, i) for i in selected_filenames]
    else:
        return st.write("No PDB file in this directory")


def pdbqt_selector(folder_path):
    filenames = [i for i in os.listdir(folder_path) if i.endswith(".pdbqt")]
    if filenames:
        select_all = st.checkbox("Select all files")
        if select_all:
            return [os.path.join(folder_path, i) for i in filenames]
        selected_filenames = st.multiselect("Select a file", filenames)
        return [os.path.join(folder_path, i) for i in selected_filenames]
    else:
        return st.write("No PDB file in this directory")


def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            ziph.write(
                os.path.join(root, file),
                os.path.relpath(os.path.join(root, file), os.path.join(path, "..")),
            )
        for directory in dirs:
            ziph.write(os.path.join(root, directory))
