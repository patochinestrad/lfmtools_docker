import os
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


# def zipdir(path, ziph):
#     # ziph is zipfile handle
#     for root, dirs, files in os.walk(path):
#         for file in files:
#             ziph.write(
#                 os.path.join(root, file),
#                 os.path.relpath(os.path.join(root, file), os.path.join(path, "..")),
#             )
#         for directory in dirs:
#             ziph.write(os.path.join(root, directory))
def zipdir(path, ziph):
    # ziph is zipfile handle
    for root, dirs, files in os.walk(path):
        for file in files:
            file_path = os.path.join(root, file)
            rel_path = os.path.relpath(file_path, path)
            ziph.write(file_path, rel_path)
        for directory in dirs:
            dir_path = os.path.join(root, directory)
            rel_path = os.path.relpath(dir_path, path)
            ziph.write(dir_path, rel_path)
