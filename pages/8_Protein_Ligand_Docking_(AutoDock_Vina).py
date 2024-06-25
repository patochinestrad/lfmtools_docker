import streamlit as st
import tempfile
from streamlit_molstar import st_molstar
from src import autodockFuncs, helper_funcs
from src import uniprotFuncs

st.title("Protein Ligand Docking (AutoDock Vina)")

st.subheader(
    "Input one receptor file in PDBQT format and one or more ligands in PDBQT format to perform the docking"
)

receptor = st.file_uploader("Upload receptor file", type="pdb")
ligands = st.file_uploader(
    "Upload ligands files", type="pdbqt", accept_multiple_files=True
)

center = [
    st.number_input("Center x", value=0.000),
    st.number_input("Center y", value=0.000),
    st.number_input("Center z", value=0.000),
]

offset = [
    st.slider("Offset x", -10.0, 10.0, 0.0),
    st.slider("Offset y", -10.0, 10.0, 0.0),
    st.slider("Offset z", -10.0, 10.0, 0.0),
]

center = [sum(x) for x in zip(center, offset)]
size = [
    st.number_input("Size x", value=1),
    st.number_input("Size y", value=1),
    st.number_input("Size z", value=1),
]


box = autodockFuncs.generate_pdb(center, size)
if receptor:
    if st.button("Render receptor and box"):
        with tempfile.NamedTemporaryFile(
            delete=False, mode="w", encoding="utf-8", suffix=".pdb"
        ) as temp_file:
            temp_file.write(box)
            temp_file.write(receptor.read().decode("utf-8"))
            temp_file_path = temp_file.name

        with tempfile.NamedTemporaryFile(
            delete=False, mode="w", encoding="utf-8", suffix=".pdb"
        ) as receptorfile:
            receptorfile.write(receptor.read().decode("utf-8"))
            receptor_file_path = receptorfile.name

        st_molstar([temp_file_path], height="1080px")

        def format_center_and_size(center, size):
            labels = ["center_x", "center_y", "center_z", "size_x", "size_y", "size_z"]
            values = center + size
            return "\n".join(
                f"{label} = {value}" for label, value in zip(labels, values)
            )

        st.download_button(
            "Download configuration file",
            format_center_and_size(center, size),
            "config.txt",
            "txt",
        )
