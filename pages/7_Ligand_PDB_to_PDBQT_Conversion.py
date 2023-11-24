import streamlit as st
from src import autodockFuncs, helper_funcs
import os
from tempfile import TemporaryDirectory
import io
import zipfile

st.title("Ligand PDB to PDBQT conversion")

st.subheader(
    "Input one or multiple ligand in pdb format to conver to PDBQT for posterior docking"
)

files = st.file_uploader("Upload PDB files", type="pdb", accept_multiple_files=True)

if files:
    with TemporaryDirectory() as tempDir:
        # '''
        # Save input files to tempDir to work on them
        # '''

        for file in files:
            with open(os.path.join(tempDir, file.name), "wb") as f:
                f.write(file.getbuffer())

        if st.button("Convert all to PDBQT format"):
            with st.spinner("Converting ligands. Please wait"):
                outdir = os.path.join(tempDir, "converted_ligands")
                os.mkdir(outdir)
                autodockFuncs.ligandToPDBQT(
                    tempDir,
                    outdir,
                )
                st.success("Done")

            with io.BytesIO() as buffer:
                with zipfile.ZipFile(buffer, mode="w") as archive:
                    helper_funcs.zipdir(tempDir, archive)
                buffer.seek(0)
                st.download_button(
                    "Download results",
                    data=buffer,
                    file_name="pdbqt_conversion.zip",
                    mime="application/zip",
                )
