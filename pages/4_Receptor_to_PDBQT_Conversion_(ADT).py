import streamlit as st
import os
from src import autodockFuncs, helper_funcs
from tempfile import TemporaryDirectory
import zipfile
import io


st.title("Protein preparation for docking")

st.subheader(
    "Upload one or several PDB files to be converted to PDBQT with AutoDock Tools"
)

files = st.file_uploader("Upload PDB files", type="pdb",
                         accept_multiple_files=True)
if files:
    with TemporaryDirectory() as tempDir:
        # '''
        # Save input files to tempDir to work on them
        # '''

        for file in files:
            with open(os.path.join(tempDir, file.name), "wb") as f:
                f.write(file.getbuffer())

        if st.button("Convert all to PDBQT format"):
            with st.spinner("Converting structures. Please wait"):
                outdir = os.path.join(tempDir, "converted_structures")
                os.mkdir(outdir)
                autodockFuncs.receptorToPDBQT(
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
