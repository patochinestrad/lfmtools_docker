import streamlit as st
import os
from src import autodockFuncs, helper_funcs
from tempfile import TemporaryDirectory
from random import sample
import zipfile
import io

st.title("Novel ligand generation with CReM")

st.subheader(
    "Input one ligand in pdbqt format to generate novel ligands using [CReM Database SA2](http://www.qsar4u.com/pages/crem.php) "
)

file = st.file_uploader("Upload PDBQT file", type="pdbqt", accept_multiple_files=False)


if file:
    grow_mode = st.radio("Select CReM generation method", ("Mutate", "Grow"))
    radius = st.select_slider("Select CReM radius", options=[1, 2, 3], value=2)
    n_to_keep = st.slider(
        "Select the number of new molecules to convert to .pdbqt. If selected number is greater than those generated, all generated molecules will be made.",
        min_value=0,
        max_value=100,
        value=0,
    )

    with TemporaryDirectory() as tempDir:
        st.write(tempDir)

        # '''
        # Save input files to tempDir to work on them
        # '''

        with open(os.path.join(tempDir, file.name), "wb") as f:
            f.write(file.getbuffer())
        st.write(os.listdir(tempDir))
        if st.button("Generate novel ligands"):
            with st.spinner("Generating new molecules. Please wait."):
                outpath = os.path.join(tempDir, file.name[:-6])
                os.mkdir(outpath)
                pdbqt_output = os.path.join(outpath, "pdbqts")
                os.mkdir(pdbqt_output)
                pdbqt_path = os.path.join(tempDir, file.name)

                mol = autodockFuncs.readPDBQT(pdbqt_path)
                new_mols = autodockFuncs.generateNovelCompoundsList(
                    mol, grow_mode, radius
                )
                n_newmols = len(new_mols)

                if n_to_keep >= n_newmols:
                    st.warning(
                        f"Requested to make {n_to_keep} new ligands but was able to generate only {n_newmols}. Making {n_newmols} new molecules. Try changing radius to a lower value to generate more new ligands"
                    )
                    mols_to_make = new_mols
                else:
                    mols_to_make = sample(new_mols, n_to_keep)

                autodockFuncs.compoundListToPDB(mols_to_make, outpath)

                autodockFuncs.ligandToPDBQT(outpath, pdbqt_output)
                st.success("Done!")

        with io.BytesIO() as buffer:
            with zipfile.ZipFile(buffer, mode="w") as archive:
                helper_funcs.zipdir(tempDir, archive)
            buffer.seek(0)
            st.download_button(
                "Download results",
                data=buffer,
                file_name="novel_ligands.zip",
                mime="application/zip",
            )
