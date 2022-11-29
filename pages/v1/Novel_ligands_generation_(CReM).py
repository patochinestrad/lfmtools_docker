import streamlit as st
import os
from random import sample
from src import helper_funcs, autodockFuncs

st.title("Novel ligand generation with CReM")

st.subheader(
    "Input one or several ligands in pdbqt format and returns novel ligands in pdbqt format"
)

path = st.text_input("Input folder with ligands PDBQT files")
if os.path.exists(path) == False:
    st.warning("Please enter a valid path")

else:
    pdbqt = st.selectbox(
        "Choose seed ligand",
        [i for i in os.listdir(path) if i.endswith(".pdbqt")],
    )
    if pdbqt:
        grow_mode = st.radio("Select CReM generation method", ("Mutate", "Grow"))
        radius = st.select_slider("Select CReM radius", options=[1, 2, 3], value=2)
        n_to_keep = st.slider(
            "Select the number of new molecules to convert to .pdbqt. If selected number is greater than those generated, all generated molecules will be made.",
            min_value=0,
            max_value=100,
            value=0,
        )
        if st.button("Generate new ligands"):
            with st.spinner("Generating new molecules. Please wait."):
                outpath = os.path.join(path, pdbqt[:-6])
                pdbqt_path = os.path.join(path, pdbqt)
                try:
                    os.mkdir(outpath)
                except FileExistsError:
                    st.warning(f"New ligands already generated for ligand {pdbqt[:-6]}")

                mol = readPDBQT(pdbqt_path)
                new_mols = set(
                    generateNovelCompoundsList(mol, grow_mode=grow_mode, radius=radius)
                )
                n_newmols = len(new_mols)
                st.success(f"Generated {n_newmols} from {pdbqt[:-6]}")

            if n_to_keep >= n_newmols:
                st.warning(
                    f"Requested to make {n_to_keep} new ligands but was able to generate only {n_newmols}. Making {n_newmols} new molecules. Try changing radius to a lower value to generate more new ligands"
                )
                mols_to_make = new_mols
            else:
                mols_to_make = sample(new_mols, n_to_keep)

            compoundListToPDB(mols_to_make, outpath)

            outpdbqts = os.path.join(outpath, "pdbqt")
            try:
                os.mkdir(outpdbqts)
            except FileExistsError:
                st.warning(f"New ligands already generated for ligand {pdbqt[:-6]}")

            st.success(f"Converted {len(mols_to_make)} new_ligands to .pdbqt format.")
