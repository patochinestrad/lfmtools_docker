import streamlit as st
import os, zipfile, io
from tempfile import TemporaryDirectory
from src.fingerprintSimilarity import *
import matplotlib.pyplot as plt

st.title(
    "Ligand similarity comparison using RDKit and the Tanimoto coefficient of similarity"
)

st.subheader(
    "Input any number of ligands in pdb format to generate different types of similarity matrices"
)

files = st.file_uploader("Upload PDB file", type="pdb", accept_multiple_files=True)

if files:
    with TemporaryDirectory() as tempDir:

        # """
        # Save input files to tempDir to work on them
        # """

        for file in files:
            with open(os.path.join(tempDir, file.name), "wb") as f:
                f.write(file.getbuffer())

        fp = st.radio(
            "Choose the type of fingerprint to calculate.",
            [
                "Morgan 2D Fingerprints",
                "Morgan 2D Fingerprints (Murcko Scaffold)",
                "MACCS Keys Fignerprints",
            ],
        )
        annotate = st.checkbox("Show values in cells?")
        if st.button("Compute fingerprints"):
            with st.spinner("Computing fingerprints. Please wait."):
                moldict = molsToDict(
                    [os.path.join(tempDir, i) for i in os.listdir(tempDir)]
                )
                match fp:
                    case "Morgan 2D Fingerprints":
                        fpdict = morganFPDict(moldict)
                        fptype = fp
                    case "Morgan 2D Fingerprints (Murcko Scaffold)":
                        fpdict = murckoScaffoldFPDict(moldict)
                        fptype = fp
                    case "MACCS Keys Fignerprints":
                        fpdict = maccsKeysFPDict(moldict)
                        fptype = fp

                sim_list = tanimotoSimilarity(fpdict)
                fig, ax = plt.subplots()
                ax = plotSimilarityHeatmap(sim_list, fptype=fptype, annotate=annotate)
                st.write(fig)

                if fig:
                    img = io.BytesIO()
                    fn = f"{fptype} heatmap.png"

                    plt.savefig(fn, dpi=300, bbox_inches="tight")

                    with open(fn, "rb") as img:
                        st.download_button(
                            "Download image", data=img, file_name=fn, mime="image/png"
                        )
