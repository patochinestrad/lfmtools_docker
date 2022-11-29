import streamlit as st
import os
from src import fpocketAnalysis, helper_funcs
from tempfile import TemporaryDirectory
import zipfile
import io

st.title("Fpocket Pocket Analysis")

st.subheader("Upload one or several PDB files to be analyzed with Fpocket.")
st.caption("The recommended values for each parameter are set as defaults.")


files = st.file_uploader("Upload PDB files", type="pdb", accept_multiple_files=True)
if files:
    with TemporaryDirectory() as tempDir:

        alpha_sphere_values = st.slider(
            "Select Alfa Sphere minimum and maximum radius.",
            0.0,
            10.0,
            value=(3.2, 6.4),
        )
        clustering_distance = st.number_input("Select Clustering Distance", 2.4)
        clustering_algorithm = st.radio(
            "Select Clustering Method for grouping voronoi vertices",
            [
                "Single Linkage",
                "Complete Linkange",
                "Average Linkage",
                "Centroid Linkage",
            ],
        )
        clus_algo_dict = {
            "Single Linkage": "s",
            "Complete Linkage": "m",
            "Average Linkage": "a",
            "Centroid Linkage": "c",
        }

        # """
        # Save input files to tempDir to work on them
        # """

        for file in files:
            with open(os.path.join(tempDir, file.name), "wb") as f:
                f.write(file.getbuffer())

        # """
        # Zip temDir and download results
        # """

        if st.button("Run pocket detection"):
            with st.spinner("Detecting pockets. Please wait."):
                for i in os.listdir(tempDir):

                    fpocketAnalysis.fpocketRun(
                        os.path.join(tempDir, i),
                        alpha_sphere_values[0],
                        alpha_sphere_values[1],
                        clustering_distance,
                        clus_algo_dict[clustering_algorithm],
                    )
                fpocketAnalysis.fpocketAnalysis(tempDir)
                fpocketAnalysis.orderFpocketRes(tempDir)
                st.success("Done")

        with io.BytesIO() as buffer:

            with zipfile.ZipFile(buffer, mode="w") as archive:
                helper_funcs.zipdir(tempDir, archive)
            buffer.seek(0)
            st.download_button(
                "Download results",
                data=buffer,
                file_name="fpocket_results.zip",
                mime="application/zip",
            )

    if st.checkbox("Check to analyze your results if you already have them"):
        csv_files = st.file_uploader(
            "Upload CSV files to analyse", type="csv", accept_multiple_files=True
        )

        if csv_files:
            names = [i.name[:-4] for i in csv_files]
            csvfile = st.radio("Choose file to show analysis", names)
            for i in csv_files:
                if csvfile + ".csv" == i.name:
                    fpocketAnalysis.fpocketDisplayResults(i)
