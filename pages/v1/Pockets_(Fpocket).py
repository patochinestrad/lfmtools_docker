import streamlit as st
import os
import subprocess
import shutil
import pandas as pd
import pandas_profiling
from streamlit_pandas_profiling import st_profile_report
from src import fpocketAnalysis, helper_funcs


st.title("Fpocket Pocket Analysis")

st.subheader("Upload one or several PDB files to be analyzed with Fpocket.")
st.caption("The recommended values for each parameter are set as defaults.")
path = st.text_input("Input folder with PDB files")
if os.path.exists(path) == False:
    st.warning("Please enter a valid path")

else:
    filenames = helper_funcs.pdb_selector(path)

    if filenames:
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
        if st.button("Run"):
            with st.spinner("Analyzing pockets. Please wait."):
                for i in filenames:
                    subprocess.call(
                        [
                            "fpocket",
                            "-f",
                            str(i),
                            "-m",
                            str(alpha_sphere_values[0]),
                            "-M",
                            str(alpha_sphere_values[1]),
                            "-D",
                            str(clustering_distance),
                            "-C",
                            str(clus_algo_dict[clustering_algorithm]),
                        ]
                    )
            st.success("Done")
            res_fpocket = os.path.join(path, "res_fpocket")
            try:
                os.mkdir(res_fpocket)
            except FileExistsError:
                pass
            for i in os.listdir(path):
                if i.endswith("_out"):
                    shutil.move(os.path.join(path, i), res_fpocket)
            fpocketAnalysis.fpocketAnalysis(res_fpocket)

        if st.checkbox("Display results"):
            res_parsed = os.path.join(path, "res_fpocket/res_parsed")
            if os.path.exists(res_parsed):
                show = st.radio("Choose protein to show", os.listdir(res_parsed))
                df = pd.read_csv(os.path.join(res_parsed, show), sep=",")
                st.dataframe(df.drop(["Unnamed: 0", "ProteinName"], axis=1))
                pr = df.drop(
                    ["Unnamed: 0", "ProteinName", "PocketResidues", "PocketName"],
                    axis=1,
                ).profile_report()

                if st.checkbox("Display profile report"):
                    st_profile_report(pr)
            else:
                st.warning("No Fpocket results present in the current path.")
