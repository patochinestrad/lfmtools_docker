import streamlit as st
import os
import subprocess
from src import mTMAnalysis, helper_funcs
from stmol import showmol
import py3Dmol


st.title("mTM-Align RMSD Analysis")

st.subheader("Select one PDB file to be compared to the rest with mTM-Align.")

path = st.text_input("Input folder with PDB files")

if os.path.exists(path) == False:
    st.warning("Please enter a valid path")

else:
    with st.expander("PDB files in folder"):
        for i in [i for i in os.listdir(path) if i.endswith(".pdb")]:
            st.markdown(i[:-4])
    res_mtmalign = os.path.join(path, "res_mtmalign")
    output_dir = os.path.join(res_mtmalign, "all_vs_all_comparison")

    if st.button("Run All vs All comparison"):
        pdb_list = [i for i in os.listdir(path) if i.endswith(".pdb")]
        input_file = os.path.join(res_mtmalign, "input_file.txt")

        if os.path.exists(res_mtmalign):
            st.warning("RMSD Analysis already made for this folder")

        else:
            try:
                os.mkdir(res_mtmalign)
            except FileExistsError:
                st.warning("Output folder already exists")
            try:
                os.mkdir(output_dir)
            except FileExistsError:
                st.warning("Output folder for this comparison already")

            with open(input_file, "w") as file:
                for path in [os.path.join(path, i) for i in pdb_list]:
                    file.write(path + "\n")

            subprocess.run(
                [st.session_state["mtmAlign"], "-i", input_file, "-outdir", output_dir]
            )

    if st.checkbox("Show results"):
        rmsd_df, TM_df = mTMAnalysis.mTMAlign(output_dir)
        st.write("RMSD Matrix")
        st.dataframe(rmsd_df, use_container_width=True)
        st.write("TM Score Matrix")
        st.dataframe(TM_df, use_container_width=True)

        if st.checkbox("View superpositions"):
            choice = st.selectbox(
                "Choose PDB to display",
                tuple(
                    [
                        i
                        for i in os.listdir(
                            os.path.join(path + "/res_mtmalign/all_vs_all_comparison")
                        )
                        if i.endswith(".pdb")
                    ]
                ),
            )

            with open(
                os.path.join(path + "/res_mtmalign/all_vs_all_comparison", choice)
            ) as file:
                system = "".join([x for x in file])
                view = py3Dmol.view()
                view.addModelsAsFrames(system)
                view.setStyle({"model": -1}, {"cartoon": {"color": "spectrum"}})
                view.zoomTo()
                showmol(view)
