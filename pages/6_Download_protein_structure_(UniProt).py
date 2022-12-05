import streamlit as st
from src import uniprotFuncs
from tempfile import TemporaryDirectory


st.title("Download protein structures from UniProt")

st.subheader(
    "Download protein structures from UniProt according to different parameters set by the user"
)

upid = st.text_input("UniProt ID")
if upid:
    with TemporaryDirectory() as tempDir:
        response = uniprotFuncs.getEntryTXT(upid)
        if not isinstance(response, Exception):
            pdbdf = uniprotFuncs.getPDBTable(response)
            st.write("PDB entries listed under this UniProtID:")
            st.dataframe(pdbdf, use_container_width=True)
            select_pdb = st.selectbox(
                "Select PDB entry",
                pdbdf.PDB,
            )
            if select_pdb:
                select_chain = st.selectbox(
                    "Select entry Chain",
                    pdbdf.loc[pdbdf["PDB"] == select_pdb]["Chains"].explode(),
                )
            st.write(f"You chose entry {select_pdb}, chain {select_chain}")

            ligands = uniprotFuncs.listLigands(select_pdb, select_chain)
            if ligands:
                hetatoms_to_keep = st.multiselect(
                    "Select Non-protein elements to keep",
                    ligands,
                )
            else:
                hetatoms_to_keep = False

            if hetatoms_to_keep:
                list_of_strings = [f"_{i}" for i in hetatoms_to_keep]
                string = "".join(list_of_strings)
                filename = f"{select_pdb}{select_chain}{string}.pdb"
            else:
                filename = f"{select_pdb}{select_chain}.pdb"
            st.download_button(
                "Download PDB Chain",
                data=open(
                    uniprotFuncs.downloadPDB(
                        select_pdb, select_chain, hetatoms_to_keep, tempDir
                    )
                ).read(),
                file_name=filename,
                mime="text/plain",
            )
