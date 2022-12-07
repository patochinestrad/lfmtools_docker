import streamlit as st
from src import uniprotFuncs, helper_funcs
from tempfile import TemporaryDirectory
import zipfile
import io

st.title("Download protein structures from UniProt")

st.subheader(
    "Download protein structures from UniProt according to different parameters set by the user"
)

upid = st.text_input("UniProt ID")
if upid:
    st.write(
        f"You can find the AlphaFold model for this UniProt entry [here](https://alphafold.ebi.ac.uk/entry/{upid})."
    )

    with TemporaryDirectory() as tempDir:
        response = uniprotFuncs.getEntryTXT(upid)
        if not isinstance(response, Exception):
            pdbdf = uniprotFuncs.getPDBTable(response)
            if not pdbdf.empty:
                st.write("PDB entries listed under this UniProtID:")
                st.dataframe(pdbdf, use_container_width=True)
                if st.button("Download all structures with their ligands"):
                    with st.spinner(
                        "Preparing files, it might take some time depending on the number of structures."
                    ):
                        for index, row in pdbdf.iterrows():
                            pdb = row["PDB"]
                            chains = row["Chains"]
                            for chain in chains:
                                uniprotFuncs.downloadPDB(
                                    pdb, chain, hetatoms=[], dir=tempDir
                                )

                    with io.BytesIO() as buffer:
                        with zipfile.ZipFile(buffer, mode="w") as archive:
                            helper_funcs.zipdir(tempDir, archive)
                        buffer.seek(0)
                        st.download_button(
                            "Download structures",
                            data=buffer,
                            file_name=f"{upid}_structures.zip",
                            mime="application/zip",
                        )

                else:
                    pdblist = pdbdf["PDB"].tolist()
                    pdblist.insert(0, "Select PDB")
                    select_pdb = st.selectbox(
                        "Select PDB entry",
                        pdblist,
                    )

                    if select_pdb == "Select PDB":
                        pass
                    else:
                        select_chain = st.selectbox(
                            "Select entry Chain",
                            pdbdf.loc[pdbdf["PDB"] == select_pdb]["Chains"].explode(),
                        )
                        st.write(f"You chose entry {select_pdb}, chain {select_chain}")
                        if st.checkbox("Display structure in visualizer"):
                            assembly_number = ord(select_chain) - 96
                            uniprotFuncs.displayPDB(select_pdb, assembly_number)
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
            else:
                st.warning("This UniProt ID has no X-ray or NMR structure solved.")
