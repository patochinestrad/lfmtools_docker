import os
import glob
import re
import pandas as pd
import subprocess
import streamlit as st
from pandas_profiling import ProfileReport
from streamlit_pandas_profiling import st_profile_report
import shutil


def orderFpocketRes(path):
    res_fpocket = os.path.join(path, "res_fpocket")
    try:
        os.mkdir(res_fpocket)
    except FileExistsError:
        return st.warning("File exists error")

    for i in os.listdir(path):
        if i.endswith("_out"):
            shutil.move(os.path.join(path, i), res_fpocket)


def fpocketRun(pdb_file, alpha_min, alpha_max, clus_dist, clus_algo):
    subprocess.call(
        [
            "fpocket",
            "-f",
            str(pdb_file),
            "-m",
            str(alpha_min),
            "-M",
            str(alpha_max),
            "-D",
            str(clus_dist),
            "-C",
            str(clus_algo),
        ]
    )


def fpocketAnalysis(path):

    for i in os.listdir(path):
        res_parsed = os.path.join(path, "res_parsed/")
        try:
            os.mkdir(res_parsed)
        except FileExistsError:
            pass
        if i.endswith("_out"):
            res_folder = os.path.join(path, i)
            pockets = os.path.join(res_folder, "pockets/")
            proteinName = i[:-4]
            proteinPockets = []
            for pocket in glob.glob(pockets + "*.pdb"):
                pocketData = []
                pocketData.append(proteinName)
                pocketName = pocket[len(pockets) : -8]
                pocketData.append(pocketName)
                with open(pocket, "r") as file:
                    pocketResidues = []
                    for line in file.readlines():
                        if re.search("Pocket Score", line):
                            pocketScore = float(
                                re.findall(r"[-+]?(?:\d*\.\d+|\d+)", line)[-1]
                            )
                            pocketData.append(pocketScore)
                        if re.search("Drug Score", line):
                            drugScore = float(
                                re.findall(r"[-+]?(?:\d*\.\d+|\d+)", line)[-1]
                            )
                            pocketData.append(drugScore)
                        if re.search("Hydrophobicity Score", line):
                            hydrophobicityScore = float(
                                re.findall(r"[-+]?(?:\d*\.\d+|\d+)", line)[-1]
                            )
                            pocketData.append(hydrophobicityScore)
                        if re.search(" Pocket volume", line):
                            pocketVolume = float(
                                re.findall(r"[-+]?(?:\d*\.\d+|\d+)", line)[-1]
                            )
                            pocketData.append(pocketVolume)
                        if line.startswith("ATOM"):
                            columns = line.split()
                            resNum = columns[3] + "-" + columns[5]
                            pocketResidues.append(resNum)
                    pocketData.append(list(set(pocketResidues)))
                proteinPockets.append(pocketData)

            df = pd.DataFrame(
                columns=[
                    "ProteinName",
                    "PocketName",
                    "PocketScore",
                    "DrugScore",
                    "HydrophobicityScore",
                    "PocketVolume",
                    "PocketResidues",
                ],
                data=proteinPockets,
            )
            df.to_csv(res_parsed + proteinName + ".csv")


def fpocketDisplayResults(csv):
    df = pd.read_csv(csv, sep=",")
    st.dataframe(df.drop(["Unnamed: 0", "ProteinName"], axis=1))
    pr = df.drop(
        ["Unnamed: 0", "ProteinName", "PocketResidues", "PocketName"],
        axis=1,
    ).profile_report()
    st_profile_report(pr)
