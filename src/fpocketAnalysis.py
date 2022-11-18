import os
import glob
import re
import pandas as pd

def fpocketAnalysis(path):

    for i in os.listdir(path):
        res_parsed = os.path.join(path, 'res_parsed/')
        try:
            os.mkdir(res_parsed)
        except FileExistsError:
            pass
        if i.endswith('_out'):
            res_folder = os.path.join(path, i)
            pockets = os.path.join(res_folder, 'pockets/')
            proteinName = i[:-4]
            proteinPockets = []
            for pocket in glob.glob(pockets+'*.pdb'):
                pocketData = []
                pocketData.append(proteinName)
                pocketName = pocket[len(pockets):-8]
                pocketData.append(pocketName)
                with open(pocket, 'r') as file:
                    pocketResidues = []
                    for line in file.readlines():
                        if re.search('Pocket Score', line):
                            pocketScore = float(re.findall(r'[-+]?(?:\d*\.\d+|\d+)', line)[-1])
                            pocketData.append(pocketScore)
                        if re.search('Drug Score', line):
                            drugScore = float(re.findall(r'[-+]?(?:\d*\.\d+|\d+)', line)[-1])
                            pocketData.append(drugScore)
                        if re.search('Hydrophobicity Score', line):
                            hydrophobicityScore = float(re.findall(r'[-+]?(?:\d*\.\d+|\d+)', line)[-1])
                            pocketData.append(hydrophobicityScore)
                        if re.search(' Pocket volume', line):
                            pocketVolume = float(re.findall(r'[-+]?(?:\d*\.\d+|\d+)', line)[-1])
                            pocketData.append(pocketVolume)
                        if line.startswith('ATOM'):
                            columns = line.split()
                            resNum = columns[3]+'-'+columns[5]
                            pocketResidues.append(resNum)
                    pocketData.append(list(set(pocketResidues)))
                proteinPockets.append(pocketData)
            
            
            df = pd.DataFrame(columns=['ProteinName', 'PocketName', 'PocketScore', 'DrugScore', 'HydrophobicityScore', 'PocketVolume', 'PocketResidues'], data=proteinPockets)     
            df.to_csv(res_parsed+proteinName+'.csv')
