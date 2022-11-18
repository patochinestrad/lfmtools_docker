import re
import pandas as pd



def mTMAlign(output_dir):
    

    rmsd_data = []
    rmsd_matrix = output_dir+'/pairwise_rmsd.txt'
    with open(rmsd_matrix) as rmsd_matrix:
        lines = rmsd_matrix.readlines()
        lines[0] = 'proteins.pdb' + lines[0]
        for i in lines:
            rmsd_data.append(list(filter(None,re.split('.pdb|\s+', i))))
    rmsd_df = pd.DataFrame(data=rmsd_data)
    rmsd_df = rmsd_df.set_index(0)

    TM_data = []
    TM_matrix = output_dir+'/pairwise_TMscore.txt'
    with open(TM_matrix) as TM_matrix:
        lines = TM_matrix.readlines()
        lines[0] = 'proteins.pdb' + lines[0]
        for i in lines:
            TM_data.append(list(filter(None,re.split('.pdb|\s+', i))))
    TM_df = pd.DataFrame(data=TM_data)
    TM_df = TM_df.set_index(0)

    return rmsd_df, TM_df
    