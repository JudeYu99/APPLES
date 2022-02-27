#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@ author: Yu Zhu

@ Email: yzhu99@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

"""

### Introduction of NACEN_R.py
#
# > Invoke NACEN.R programme to execute the AAN construction and analysis in R. Then, this python script reads the output file for further analysis.
#
# @ This part of program is delicated for generating amino acid network and performing topology analysis with R package NACEN.
# @ A working installation of R >= 3.2.0 and NACEN package are required when running this sub-program.
#
# @ Note: NACEN package can be downloaded in the website: http://sysbio.suda.edu.cn/NACEN/NACEN.rar
# @ Note: NACEN requires DSSP programme for building AAN which can be installed with command: 
#   
#    conda install -c salilab dssp
#
# @ Reference: 
#    1. Yan W, Hu G, Liang Z, et al. Node-Weighted Amino Acid Network Strategy for Characterization and Identification of Protein Functional Residues. J Chem Inf Model. 2018;58(9):2024-2032. 
#    2. Touw, W. G., Baakman, C., Black, J., te Beek, T. A., Krieger, E., Joosten, R. P., & Vriend, G. (2015). A series of PDB-related databanks for everyday needs. Nucleic acids research, 43(Database issue), D364–D368.
#   
# @ Python package in need: pandas
#
#############################################


import pandas as pd


def get_NACEN_betweenness(NACEN_output_dir):
    
    # Read NACEN.R programme output.
    NACEN = pd.read_csv(NACEN_output_dir, sep = "\t")
    
    # Add a column of Residues.
    NACEN["Residue"] = NACEN["ResName"].map(str) + NACEN["ResID"].map(str)
    #NACEN["Chain_Residue"] = NACEN["Chain"].map(str) + "_" + NACEN["ResName"].map(str) + NACEN["ResID"].map(str)
    
    # Get all chains in the progress.
    chains = set(NACEN["Chain"])
    
    # Summarize the betweenness and residues.
    NACEN_01 = NACEN[NACEN["Chain"] == sorted(list(chains))[0]]
    NACEN_02 = NACEN[NACEN["Chain"] == sorted(list(chains))[1]]
    
    sorted_NACEN_01 = NACEN_01.sort_values(by = "Betweenness")[::-1]
    protein_01_res = sorted_NACEN_01.iloc[:, 5]
    betweenness_01 = sorted_NACEN_01.iloc[:, 1]
    chain_01 = sorted_NACEN_01.iloc[:, 2]
    
    sorted_NACEN_02 = NACEN_02.sort_values(by = "Betweenness")[::-1]
    protein_02_res = sorted_NACEN_02.iloc[:, 5]
    betweenness_02 = sorted_NACEN_02.iloc[:, 1]
    chain_02 = sorted_NACEN_02.iloc[:, 2]
    
    # Summary of betweenness
    betweenness_top_df1 = pd.DataFrame()
    betweenness_top_df1["betweenness_rank"] = [i + 1 for i in range(len(betweenness_01))]
    betweenness_top_df1["chain_01"] = list(chain_01)
    betweenness_top_df1["betweenRes_01"] = list(protein_01_res)
    betweenness_top_df1["betweenness_01"] = list(betweenness_01)
    
    betweenness_top_df2 = pd.DataFrame()
    betweenness_top_df2["betweenness_rank"] = [i + 1 for i in range(len(betweenness_02))]
    betweenness_top_df2["chain_02"] = list(chain_02)
    betweenness_top_df2["betweenRes_02"] = list(protein_02_res)
    betweenness_top_df2["betweenness_02"] = list(betweenness_02)
    
    print("\n> *** Betweenness calculation finished! *** <\n")
    
    return (betweenness_top_df1, betweenness_top_df2)

