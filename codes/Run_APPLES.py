#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@ author: Yu Zhu

@ Email: yzhu99@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

"""


### Introduction of Run_APPLES.py
#
# @ This part of program is delicated for running the APPLES program.
#
# @ Python package in need: time, os
#
#############################################


import time
import os
from pdb_Parsing import *
from Fpocket_calculation import *
from get_BLAST_from_UniProt import *
from merge_MSAs import *
from get_aligned_MSA_from_Clustal_Omega import *
from MI_coevolution import *
from PRS_effectiveness import *
from NACEN_R import *
from get_Pocket_residues import *
from P_Value import *


# Parameters for modifying.
selection_01 = "A"
selection_02 = "B"
pdb_file = "../PPI.pdb"


if is_chain_valid(pdb_file, selection_01, selection_02):
    
    print("\n========== ****** @ Start Computing. ****** ==========\n")
    
    os.mkdir("../APPLES_OUT")
    sequences = get_Sequence(pdb_file, selection_01, selection_02)
    new_pdb_name = get_Sub_Structure(pdb_file, selection_01, selection_02)
    seq_01 = sequences[0]
    seq_02 = sequences[1]
    time.sleep(3)

    print("\n========== ****** @ Running Pockets Detecting. ****** ==========\n")
    Fpocket_in = pdb_file.split(".pdb")[0] + "_" + selection_01 + "_" + selection_02
    Fpocket_calculation(Fpocket_in + ".pdb")
    
    os.chdir(Fpocket_in + "_out")
    PPI_flag = PPI_pockets_detection(selection_01, selection_02)
    time.sleep(3)

    os.chdir("../../codes")
    
    if "Detected" in PPI_flag:

        print("\n========== ****** @ Running BLAST. ****** ==========\n")
        BLAST_01 = BLAST_MSA(seq_01, "01")
        BLAST_02 = BLAST_MSA(seq_02, "02")
        
        fasta_01 = "../APPLES_OUT/01_blast.fasta"  
        fasta_02 = "../APPLES_OUT/02_blast.fasta"
        write_aligned_fasta(fasta_01, fasta_02)
        
        
        print("\n========== ****** @ Running Clustal Omega. ****** ==========\n")
        Clustal_Omega("../APPLES_OUT/01_02_blast.fasta")
        Clustal_Omega("../APPLES_OUT/02_01_blast.fasta")
        time.sleep(3)
    
    
        print("\n========== ****** @ Running Co-evolution Calculation. ****** ==========\n")
        msa_file_01 = "../APPLES_OUT/01_02_aligned_blast.fasta"
        msa_file_02 = "../APPLES_OUT/02_01_aligned_blast.fasta"
        
        time.sleep(1)
        coevo_result = MI_coevolution(msa_file_01, msa_file_02, prefix = "../APPLES_OUT/protein_protein")
        coevo_01_df = coevo_result[0]
        coevo_01_df["chain_01"] = selection_01
        coevo_01_df = coevo_01_df[["coevo_rank", "chain_01", "coevoRes_01", "coevolution_01"]]
        coevo_01_df.to_csv("../APPLES_OUT/Top_coevolve_residues_01.csv", index = False, sep = "\t")
    
        coevo_02_df = coevo_result[1]
        coevo_02_df["chain_02"] = selection_02
        coevo_02_df = coevo_02_df[["coevo_rank", "chain_02", "coevoRes_02", "coevolution_02"]]
        coevo_02_df.to_csv("../APPLES_OUT/Top_coevolve_residues_02.csv", index = False, sep = "\t")
        time.sleep(3)
        

        print("\n========== ****** @ Running PRS Effectiveness Calculation. ****** ==========\n")
        PRS_result = PRS_effectiveness(pdb_file, selection_01, selection_02)
        PRS_result[0].to_csv("../APPLES_OUT/Top_PRS_effectiveness_residues_01.csv", index = False, sep = "\t")
        PRS_result[1].to_csv("../APPLES_OUT/Top_PRS_effectiveness_residues_02.csv", index = False, sep = "\t")
        time.sleep(3)
        
    
        print("\n========== ****** @ Running NACEN Betweenness Calculation. ****** ==========\n")
        os.system("Rscript NACEN.R -d /opt/anaconda3/bin/mkdssp -p " + Fpocket_in + ".pdb")
        NACEN_output_dir = "../APPLES_OUT/NACEN_output.csv"
        NACEN_result = get_NACEN_betweenness(NACEN_output_dir)
        NACEN_result[0].to_csv("../APPLES_OUT/Top_Betweenness_residues_01.csv", index = False, sep = "\t")
        NACEN_result[1].to_csv("../APPLES_OUT/Top_Betweenness_residues_02.csv", index = False, sep = "\t")
        time.sleep(3)
    
    
        print("\n========== ****** @ Running Pocket Residues Analysis. ****** ==========\n")
        path = Fpocket_in + "_out"

        pocket_res = Pockets_Residues_all(path, selection_01, selection_02)
        PPI_pocket_res = PPI_Pockets_Residues_all(selection_01, selection_02)
        PPIs = sum_PPIs(selection_01, selection_02)
        PPIs_P = []
        for i in range(PPIs.shape[0]):
            temp_p = calculate_P(len(pocket_res), PPIs.iloc[i, 2], len(PPI_pocket_res))
            PPIs_P.append(temp_p)
        PPIs["PPI_P_Value"] = PPIs_P
        print("\n> *** PPI P values calculated! *** <\n")
        
        Allo_top = Allo_res_top(threshold = 0.5)
        APPLES_P = APPLES_P(PPIs, Allo_top)
        PPIs["APPLES_P_Value"] = APPLES_P
        print("\n> *** APPLES P values calculated! *** <\n")
                
        APPLES_Score = APPLES_Score(selection_01, selection_02, PPIs, coevo_result, PRS_result, NACEN_result)
        PPIs["APPLES_Score"] = APPLES_Score
        print("\n> *** APPLES score calculated! *** <\n")
        PPIs.to_csv("../APPLES_OUT/APPLES_PPI.csv", index = False, sep = "\t")
        
        print("\n========== ****** @ APPLES Calculation finished! ****** ==========\n")

else:
    
    print("\n> *** Invalid input! *** <\n")
  