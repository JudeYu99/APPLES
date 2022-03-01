#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@ author: Yu Zhu

@ Email: yzhu99@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

"""

### Introduction of get_Pocket_residues.py
#
# @ This part of program is delicated for getting pocket residues from Fpocket outputs.
# @ Reference: Le Guilloux V, Schmidtke P, Tuffery P. Fpocket: an open source platform for ligand pocket detection. BMC Bioinformatics. 2009;10:168. Published 2009 Jun 2.
#
# @ Python package in need: os, pandas
#
#############################################


import os
import pandas as pd
from Fpocket_calculation import *


# Get pocket residues.
def Pockets_Residues_all(path, selection_01, selection_02):
    
    os.chdir(path)
    
    if check_path(path):
        
        files = os.listdir()
        df = pd.DataFrame(columns = ["Fpocket_Number", "Chain_Residues"])
       
        all_res = []
        for file in files:
            if file.endswith(".pdb"):
                pocket_number = file.split("_")[0]
                temp_dic = find_pocket_residues_from_pdb(file)
                
                for chain in list(temp_dic.keys()):
                    for res in temp_dic[chain]:
                        chain_res_num = str(chain) + "_" + str(res)
                        all_res.append(chain_res_num)
                
                df.loc[pocket_number] = [pocket_number, temp_dic]

        df.sort_index().to_csv("../../APPLES_OUT/Fpockets_all_residues.csv", sep = "\t", index = False)

        all_res = list(set(all_res))
        
        os.chdir("../../codes")
        
        print("\n> *** Jobs Done! *** <\n")
    else:
        print("\n> *** Invalid input! *** <\n")

    return all_res


# Get PPI pocket residues.
def PPI_Pockets_Residues_all(selection_01, selection_02):
    
    PPIs = pd.read_csv("../APPLES_OUT/Fpockets_PPI_residues.csv", sep = "\t")
    PPI_all_res = []
    
    for i in range(PPIs.shape[0]):
        temp_dic = eval(PPIs["Chain_Residues"][i])
           
        for res in temp_dic[selection_01]:
            temp_res = selection_01 + "_" + str(res)
            PPI_all_res.append(temp_res)
    
        for res in temp_dic[selection_02]:
            temp_res = selection_02 + "_" + str(res)
            PPI_all_res.append(temp_res)
        
    PPI_all_res = list(set(PPI_all_res))
        
    return PPI_all_res


# Summary for each PPI pockets.
def sum_PPIs(selection_01, selection_02):
    
    PPIs = pd.read_csv("../APPLES_OUT/Fpockets_PPI_residues.csv", sep = "\t")
    res_num = []
    chain_res = []
    
    for i in range(PPIs.shape[0]):
        temp_dic = eval(PPIs["Chain_Residues"][i])
        pocket_res = []
        
        for res in temp_dic[selection_01]:
            temp_res = selection_01 + "_" + str(res)
            pocket_res.append(temp_res)
        
        for res in temp_dic[selection_02]:
            temp_res = selection_02 + "_" + str(res)
            pocket_res.append(temp_res)
        
        res_num.append(len(pocket_res))
        chain_res.append(pocket_res)
        
    PPIs["Chain_Res"] = chain_res
    PPIs["Residues_Counts"] = res_num
    PPIs.drop(labels = 'Chain_Residues', axis = 1, inplace = True)
    
    return PPIs

  