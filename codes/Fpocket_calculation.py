#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@ author: Yu Zhu

@ Email: yu.zhu.22@ucl.ac.uk

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

"""

### Introduction of Fpocket_calculation.py
#
# @ This part of program is delicated for calculate surface cavity by Fpocket programme and detect pockets on the interface given a protein complex pdb file.
# @ Fpocket can be downloaded from http://fpocket.sourceforge.net or https://github.com/Discngine/fpocket.
# @ Note: Fpocket must be installed in advance.
# @ Reference: Le Guilloux V, Schmidtke P, Tuffery P. Fpocket: an open source platform for ligand pocket detection. BMC Bioinformatics. 2009;10:168.
#
# @ Python package in need: os, pandas
#
#############################################


import os
import pandas as pd


# Perform Fpocket calculation. 
def Fpocket_calculation(pdb_file):
    
    if pdb_file.endswith(".pdb"):
        os.system(("fpocket -f {}").format(pdb_file))
    else:
        print("Invalid input!")


# Check the path of Fpocket outputs folder.
def check_path(path):
    
    flag = False
    
    if path.endswith("_out") or path.endswith("_out/"):
        os.chdir("./pockets")
        flag = True
        print("\n> *** Path has been changed to pockets folder. *** <")
    
    elif path.endswith("pockets") or path.endswith("pockets/"):
        flag = True
        print("\n> *** Path has already been in pockets folder. *** <")
        
    else:
        print("\n> *** Please check the path! *** <")

    return flag


# Get pocket residues from one pdb file.
# Return a dictionary whose keys are chain identifiers and values are residues. 
def find_pocket_residues_from_pdb(pdb_file):
    
    with open(pdb_file, 'r+') as input_pdb:
        lines = input_pdb.readlines()
    
    residues = []
    chains = []
    
    for line in lines:
        if line.startswith('ATOM'):
            resname = line[17:20].strip()
            resnum = line[23:26].strip()
            chain = line[21].strip()
            chains.append(chain)
            residue = chain + "_" + resname + str(resnum)
            residues.append(residue)
    
    chains = list(set(chains))
    
    chain_residue = {}
    
    for chain in chains:
        chain_residue[chain] = []
    
        for residue in residues:
            if residue.startswith(chain):
                chain_residue[chain].append(residue[2:])
                
        chain_residue[chain] = list(set(chain_residue[chain]))
        
    return chain_residue


# Detect whether exists PPI pockets.
def PPI_pockets_detection(selection_01, selection_02):
    path = os.path.abspath(".")
    PPI_flag = []
    
    if check_path(path):
        
        files = os.listdir()
        df = pd.DataFrame(columns = ["Fpocket_Number", "Chain_Residues"])
            
        for file in files:
            if file.endswith(".pdb"):
                pocket_number = file.split("_")[0]
                temp_dic = find_pocket_residues_from_pdb(file)
                
                if len(temp_dic) >= 2 and selection_01 in temp_dic.keys() and selection_02 in temp_dic.keys():
                    df.loc[pocket_number] = [pocket_number, temp_dic]
                    df.sort_index().to_csv("../../APPLES_OUT/Fpockets_PPI_residues.csv", sep = "\t", index = False)
                    PPI_flag.append("Detected")
                
                else:
                    
                    PPI_flag.append("Undetected")
                    
        if "Detected" in PPI_flag:
            print("\n> *** PPI pockets detected! *** <")
        else:
            print("\n> *** No PPI pockets detected! *** <")
        
        print("\n> *** Jobs Done! *** <\n")
    
    else:
        
        print("\n> *** Invalid input! *** <\n")
    
    return PPI_flag
  
