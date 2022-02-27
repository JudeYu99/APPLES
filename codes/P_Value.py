#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@ author: Yu Zhu

@ Email: yzhu99@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

"""

### Introduction of P_Value.py
#
# @ This part of program is delicated for calculating P values from all output files.
# @ Reference: Le Guilloux V, Schmidtke P, Tuffery P. Fpocket: an open source platform for ligand pocket detection. BMC Bioinformatics. 2009;10:168. Published 2009 Jun 2.
#
# @ Python package in need: scipy, pandas
#
#############################################


from scipy.special import comb
import pandas as pd


# Calculate P values.
def calculate_P(N, M, n):
    temp = []
    for i in range(M):
        temp.append(comb(M, i) * comb(N-M, n-i) / comb(N, n))
        
    P_Value = 1 - sum(temp)
    
    if P_Value <= 0:
        P_Value = 0
    
    return (format(P_Value,'.2E'))


# Find out residues with top performance for APPLES P values.
def Allo_res_top(threshold):
    
    NACEN_01 = pd.read_csv("../APPLES_OUT/Top_Betweenness_residues_01.csv", sep = "\t")
    NACEN_01_res = list(NACEN_01["chain_01"] + "_" + NACEN_01["betweenRes_01"])
    NACEN_01_res_top = NACEN_01_res[:round(len(NACEN_01_res) * threshold)]
    
    NACEN_02 = pd.read_csv("../APPLES_OUT/Top_Betweenness_residues_02.csv", sep = "\t")
    NACEN_02_res = list(NACEN_02["chain_02"] + "_" + NACEN_02["betweenRes_02"])
    NACEN_02_res_top = NACEN_02_res[:round(len(NACEN_02_res) * threshold)]
    
    NACEN_top = NACEN_01_res_top + NACEN_02_res_top
    
    PRS_01 = pd.read_csv("../APPLES_OUT/Top_PRS_effectiveness_residues_01.csv", sep = "\t")
    PRS_01_res = list(PRS_01["chain_01"] + "_" + PRS_01["effeRes_01"])
    PRS_01_res_top = PRS_01_res[:round(len(PRS_01_res) * threshold)]
    
    PRS_02 = pd.read_csv("../APPLES_OUT/Top_PRS_effectiveness_residues_02.csv", sep = "\t")
    PRS_02_res = list(PRS_02["chain_02"] + "_" + PRS_02["effeRes_02"])
    PRS_02_res_top = PRS_02_res[:round(len(PRS_02_res) * threshold)]
    
    PRS_top = PRS_01_res_top + PRS_02_res_top
    
    all_res_counts = len(NACEN_01_res) + len(NACEN_02_res) + len(PRS_01_res) + len(PRS_02_res)
    
    return (NACEN_top, PRS_top, all_res_counts)


# Calculate APPLES P values.
def APPLES_P(PPIs, Allo_top):
    
    NACEN_top = Allo_top[0]
    PRS_top = Allo_top[1]
    all_res_counts = Allo_top[2]
    
    Allo_res_top = list(set(NACEN_top).union(set(PRS_top)))
    Allo_P = []
    
    for i in range(PPIs.shape[0]):
        temp_res = PPIs.iloc[i, 1]
        temp_counts = len(list(set(temp_res).intersection(set(Allo_res_top))))
        Allo_P.append(calculate_P(all_res_counts, temp_counts, PPIs.iloc[i, 2]))
        
    return Allo_P

