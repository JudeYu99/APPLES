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
#
# @ Python package in need: scipy, pandas, numpy
#
#############################################


from scipy.special import comb
import pandas as pd
import numpy as np


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
    
    coevo_01 = pd.read_csv("../APPLES_OUT/Top_coevolve_residues_01.csv", sep = "\t")
    coevo_01_res = list(coevo_01["chain_01"] + "_" + coevo_01["coevoRes_01"])
    coevo_01_res_top = coevo_01_res[:round(len(coevo_01_res) * threshold)]
    
    coevo_02 = pd.read_csv("../APPLES_OUT/Top_coevolve_residues_02.csv", sep = "\t")
    coevo_02_res = list(coevo_02["chain_02"] + "_" + coevo_02["coevoRes_02"])
    coevo_02_res_top = coevo_02_res[:round(len(coevo_02_res) * threshold)]
    
    coevo_top = coevo_01_res_top + coevo_02_res_top
    
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
    
    PRS_top = coevo_top + PRS_01_res_top + PRS_02_res_top
    
    all_res_counts = len(NACEN_01_res) + len(NACEN_02_res) + len(PRS_01_res) + len(PRS_02_res)
    
    return (coevo_top, NACEN_top, PRS_top, all_res_counts)


# Calculate APPLES P values.
def APPLES_P(PPIs, Allo_top):
    
    coevo_top = Allo_top[0]
    NACEN_top = Allo_top[1]
    PRS_top = Allo_top[2]
    all_res_counts = Allo_top[3]
    
    temp_top = list(set(NACEN_top).union(set(PRS_top)))
    Allo_res_top = list(set(temp_top).union(set(coevo_top)))
    Allo_P = []
    
    for i in range(PPIs.shape[0]):
        temp_res = PPIs.iloc[i, 1]
        temp_counts = len(list(set(temp_res).intersection(set(Allo_res_top))))
        Allo_P.append(calculate_P(all_res_counts, temp_counts, PPIs.iloc[i, 2]))
        
    return Allo_P



def mean(input_list):
    list_len = len(input_list)
    list_sum = sum(input_list)
    list_mean = list_sum / list_len
    
    return list_mean
    

# Calculate APPLES Score.
def APPLES_Score(selection_01, selection_02, PPIs, coevo_result, PRS_result, NACEN_result):
    
    df1 = coevo_result[0]
    df1["chain_res"] = df1["chain_01"] + "_" + df1["coevoRes_01"]
    df2 = PRS_result[0]
    df2["chain_res"] = df2["chain_01"] + "_" + df2["effeRes_01"]
    df3 = NACEN_result[0]
    df3["chain_res"] = df3["chain_01"] + "_" + df3["betweenRes_01"]
    protein_01_sum = df1.merge(df2, on = "chain_res").merge(df3, on = "chain_res")[["chain_res", "coevolution_01", "effectiveness_01", "betweenness_01"]]
    protein_01_sum.columns = ['chain_res', 'coevolution', 'effectiveness', 'betweenness']
    
    df4 = coevo_result[1]
    df4["chain_res"] = df4["chain_02"] + "_" + df4["coevoRes_02"]
    df5 = PRS_result[1]
    df5["chain_res"] = df5["chain_02"] + "_" + df5["effeRes_02"]
    df6 = NACEN_result[1]
    df6["chain_res"] = df6["chain_02"] + "_" + df6["betweenRes_02"]
    protein_02_sum = df4.merge(df5, on = "chain_res").merge(df6, on = "chain_res")[["chain_res", "coevolution_02", "effectiveness_02", "betweenness_02"]]
    protein_02_sum.columns = ['chain_res', 'coevolution', 'effectiveness', 'betweenness']

    APPLES_sum = pd.concat([protein_01_sum, protein_02_sum], axis = 0)
    APPLES_sum.index = APPLES_sum["chain_res"]

    std_cols = ["coevolution", "effectiveness", "betweenness"]

    # Rescaling (max-min normalization)
    for item in std_cols:
        max_tmp = np.max(np.array(APPLES_sum[item]))
        min_tmp = np.min(np.array(APPLES_sum[item]))
        if (max_tmp != min_tmp):
            APPLES_sum[item] = APPLES_sum[item].apply(lambda x: (x - min_tmp) / (max_tmp - min_tmp))

    APPLES_Score = []

    for i in range(PPIs.shape[0]):
        
        coevo_temp = []
        PRS_temp = []
        NACEN_temp = []
        
        for temp_res in PPIs.iloc[i, 1]:
            coevo_temp.append(APPLES_sum.loc[temp_res][1])
            PRS_temp.append(APPLES_sum.loc[temp_res][2])
            NACEN_temp.append(APPLES_sum.loc[temp_res][3])
        
        APPLES_out = sum([mean(coevo_temp), mean(PRS_temp), mean(NACEN_temp)])

        APPLES_Score.append(format(APPLES_out,'.3f'))
        
    return APPLES_Score

  