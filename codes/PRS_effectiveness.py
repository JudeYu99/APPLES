#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@ author: Yu Zhu

@ Email: yzhu99@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

"""

### Introduction of PRS_effectiveness.py
#
# @ This part of program is delicated for calculating PRS effectiveness based on merged MSAs.
# @ Two separate MSA files are required when running this sub-program.
# @ Reference: Atilgan C, Atilgan AR. Perturbation-response scanning reveals ligand entry-exit mechanisms of ferric binding protein. PLoS Comput Biol 2009 5(10):e1000544.
#
# @ Python package in need: prody, numpy, pandas
#
#############################################


from prody import *
from pylab import *
import pandas as pd
import time


def PRS_effectiveness(pdb, selection_01, selection_02):
    
    time.sleep(1)
    # Parse pdb structure.
    protein_complex = parsePDB(pdb)
    
    selection = "calpha " + "and chain " + selection_01 + " " + selection_02

    protein_complex_ca = protein_complex.select(selection)
    protein_01 = protein_complex.select("calpha and chain " + selection_01)
    protein_02 = protein_complex.select("calpha and chain " + selection_02)
    
    chain_01_len = len(protein_01.getResindices())
    
    # ANM calculation.
    anm_complex = ANM('Protein Complex')
    anm_complex.buildHessian(protein_complex_ca)
    anm_complex.calcModes('all')
    
    # PRS calculation.
    prs_mat, effectiveness, sensitivity = calcPerturbResponse(anm_complex)
    
    # Effectiveness of protein 1 / chain 1.
    protein_01_Resnums = protein_01.getResnums()
    protein_01_Resnames = protein_01.getResnames()
    
    protein_01_res = []
    
    for i in range(len(protein_01_Resnums)):
        temp = protein_01_Resnames[i] + str(protein_01_Resnums[i])
        protein_01_res.append(temp)
        
    effective_01_index = effectiveness[:chain_01_len].argsort()[::-1]
    protein_01_effective = effectiveness[effective_01_index]
    
    # Effectiveness of protein 2 / chain 2.
    protein_02_Resnums = protein_02.getResnums()
    protein_02_Resnames = protein_02.getResnames()
    
    protein_02_res = []
    
    for i in range(len(protein_02_Resnums)):
        temp = protein_02_Resnames[i] + str(protein_02_Resnums[i])
        protein_02_res.append(temp)
    
    effective_02_index = (effectiveness[chain_01_len:].argsort() + chain_01_len)[::-1]
    protein_02_effective = effectiveness[effective_02_index]
    
    
    # Summary of conservation
    effective_top_df1 = pd.DataFrame()
    effective_top_df1["effectiveness_rank"] = [i + 1 for i in range(len(protein_01_Resnums))]
    effective_top_df1["chain_01"] = [selection_01] * len(protein_01_Resnums)
    effective_top_df1["effeRes_01"] = np.array(protein_01_res)[effective_01_index]
    effective_top_df1["effectiveness_01"] = protein_01_effective
    
    effective_top_df2 = pd.DataFrame()
    effective_top_df2["effectiveness_rank"] = [i + 1 for i in range(len(protein_02_Resnums))]
    effective_top_df2["chain_02"] = [selection_02] * len(protein_02_Resnums)
    effective_top_df2["effeRes_02"] = np.array(protein_02_res)[effective_02_index - chain_01_len]
    effective_top_df2["effectiveness_02"] = protein_02_effective
    
    print("\n> *** Effectiveness calculation finished! *** <\n")
    
    # Map PRS effectiveness to pdb structures.
    #anm_complex.setBetas(effectiveness)
    #writePDB('proteins_ca_effectiveness.pdb', anm_complex)
    
    return (effective_top_df1, effective_top_df2)
  