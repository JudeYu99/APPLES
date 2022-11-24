#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@ author: Yu Zhu

@ Email: yu.zhu.22@ucl.ac.uk

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

"""

### Introduction of MI_coevolution.py
#
# @ This part of program is delicated for calculating conservation and co-evolution score based on merged MSAs.
# @ Two separate MSA files are required when running this sub-program.
# @ A merged MSA file will be created in fasta format.
# @ Reference: Liu Y, Bahar I. Sequence Evolution Correlates with Structural Dynamics 2012 Mol Biol Evol 29(9):2253-2263; Liu Y, Gierasch LM, Bahar I Role of Hsp70 ATPase domain intrinsic dynamics and sequence evolution in enabling its functional interactions with NEFs 2010 PLoS Comput Biol 6(9)
#
# @ Python package in need: prody, numpy, pandas
#
#############################################


from prody import *
import numpy as np
import pandas as pd


def MI_coevolution(msa_file_01, msa_file_02, prefix = "../protein_protein"):

    # Load MSA files.
    msa_01 = parseMSA(msa_file_01)
    msa_02 = parseMSA(msa_file_02)

    # Merge MSA.
    msa = mergeMSA(msa_01, msa_02)
    
    print("\n> *** Merge MSAs finished! *** <\n")

    # Refine MSAs.
    msa_refine_01 = refineMSA(msa_01, label = "9606", rowocc = 0.4, seqid = 0.98)
    msa_refine_02 = refineMSA(msa_02, label = "9606", rowocc = 0.4, seqid = 0.98)
    msa_refine = refineMSA(msa, label = "9606", rowocc = 0.4, seqid = 0.98)

    # Get sequence length.
    protein_01_len = msa_refine_01.numResidues()
    protein_02_len = msa_refine_02.numResidues()
    protein_len = msa_refine.numResidues()

    # Write merged MSA files.
    writeMSA(prefix + '_msa_aligned.fasta', msa_refine)

    # Co-evolution calculation with mutual information.
    mutinfo = buildMutinfoMatrix(msa_refine)
    mutinfo_inter_protein = mutinfo[:protein_01_len, protein_01_len:]
    mutinfo_inter_df = pd.DataFrame(mutinfo_inter_protein)
    
    coevo_01 = list(mutinfo_inter_df.mean(1))
    coevo_02 = list(mutinfo_inter_df.mean(0))
    
    # Get the whole sequence of HUMAN proteins from merged MSA.
    index_9606 = msa_refine.getIndex("9606")
    seq_9606 = msa_refine[index_9606]
    seq_9606_Res = seq_9606.getArray()
    
    
    # Get the corresponding residue name and number (e.g. H11, M20) for two proteins.
    seq_01 = []
    for i in range(protein_01_len):
        tmp_res = AA_convert(str(seq_9606_Res[i])[2])
        Res = tmp_res + str(i + 1)
        seq_01.append(Res)
        
    seq_02 = []
    for i in range(protein_02_len):
        tmp_res = AA_convert(str(seq_9606_Res[i + protein_01_len])[2])
        Res = tmp_res + str(i + 1)
        seq_02.append(Res)

    coevo_df1 = pd.DataFrame()
    #coevo_df1["coevolution_index"] = [i + 1 for i in range(protein_01_len)]
    coevo_df1["coevoRes_01"] = seq_01
    coevo_df1["coevolution_01"] = coevo_01
    coevo_top_df1 = coevo_df1.sort_values(by = "coevolution_01", ascending = False)
    coevo_top_df1["coevo_rank"] = [i + 1 for i in range(protein_01_len)]
    coevo_top_df1.reindex(columns = ["coevo_rank", "coevoRes_01", "coevolution_01"])    

    coevo_df2 = pd.DataFrame()
    #coevo_df2["coevolution_index"] = [i + 1 for i in range(protein_02_len)]
    coevo_df2["coevoRes_02"] = seq_02
    coevo_df2["coevolution_02"] = coevo_02
    coevo_top_df2 = coevo_df2.sort_values(by = "coevolution_02", ascending = False)
    coevo_top_df2["coevo_rank"] = [i + 1 for i in range(protein_02_len)]
    coevo_top_df2.reindex(columns = ["coevo_rank", "coevoRes_02", "coevolution_02"])

    print("\n> *** MI coevolution finished! *** <\n")

    return (coevo_top_df1, coevo_top_df2)


def AA_convert(RES):
    # Input 3-letter AA code.
    AA = {"R":"ARG", "H":"HIS", "K":"LYS", "D":"ASP", "E":"GLU", "S":"SER", "T":"THR", "N":"ASN", "Q":"GLN", "C":"CYS", "G":"GLY", "P":"PRO", "A":"ALA", "V":"VAL", "I":"ILE", "L":"LEU", "M":"MET", "F":"PHE", "Y":"TYR", "W":"TRP"}

    return AA[RES]
  
