#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@ author: Yu Zhu

@ Email: yu.zhu.22@ucl.ac.uk

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

"""

### Introduction of pdb_Parsing.py
#
# @ This part of program is delicated for getting two sequences and corresponding structure given a protein complex pdb file.
#
# @ Python package in need: os, pandas
#
#############################################

def AA_dict(RES):
    # Input 3-letter AA code.
    AA = {"ARG":"R", "HIS":"H", "LYS":"K", "ASP":"D", "GLU":"E", "SER":"S", "THR":"T", "ASN":"N", "GLN":"Q", "CYS":"C", "GLY":"G", "PRO":"P", "ALA":"A", "VAL":"V", "ILE":"I", "LEU":"L", "MET":"M", "PHE":"F", "TYR":"Y", "TRP":"W"}

    return AA[RES]


def is_chain_valid(pdb_file, selection_01, selection_02):
    
    input_pdb = open(pdb_file, 'r+')
    chains = []
    
    for lines in input_pdb:
        if lines.startswith('ATOM'):
            chains.append(lines[21].strip())
    
    Chains = list(set(chains))
    flag = False
    
    if selection_01 in Chains and selection_02 in Chains and selection_01 != selection_02:
        flag = True
    
    return flag


# Get the AA sequence given chain identifiers.
def get_Sequence(pdb_file, selection_01, selection_02):
    
    input_pdb = open(pdb_file, 'r+')
    
    chain_01 = {}
    chain_02 = {}
    
    for lines in input_pdb:
        if lines.startswith('ATOM') and lines[21].strip() == selection_01:
            resnum = lines[22:26].strip()
            resname = AA_dict(lines[17:21].strip())
            chain_01[resnum] = resname
            
        elif lines.startswith('ATOM') and lines[21].strip() == selection_02:
            resnum = lines[22:26].strip()
            resname = AA_dict(lines[17:21].strip())
            chain_02[resnum] = resname

    seq_01 = "".join(chain_01.values())
    seq_02 = "".join(chain_02.values())
    
    return (seq_01, seq_02)


# Get the sub structure given chain identifiers.
def get_Sub_Structure(pdb_file, selection_01, selection_02):
    
    input_pdb = open(pdb_file, 'r+')
    new_pdb_file = []
    
    for lines in input_pdb:
        
        if lines.startswith('ATOM') and lines[21].strip() == selection_01:
            new_pdb_file.append(lines)
            
        elif lines.startswith('ATOM') and lines[21].strip() == selection_02:
            new_pdb_file.append(lines)

    pdb_identifier = pdb_file.split("/")[-1][:-4]
    new_pdb_filename = "../" + pdb_identifier + "_" + selection_01 + "_" + selection_02 + ".pdb"
    
    with open(new_pdb_filename, "w+") as output:
        output.writelines(new_pdb_file)
        
    return new_pdb_filename
  
