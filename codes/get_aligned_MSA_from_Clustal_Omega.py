#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@ author: Yu Zhu

@ Email: yzhu99@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

"""

### Introduction of get_aligned_from_Clustal_Omega.py
#
# @ This part of program is delicated for getting aligned MSAs of protein dimers from Clustal Omega Web Server (https://www.ebi.ac.uk/Tools/msa/clustalo/).
# @ Network connection is required when running this sub-program.
# @ Reference: Madeira F, Park YM, Lee J, et al. The EMBL-EBI search and sequence analysis tools APIs in 2019. Nucleic Acids Research. 2019 Jul;47(W1):W636-W641. 
#
# @ Python package in need: selenium, wget, time
# @ BLAST parameters: Enter or paste a set of - PROTEIN
#                     OUTPUT FORMAT - ClustalW
#                     Other parameters - Default
#
#############################################


from selenium import webdriver
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.common.by import By
import time
import wget


# Get Clustal MSA result url.
def get_aligned_MSA_from_Clustal_Omega(merged_sequence):
    
    # Running the programme silently.
    option = webdriver.ChromeOptions()
    option.add_argument('headless')

    # Create a new Chrome session.
    driver = webdriver.Chrome(chrome_options=option)
    
    # Navigate to the UniProt BLAST page.
    driver.get("https://www.ebi.ac.uk/Tools/msa/clustalo/")
    
    # Hide the banner.
    driver.find_element(By.ID, "data-protection-agree").click()
    
    # Get the search textbox and fill it.
    search_field = driver.find_element(By.ID, "sequence")
    search_field.send_keys(merged_sequence)
    
    # Set OUTPUT FORMAT as Pearson/FASTA.
    driver.find_element(By.XPATH, '//*[@id="outfmt"]/option[3]').click()
    
    time.sleep(3)
    
    # Submit the job.
    driver.find_element(By.XPATH, '//*[@id="jd_submitButtonPanel"]/input').click()
    
    # Wait until the Download button shows up (job finished).
    # Time limitation: 3600s (1 hour).
    WebDriverWait(driver, 3600).until(EC.element_to_be_clickable((By.XPATH, '//*[@id="alnFile"]'))).click()
    
    Clustal_MSA_url = driver.current_url
    
    driver.quit()
    
    return Clustal_MSA_url


# Download the MSA file in fasta format.
def get_Clustal_MSA_in_fasta(Clustal_MSA_url, filename):

    wget.download(Clustal_MSA_url, out = filename, bar = None)


# Running Part
def Clustal_Omega(fasta_file):
    
    with open(fasta_file, "r+") as MSAs:
        merged_sequence = MSAs.read()
	
    Clustal_MSA_url = get_aligned_MSA_from_Clustal_Omega(merged_sequence)
    
    out_prefix = fasta_file.split("blast")[0]
    
    get_Clustal_MSA_in_fasta(Clustal_MSA_url, out_prefix + "aligned_blast.fasta")

    print("\n> *** Clustal Alignment " + str(fasta_file[14:16]) + " finished! *** <\n")
  