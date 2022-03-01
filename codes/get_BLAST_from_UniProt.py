#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""

@ author: Yu Zhu

@ Email: yzhu99@stu.suda.edu.cn

@ Address: Center for Systems Biology, Department of Bioinformatics, School of Biology and Basic Medical Sciences, Soochow University, Suzhou 215123, China.

"""

### Introduction of get_BLAST_from_UniProt.py
#
# @ This part of program is delicated for getting MSAs of protein monomer from UniProt BLAST (https://www.uniprot.org/blast/).
# @ Network connection is required when running this sub-program.
# @ Reference: Stephen F. Altschul, Thomas L. Madden, Alejandro A.Schaffer, Jinghui Zhang, Zheng Zhang, Webb Miller, and David J.Lipman (1997), "Gapped BLAST and PSI-BLAST: a new generation of protein database search programs", Nucleic Acids Res. 25:3389-3402.
#
# @ Python package in need: selenium, wget, time
# @ BLAST parameters: Target database - UniProtKB reference proteomes plus Swiss-prot
#                     E-Threshold - 10
#                     Matrix - Auto
#                     Filtering - None
#                     Gapped - yes
#                     Hits - 1000
#############################################


from selenium import webdriver
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.wait import WebDriverWait
from selenium.webdriver.common.by import By
import time
import wget


# Get MSA file url and top identity identifier.
def get_BLAST_from_UniProt(sequence):
    
    # Running the programme silently.
    option = webdriver.ChromeOptions()
    option.add_argument('headless')

    # Create a new Chrome session.
    driver = webdriver.Chrome(chrome_options=option)

    # Navigate to the UniProt BLAST page.
    driver.get("https://www.uniprot.org/blast/")
 
    # Get the search textbox and fill it.
    search_field = driver.find_element(By.ID, "blastQuery")
    search_field.send_keys(sequence)

    # Set BLAST hits to 1000.
    driver.find_element(By.XPATH, '//*[@id="blast-options"]/div[6]/p[2]/select/option[6]').click()

    # Hide the banner and submit BLAST.
    driver.find_element(By.ID, "privacy-panel-accept").click()
    driver.find_element(By.ID, "blast-submit").click()

    # Wait until the Download button shows up (job finished).
    # Time limitation: 3600s (1 hour).
    WebDriverWait(driver, 3600).until(EC.element_to_be_clickable((By.XPATH, '//*[@id="download-button"]'))).click()

    # Download all 999 sequences.
    driver.find_element(By.XPATH, '//*[@id="all"]').click()

    # Download in fasta format.
    driver.find_element(By.XPATH, '//*[@id="format"]/option[1]').click()

    # Download in uncompressed state.
    driver.find_element(By.XPATH, '//*[@id="option-uncompressed"]').click()

    # Execute downloads.
    driver.find_element(By.XPATH, '//*[@id="menu-go"]').click()
    
    time.sleep(3)
    
    # Record the MSA url.
    MSA_url = driver.current_url

    # Close the browser window.
    driver.quit()

    return (MSA_url)


# Download the MSA file in fasta format.
def get_MSA_in_fasta(MSA_url, seq_order):
    
    filename = "../APPLES_OUT/" + str(seq_order) + "_blast.fasta"
    
    if MSA_url.endswith(".fasta"):
        wget.download(MSA_url, out = filename, bar = None)
    else:
        wget.download(MSA_url + ".fasta", out = filename, bar = None)
        

# Main function for running the programme.
def BLAST_MSA(sequence, seq_order):
    
    MSA_url = get_BLAST_from_UniProt(sequence)
    get_MSA_in_fasta(MSA_url, seq_order)
    
    print("\n> *** BLAST " + str(seq_order) + " finished! *** <\n")
    
    return (MSA_url, seq_order)
  
