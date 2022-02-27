# APPLES
**A**llosteric **P**rotein-protein interaction **P**ocket **L**ikelihood **E**valuation **S**ystem (APPLES)

## 1.  Introduction

Allosteric effect, as the "second code of life", is a direct and effective way for protein macromolecular machines to regulate their own functions, and will become a research hotspot in the post-AlphaFold era. Based on our previous research results, this project intends to focus on the allosteric mechanism of protein complexes  at the level of structural systems biology, and comprehensively use a variety of theoretical calculation methods to study the allosteric mechanism for targeting PPIs. Integrating the sequence conservation and coevolutionary information, as well as the topological and dynamic properties of network, we built **APPLES** - a computing platform for the systematic study of allosteric mechanism. We hoped that the completion of this project can provide theoretical guidance and bioinformatics tools for the design of allosteric drugs targeting the protein-protein interactions.


## 2. Usage
  
  1) Anaconda is recommended for the APPLES and only Unix/Linux/MacOS operating systems are allowed.
  2) Python packages in need: os, re, string, sys, time, numpy, scipy, pandas, prody, wget, selenium
  3) R packages in need: NACEN
  4) Python language version: >= 3.7
  5) R language version: >= 3.2.0
  6) NACEN package can be downloaded in the website: http://sysbio.suda.edu.cn/NACEN/NACEN.rar. NACEN requires DSSP programme for building AAN which can be installed with command: 
  ```
  conda install -c salilab dssp
  ```  
  7) The input of the program is a PDB file containing multiple chains and users can choose the chains they need. Before running ***Run_APPLES.py***, please modify line 35, line 37 for PDB chain identifiers and line 38 for PDB files. Then, change the working directory into codes folder and run the following command.
  ```
  python3 Run_APPLES.py
  ```

## 3. Examples

- Example PDB input file can be referred to [***PPI.pdb***](https://github.com/JudeYu99/APPLES/blob/main/PPI.pdb).



## 4. Reference
- Anaconda Software Distribution. (2020). Anaconda Documentation. Anaconda Inc.
- Yan W, Hu G, Liang Z, et al. Node-Weighted Amino Acid Network Strategy for Characterization and Identification of Protein Functional Residues. J Chem Inf Model. 2018;58(9):2024-2032.
- Touw, W. G., Baakman, C., Black, J., te Beek, T. A., Krieger, E., Joosten, R. P., & Vriend, G. (2015). A series of PDB-related databanks for everyday needs. Nucleic acids research, 43(Database issue), D364–D368.
- Le Guilloux V, Schmidtke P, Tuffery P. Fpocket: an open source platform for ligand pocket detection. BMC Bioinformatics. 2009;10:168. Published 2009 Jun 2.
- Liu, Y., & Bahar, I. (2012). Sequence evolution correlates with structural dynamics. Molecular biology and evolution, 29(9), 2253–2263.
- Liu, Y., Gierasch, L. M., & Bahar, I. (2010). Role of Hsp70 ATPase domain intrinsic dynamics and sequence evolution in enabling its functional interactions with NEFs. PLoS computational biology, 6(9), e1000931.
- Harris, C. R., Millman, K. J., van der Walt, S. J., Gommers, R., Virtanen, P., Cournapeau, D., … Oliphant, T. E. (2020). Array programming with NumPy. Nature, 585, 357–362.
- Virtanen, P., Gommers, R., Oliphant, T. E., Haberland, M., Reddy, T., Cournapeau, D., … SciPy 1.0 Contributors. (2020). SciPy 1.0: Fundamental Algorithms for Scientific Computing in Python. Nature Methods, 17, 261–272. 
- McKinney, W., & others. (2010). Data structures for statistical computing in python. In Proceedings of the 9th Python in Science Conference (Vol. 445, pp. 51–56).


