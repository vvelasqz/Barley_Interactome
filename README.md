# Support code for Velásquez-Zapata *et.al.* (2021) An interolog-based barley interactome as an integration framework for immune signaling.   

This repository contains the support code for the manuscript Velásquez-Zapata *et.al.* (2021) An interolog-based barley interactome as an integration framework for immune signaling.   

**Overview**

This repository contains the main functions that were implemented to develop the data analysis of the manuscript. It is separated in several R functions that were used across the different analyses through the results section. 

The paper "An interolog-based barley interactome as an integration framework for immune signaling " has been submitted for peer review. Please refer to the manuscript for more details about each function. The entire pipeline has been coded in R.

**Citations**

If you use part of this code please cite  

* Velásquez-Zapata V, Elmore JM, Fuerst G, Wise RP. 2021. An interolog-based barley interactome as an integration framework for immune signaling.


**Software requirements**

* R
* R packages: reshape2, tidyverse, psych, DESeq2, mass, igraph

**Functions in this repository**

Seven R files were created with the main functions used to perform the analysis. A brief description of each file is provided: 

1. Comparison HvInt and AtInt.R contains the network comparison between the predicted barley interactome (HvInt) and the experimental Arabidospsis interactome (AtInt).

2. Node stats DiffSLC.R contains a function to measure all the topological properties of the nodes used in the manuscript. It also contains our implementation of DiffSLC with RNASeq data. 

3. GO analysis.R contains the code to perform GO analysis.

4. Hypergeometric tests.R contains several functions to do the hypergeometric tests to measure enrichment across the manuscript.

5. Permute correlation.R contains the functions to measure distance correlation among pairs of genes and the permutation tests to build the resistant and susceptible interactomes.

6. Co-expressed interactome analysis.R contains the network analyses of the co-expressed interactomes and their comparisons.

7. MLAInt construction.R contains the code for building the MLAInt network and its analysis.

 
