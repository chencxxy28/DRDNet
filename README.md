# DRDNet (Disease risk-associated pseudo-dynamic networks)
A statistical framework for recovering dynamic networks from static data 

# simulations
In the simulations folder, we provide the codes to generate the results for two cases in the main text of the paper: 1. when the index(agent) is known; 2 when the index(agent) is unobserved and imputed(estimated). For each case, we consider four combinations of sample sizes (n=100, n=200) and variances of measurement errors (0.1 and 0.5). They are named: n100v1_one, n100v5_one, n200v1_one, n200v5_one, respectively. By running these codes can reproduce the result in the manuscript.  

# real data applications
In the real data application folder, we provide three code files: the first one should run is called "dataload_clean_geneselect_seperate", which processes the data, calculates the imputed index(agent), selects important genes; after that, run the code named "change ENS to gene names" to change ENS label to gene names. Finally, run the code named "cvm_gode_07152021"  to conduct DRDNet and make network plots. The user should install necessary packages and change the directory in r script. The datasets used for the analyses described in this manuscript were obtained from dbGaP at http://www.ncbi.nlm.nih.gov/gap through dbGaP accession number phs000424.v7.p2.

