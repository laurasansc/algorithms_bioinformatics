# algorithms_bioinformatics

* read-in data
  - visualize the distribution of data to decide which threshold of IC50 to use
  - filter out non-binders (ic50 = 10,000 nM)
  - decide which alleles to use and split the peptide sequences(80/20) for train/test
  - cluster peptides with Hobohm and [heuristics](01_heuristics_LSC.ipynb) (1/r\*s)
  
* construct the PSSM for X MHC I alleles (X = 20)
* visualize the PSSMs as logos
* do [prediction](02_predictions_LSC.ipynb) with PSSM with test data
* [Pearson correlations](02_predictions_LSC.ipynb)



The code can be run all by downloading the entire repository and running the [\_\_main__.py](__main__.py) file. Be aware that in all the files there is a specified data directory which needs to be changed to the new directory where you downloaded the repository. PS: It needs to be fixed so you don't actually need to touch anything from the code and you can call several different variables. 
