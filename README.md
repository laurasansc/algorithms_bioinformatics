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

