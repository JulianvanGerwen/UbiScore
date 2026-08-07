This is code that extracts compara protein sequence alignments and combines them with multi-species ubi-omics to estimate ubiquitin site conservation. 
It is originally based off the ptmAge code (https://github.com/evocellnet/ptmAge.git).

----- Guide to files ----

01_prepare_sites.Rmd - Prepare ubiquitin sites for conservation analysis
02_conservation_pipeline.ipynb - Extract compara sequences and associated ubi-site statistics
03_conservationClassification.Rmd - Use compara alignments and the associated ubi-site statistics to classify ubi-site conservation