# Instructions for running the hotspot pipeline

# Step-by-step instructions
- Run first steps of 05_hotspotsPred.Rmd
- Move output from output/data/hotspots/domains/domainsPreAlign and output/data/hotspots/ubnSites_forHotspots to a computational cluster
- Run scripts/hotspot_scripts/02_runMAFFT_for_hotspots.py to align domain sequences (conda environment: mafft)
- Run scripts/hotspot_scripts/03_run_ptm_hotspots.py to run hotspot analysis (conda environment: ptm_hotspots)
- Move output from output/data/hotspots/hotspotOutput.csv back to local. Run more steps of 05_hotspotPred.Rmd
- Move output from output/data/hotspots/domains/representativeStrs/strRep_seqDomain.fasta to the cluster 
- Run scripts/hotspot_scripts/04_repStrsMAFFT.py to align domain sequences with the representative structure (conda environment: mafft)
- Run scripts/hotspot_scripts/05_repStrsHotspots.py to run hotspot analysis with the representative structure (conda environment: ptm_hotspots)
- Move output from output/data/hotspots/domains/representativeStrs back to local 

# Conda environments
## ptm_hotspots
conda create --name ptm_hotspots python=3 numpy pandas scipy statsmodels
conda install -c conda-forge biopython
## mafft
conda create --name mafft
conda install conda-forge::mafft
conda install -c conda-forge biopython