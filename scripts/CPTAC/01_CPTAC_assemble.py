# Here I fetch CPTAC data

###### Initialise #####
import cptac
import pandas as pd
import numpy as np
import re

data_dict = {}
transcriptomics_dict = {}
proteomics_dict = {}
phosphoproteomics_dict = {}
ubi_dict = {}
clin_dict = {}
cancer_types = ['hnscc','luad','ov','gbm','pdac',
                'lscc','coad','ucec','brca','ccrcc']

def subtract_by_median(slice):
    median = slice.median()
    return slice - median

##### Fetch data #####
# Access to the CPTAC data 
data_dict['hnscc'] = cptac.Hnscc()
data_dict['luad'] = cptac.Luad()
data_dict['ov'] = cptac.Ov()
data_dict['gbm'] = cptac.Gbm()
data_dict['pdac'] = cptac.Pdac()
data_dict['lscc'] = cptac.Lscc()
data_dict['coad'] = cptac.Coad()
data_dict['ucec'] = cptac.Ucec()
data_dict['brca'] = cptac.Brca()
data_dict['ccrcc'] = cptac.Ccrcc()

for ctype in cancer_types:
    transcriptomics_dict[ctype] = data_dict[ctype].get_dataframe('transcriptomics', 'bcm')
    proteomics_dict[ctype] = data_dict[ctype].get_dataframe('proteomics', 'bcm')
    phosphoproteomics_dict[ctype] = data_dict[ctype].get_dataframe('phosphoproteomics', 'bcm')
    clin_dict[ctype] = data_dict[ctype].get_clinical('mssm')
    
    # Aggregate the columns pSite-wise
    phosphoproteomics_dict[ctype].columns.names = ['Name','pSite','Peptide','ENSG']
    phosphoproteomics_dict[ctype] = phosphoproteomics_dict[ctype].T.groupby(level=['Name','pSite','ENSG']).median().T    
    
    # drop Non-tumor samples (.N)
    transcriptomics_dict[ctype] = transcriptomics_dict[ctype][~transcriptomics_dict[ctype].index.str.endswith('.N')]
    proteomics_dict[ctype] = proteomics_dict[ctype][~proteomics_dict[ctype].index.str.endswith('.N')]
    phosphoproteomics_dict[ctype] = phosphoproteomics_dict[ctype][~phosphoproteomics_dict[ctype].index.str.endswith('.N')]
    clin_dict[ctype] = clin_dict[ctype][~clin_dict[ctype].index.str.endswith('.N')]
    
     # Calculate log2 ratio by subtracting median of log2 values (gene/protein wise)
    transcriptomics_dict[ctype] = transcriptomics_dict[ctype].apply(subtract_by_median, axis=0)
    proteomics_dict[ctype] = proteomics_dict[ctype].apply(subtract_by_median, axis=0)
    phosphoproteomics_dict[ctype] = phosphoproteomics_dict[ctype].apply(subtract_by_median, axis=0)

# create a merged file
transcriptomics = pd.concat(transcriptomics_dict.values())
proteomics = pd.concat(proteomics_dict.values())
phosphoproteomics = pd.concat(phosphoproteomics_dict.values())
clin = pd.concat(clin_dict.values())

#Remove -Inf
phosphoproteomics = phosphoproteomics.replace([np.inf, -np.inf], np.nan)


# Extract the levels of the MultiIndex columns
levels_df = pd.DataFrame(phosphoproteomics.columns.tolist(), columns=['gene', 'site', 'ENSG'])

# Create UID
levels_df["gene_site"] = levels_df["gene"] + "_" + levels_df["site"]
print(len(levels_df["gene_site"]))
print(len(levels_df["gene_site"].unique()))
print(len(levels_df["site"].unique()))

# Assign UID to data
phosphoproteomics.columns = levels_df["gene_site"]

# Transpose and export
phosphoproteomics.T.to_csv('../../output/data/CPTAC/phosphoproteomics_merged.csv')
proteomics.T.to_csv('../../output/data/CPTAC/proteomics_merged.csv')
transcriptomics.T.to_csv('../../output/data/CPTAC/transcriptomics_merged.csv')
levels_df.to_csv('../../output/data/CPTAC/phosphoproteomics_merged_identif.csv')
clin.to_csv('../../output/data/CPTAC/clin_merged.csv')