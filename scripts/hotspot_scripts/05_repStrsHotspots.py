# Here I run hotspot analysis only on hotspot domains, along with their representative structures
# Must be run from my ptm_hotspots environment

# The positions that it gives for the representative structure is with respect to the domain sequence, not the PDB sequence that it is imbedded in

#Set up
import subprocess
import pandas as pd
import os
from Bio import SeqIO

#Functions
def sequenceIndexOverallMapper(sequence):
    # Dictionary to hold the mapping: key = residue index (ignoring gaps), value = overall index (including gaps)
    residue_mapping = {}

    # Counter for the residue index (only non-gap characters)
    res_index = 0

    # Iterate over the sequence with overall index starting at 1
    for overall_index, char in enumerate(sequence, start=1):
        if char != '-':
            res_index += 1
            residue_mapping[res_index] = overall_index
    
    #Invert
    residue_mapping_inv = {v: k for k, v in residue_mapping.items()}

    return(residue_mapping_inv)

def read_fasta_to_dict(fasta_file):
    fasta_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    return fasta_dict


#Set up p-value conversion factor
#This allows me to detect hotspots at the same threshold as without representative sequences
##Load all results
residueResAll = pd.read_csv("../../output/data/hotspots/domains/hotspotOutput/hotspotResidueOutput.csv")
##Get denominator for p-value conversion
pvalConversion = round(residueResAll["p_adjust"].div(residueResAll["pvals"]).max())
with open("../../output/data/hotspots/domains/hotspotOutput/pvalConversion.txt", "w") as f:
    f.write(str(pvalConversion))


#Run hotspot detection
#Run on sequence alignments run with representative sequences, but with represnetative sequences dropped. This ensures statistics are the same as without representative sequences
commandHotspots = [
    "python", "../ptm_hotspots/ptm_hotspots.py",
    "-o", "../../output/data/hotspots/domains/representativeStrs/hotspotOutput.csv",
    "--dir", "../../output/data/hotspots/domains/representativeStrs/domainsPostAlignUnwrap_noRep",
    "--ptmfile", "../../output/data/hotspots/ubnSites_forHotspots",
    "--aa_letters", "K",
    "--useRawPvals",
    "--threshold", str(0.01/pvalConversion)
]
commandResidues = [
    "python", "../ptm_hotspots/ptm_hotspots.py",
    "-o", "../../output/data/hotspots/domains/representativeStrs/hotspotResidueOutput.csv",
    "--dir", "../../output/data/hotspots/domains/representativeStrs/domainsPostAlignUnwrap_noRep",
    "--ptmfile", "../../output/data/hotspots/ubnSites_forHotspots",
    "--aa_letters", "K",
    "--printSites",
    "--useRawPvals",
    "--threshold", str(0.01/pvalConversion)
]
subprocess.run(commandHotspots, check=True)
subprocess.run(commandResidues, check=True)


#Map back in representative structure sequences
#Read hotspot results and subset for distinct
residueRes = pd.read_csv("../../output/data/hotspots/domains/representativeStrs/hotspotResidueOutput.csv")
cols = ["domain", "position_aln", "foreground", "pvals", "bg_means", "bg_stdev", "p_adjust", "hotspot"]
residueRes = residueRes[cols].drop_duplicates()

#Loop over domains and add representative structure residues
##Read MAFFT results including representative structure
for domainID in residueRes['domain'].unique():
    repSeqs = read_fasta_to_dict("../../output/data/hotspots/domains/representativeStrs/domainsPostAlignUnwrap/" + domainID + ".fasta")
    repSeqID = list(repSeqs.keys())[-1]
    repSeq = repSeqs[repSeqID]
    repSeqMap = sequenceIndexOverallMapper(repSeq)
    ##Subset residue data and map in representative structure index
    residueDat = residueRes[residueRes["domain"] == domainID]
    residueDat["position_representative"] = residueDat["position_aln"].map(repSeqMap)
    ##Re-adjust p-values using the proper converison
    residueDat["p_adjust_og"] = residueDat["pvals"]*pvalConversion
    residueDat.loc[residueDat["p_adjust_og"] > 1, "p_adjust_og"] = 1
    residueDat.to_csv("../../output/data/hotspots/domains/representativeStrs/repStrOutput/" + domainID + ".csv")
##Combine into one file
input_dir = "../../output/data/hotspots/domains/representativeStrs/repStrOutput"
output_file = "../../output/data/hotspots/domains/representativeStrs/hotspotResidueOutput_strRep.csv"
csv_files = [os.path.join(input_dir, f) for f in os.listdir(input_dir) if f.endswith('.csv')]
# Read each CSV file into a DataFrame and store them in a list
dataframes = [pd.read_csv(csv_file).iloc[:, 1:] for csv_file in csv_files]
# Combine all DataFrames into one
combined_df = pd.concat(dataframes, ignore_index=True)
# Export the combined DataFrame to a CSV file
combined_df.to_csv(output_file, index=False)