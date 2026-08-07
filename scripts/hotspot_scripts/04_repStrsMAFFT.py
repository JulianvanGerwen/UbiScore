#This script adds domain sequences for represntative strucutres to the corresponding hotspot domain alignments and re-runs MAFFT on it
# This must be run from my mafft environment

import os
import glob
import subprocess
import argparse
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import MAFFT_fns

#Functions
def read_fasta_to_dict(fasta_file):
    fasta_dict = {record.id: str(record.seq) for record in SeqIO.parse(fasta_file, "fasta")}
    return fasta_dict


#Load and process representative sequences
repSeqs = read_fasta_to_dict("../../output/data/hotspots/domains/representativeStrs/strRep_seqDomain.fasta")

#Prepare output directory
outputDir = "../../output/data/hotspots/domains/representativeStrs/domainsPreAlign"
if not os.path.exists(outputDir):
    os.makedirs(outputDir)

#Loop over representative sequences and prepare for MAFFT
for id in list(repSeqs.keys()):
    domain, seqId = id.split("__")
    ##Read in corresponding sequences for the domain and add this one
    domainSeqs = read_fasta_to_dict("../../output/data/hotspots/domains/domainsPreAlign/" + domain + ".fasta")
    domainSeqs[seqId] = repSeqs[id]
    ##Write to a new file
    fasta_file = outputDir + "/" + domain + ".fasta"
    with open(fasta_file, "w") as f:
        for seq_id, sequence in domainSeqs.items():
            seq_record = SeqRecord(Seq(sequence), id=seq_id, description="")
            SeqIO.write(seq_record, f, "fasta")

#Run MAFFT
MAFFT_fns.runMAFFT_and_unravel("../../output/data/hotspots/domains/representativeStrs/domainsPreAlign", "../../output/data/hotspots/domains/representativeStrs/domainsPostAlign", "../../output/data/hotspots/domains/representativeStrs/domainsPostAlignUnwrap")

#Also drop the last item, which is the representative sequence
MAFFT_fns.dropLast_fasta_files("../../output/data/hotspots/domains/representativeStrs/domainsPostAlign", "../../output/data/hotspots/domains/representativeStrs/domainsPostAlignUnwrap_noRep")