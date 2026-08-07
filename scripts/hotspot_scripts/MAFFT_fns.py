#This script stores functions for MAFFT
#It runs MAFFT to align sequences, and then unwraps the output so that one sequence is one line
#It should be used with my mafft conda anvironment

import os
import glob
import subprocess
import argparse
from Bio import SeqIO
from Bio.SeqIO.FastaIO import FastaWriter

#Function for running MAFFT
def runMAFFT(input_dir, output_dir):
    
    # Ensure the output directory exists
    os.makedirs(output_dir, exist_ok=True)
    
    # Find all .fasta files in the input directory
    fasta_files = glob.glob(os.path.join(input_dir, "*.fasta"))
    
    # Loop over each fasta file
    for file in fasta_files:
        base = os.path.basename(file)
        out_file = os.path.join(output_dir, base)
        print(f"Processing {file} -> {out_file}")
        
        # Run 'ginsi' on the file and write the output to the out_file
        with open(out_file, "w") as out_handle:
            subprocess.run(["ginsi", file], stdout=out_handle, check=True)


#Function to process output
def process_fasta_files(input_dir, output_dir):
    # Create output directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Process each FASTA file in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            in_path = os.path.join(input_dir, filename)
            out_path = os.path.join(output_dir, filename)
            # Parse the FASTA file into records
            records = list(SeqIO.parse(in_path, "fasta"))
            with open(out_path, "w") as out_handle:
                # Use FastaWriter with wrap=0 to write sequences on a single line
                writer = FastaWriter(out_handle, wrap=0)
                writer.write_file(records)
            print(f"Processed: {filename}")

#Function to process output and drop the last sequence. Useful for leaving representative sequence out of hotspot analysis
def dropLast_fasta_files(input_dir, output_dir):
    # Create output directory if it does not exist
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    # Process each FASTA file in the input directory
    for filename in os.listdir(input_dir):
        if filename.endswith(".fasta") or filename.endswith(".fa"):
            in_path = os.path.join(input_dir, filename)
            out_path = os.path.join(output_dir, filename)
            # Parse the FASTA file into records
            records = list(SeqIO.parse(in_path, "fasta"))
            # Drop the last sequence, if present
            if records:
                records = records[:-1]
            with open(out_path, "w") as out_handle:
                # Use FastaWriter with wrap=0 to write sequences on a single line
                writer = FastaWriter(out_handle, wrap=0)
                writer.write_file(records)
            print(f"Processed: {filename}")


#Function to run MAFFT and then unravel output
def runMAFFT_and_unravel(input_dir, output_dirMAFFT, output_dirUnrav):
    runMAFFT(input_dir = input_dir, output_dir = output_dirMAFFT)
    process_fasta_files(input_dir = output_dirMAFFT, output_dir = output_dirUnrav)