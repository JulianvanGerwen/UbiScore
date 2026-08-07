#!/usr/bin/python3

##### Background #####
# This file reads the big compara aa file and breaks it up into smaller ones in subdirectories
# One file and subdirectory for each gene family
# This contains the amino acid alignments for proteins in each gene family

import os
import sys

# Load arguments
version = sys.argv[1].split(".")[-4]
tree_directory = sys.argv[2]

# Load the big AA sequence alignemnt file Compara.XX.protein.aa.fasta
alignment_file = open(sys.argv[1], "r") 

# Start with first gene family (tree)
size_name = -1
tree_number = 1
tree_name = "ENSTREE_"+str(tree_number).zfill(5)
supra_folder = "ENSTREE_"+str(int(str(tree_number).zfill(5)[0:2])+1).zfill(2)+"000"
folder = "ENSTREE_"+str(tree_number).zfill(5)
os.system("mkdir -p "+tree_directory+"/"+supra_folder+"/"+folder)
aa_file = open(tree_directory+"/"+supra_folder+"/"+folder+"/"+tree_name+".aa.fasta", "w")

# Iterate over other gene families
i = 0
while 1:
    line = alignment_file.readline()
    i = i+1
    if line == "":
        break
    if line[0:2] == "//":
        aa_file.close()                                  # We close the previous file.
        tree_number = tree_number+1                      # We increase by +1.
        tree_name = "ENSTREE_"+str(tree_number).zfill(5) #
        supra_folder = "ENSTREE_"+str(int(str(tree_number).zfill(5)[0:2])+1).zfill(2)+"000"
        folder = "ENSTREE_"+str(tree_number).zfill(5)
        os.system("mkdir -p "+tree_directory+"/"+supra_folder+"/"+folder)
        aa_file = open(tree_directory + "/"+supra_folder+"/"+ \
                       folder+"/"+tree_name+".aa.fasta", "w") # We open a new file to write in.
    elif len(line) >= 1:
        aa_file.write(line)

aa_file.close()
alignment_file.close()
