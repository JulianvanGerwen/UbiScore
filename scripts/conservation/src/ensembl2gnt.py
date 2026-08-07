#!/usr/bin/env python3

##### Background #####
# This file reads the big compara gene-tree file and breaks it up into smaller ones in subdirectories
# One file and subdirectory for each gene family
# This contains the species identity and other statistics for each gene in the gene family 


import os
import sys


# Define script to parse data
def parsed_tree(tree):
    """Simple tree parser - extracts tokens from tree string"""
    tab_tree = []
    word = ""
    tag = 0
    tag_bracket = 0
    
    for char in tree:
        if tag == 0:
            if char == "(":
                tag = 1
        if char == "[":
            if len(word) > 0:
                tab_tree.append(word)
            word = char
            tag_bracket = 1
        elif char == "]":
            word = word + char
            tab_tree.append(word)
            word = ""
            tag_bracket = 0
        elif tag_bracket == 1:
            word = word + char
        elif tag == 1:
            if char in ("(", ")", ","):
                tag = 1
                if word != "":
                    tab_tree.append(word)
                word = char
                tab_tree.append(word)
                word = ""
            elif char == ":":
                if word != "":
                    tab_tree.append(word)
                word = char
            elif char == " ":
                if word != "":
                    tab_tree.append(word)
                word = ""
            elif char.isalnum() or char in (".", "_") and char != ";":
                word = word + char
            elif char == ";":
                if word != "":
                    tab_tree.append(word)
                break
            else:
                word = word + char
    
    return tab_tree

# Define version
version = sys.argv[1].split(".")[-4]

# Initiliase variables
tree_name = ""
tree_number = 0
gene_specie_dict = {}
tag = 0

tree_directory = sys.argv[2]


# Define maximum name length
size_name = -1
# if len(sys.argv) > 2:
#     size_name = int(sys.argv[2])

# Load mega tree file Compara.XX.protein.nh.emf
tree_file = open(sys.argv[1], "r")

while 1:
    line = tree_file.readline()
    if line == "":
        break
    if line[0:2] == "//":
        tag = 0
    elif line[0:3] == "SEQ": # These lines denote the species and the corresponding gene
        # Load species definition
        line = line.rstrip()
        tab = line.split()
        specie = tab[1].split("_")[0]+"_"+tab[1].split("_")[1]
        gene_name = tab[2]
        gene_specie_dict[gene_name] = specie
    elif line[0:4] == "DATA": # These lines contain the actual homology data
        tree_number = tree_number+1
        tree_name = "ENSTREE_"+str(tree_number).zfill(5)
        supra_folder = "ENSTREE_"+str(int(str(tree_number).zfill(5)[0:2])+1).zfill(2)+"000"
        folder = "ENSTREE_"+str(tree_number).zfill(5)
        os.system("mkdir -p "+tree_directory+"/"+supra_folder+"/"+folder)

        tag = 1  # Set this so that we go to the next line
    elif tag == 1:
        line = line.rstrip()
        tab_tree = parsed_tree(line)
        # Write mapping
        new_tree = ""
        mapping_name = {}
        corresponding_file = open(tree_directory+"/"+supra_folder+"/"+folder+ \
                                  "/"+tree_name+"_treemapping.txt", "w")
        i = 0
        for item in tab_tree:
            #print(item, item.isalpha(), item.isdigit(), item.isalnum())
            if item[0].isalnum():   # 
                i = i+1
                gene_name = item
                gene_token = gene_name
                if size_name != -1:  # This just deals with if names are too big
                    gene_token = gene_name[0:size_name]+str(i).zfill(size_name)
                corresponding_file.write(gene_name+"\t"+gene_token+"\n")
                mapping_name[gene_name] = gene_token

                # gene_to_tree_file.write(gene_name+"\t"+tree_name+"\t"+ \
                #                         "ENSTREE_"+version+"/"+supra_folder+"/"+folder+"\n")
                sys.stdout.write(gene_name+"\t"+tree_name+"\t"+ \
                                 tree_directory+"/"+supra_folder+"/"+folder+"\n")
                new_tree = new_tree + gene_token  
            else:
                new_tree = new_tree + item
        corresponding_file.close()

        # Sort gene list
        gene_list_sorted = list(gene_specie_dict.keys())
        gene_list_sorted.sort()

        # Write .gnt file
        gnt_file = open(tree_directory+"/"+supra_folder+"/"+folder+"/"+tree_name+".gnt", "w")
        gnt_file.write("[\n")
        for gene in gene_list_sorted:
            specie_id = gene_specie_dict[gene]
            gnt_file.write(mapping_name[gene]+" "+specie_id.upper()+"\n")
        gnt_file.write("]\n")
        gnt_file.write(new_tree+";")   # The gene tree output
        gnt_file.close()


        # Write .tree file
        subtree_file = open(tree_directory+"/"+supra_folder+"/"+folder+"/"+tree_name+".tree", "w")
        subtree_file.write(new_tree+";")
        subtree_file.close()


        # Write genes / species file
        gnt_file = open(tree_directory+"/"+supra_folder+"/"+folder+"/"+tree_name+".list", "w")
        for gene in gene_list_sorted:
            specie_id = gene_specie_dict[gene]
            gnt_file.write(mapping_name[gene]+"\t"+specie_id+"\n")
        gnt_file.close()

        # Reset gene species dict (free memory)
        gene_specie_dict = {}
        
tree_file.close()
# gene_to_tree_file.close()
