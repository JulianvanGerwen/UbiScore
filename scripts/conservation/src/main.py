#!/usr/bin/env python3
"""
Contains main scripts for processing alignments
"""

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from pathlib import Path

def create_position_mapping(aligned_sequence):
    """
    Create mapping between alignment positions and sequence positions.
    
    This handles gaps in the alignment - maps sequence positions (without gaps)
    to alignment column positions (with gaps).
    
    Returns:
        dict: {alignment_column: sequence_position} - Maps alignment positions to sequence positions
    """
    alignment_to_sequence = {}
    sequence_pos = 0
    
    for alignment_pos, amino_acid in enumerate(aligned_sequence):
        if amino_acid != "-":  # Skip gap characters
            alignment_to_sequence[alignment_pos] = sequence_pos
            sequence_pos += 1
    
    return alignment_to_sequence

def findPtmColumns(sequences, position_mappings, ptm_sites, ptm_acceptors):
    """
    Find alignment columns that contain PTM sites.

    sequences: Dictionary of aligned sequences
    position_mappings: Dictionary of position mappings for these sequences
    ptm_sites: Dictionary of all ptm sites per gene
    ptm_acceptors: List of amino acids considered PTM acceptors
    
    Returns:
        list: Alignment column numbers that contain verified PTM sites
    """
    ptm_columns = []
    
    # Iterate over sequences
    for gene_id, sequence in sequences.items():
        if gene_id not in ptm_sites:
            continue
        
        # Get position mapping for this sequence
        pos_mapping = position_mappings[gene_id]
        
        # Check each PTM site for this gene
        for site_position in ptm_sites[gene_id]:
            try:
                # Convert 1-based site position to 0-based sequence position
                seq_pos = int(site_position) - 1
                
                # Find corresponding alignment column
                alignment_col = None
                for col, seq_pos_mapped in pos_mapping.items():
                    if seq_pos_mapped == seq_pos:
                        alignment_col = col
                        break
                
                if alignment_col is not None:
                    # Verify the amino acid is a valid acceptor
                    amino_acid = sequence[alignment_col]
                    if amino_acid in ptm_acceptors:
                        ptm_columns.append(alignment_col)
            
            except (ValueError, IndexError):
                continue
    
    # Remove duplicates and sort
    return sorted(list(set(ptm_columns)))

def calculatePtmScores(sequences, position_mappings, ptm_sites, ptm_acceptors, gene_id, region_window, central_column):
    """
    Calculate PTM scores for a gene across a region around a central PTM site.
    
    Scoring rules:
    1. Known PTM site (experimental): score = 1.0
    3. Acceptor amino acid (K for ubiquitin, S/T/Y for phospho): score = 0.5
    4. Non-acceptor amino acid: score = 0.0
    
    Returns:
        float: Maximum score in the region around the central column
    """
    scores = []
    
    # Scan region around central column
    for column in range(central_column - region_window, central_column + region_window + 1):
        if column < 0 or column >= len(sequences[gene_id]):
            continue
        
        if column not in position_mappings[gene_id]:
            continue
        
        amino_acid = sequences[gene_id][column]
        sequence_position = str(position_mappings[gene_id][column] + 1)  # Convert to 1-based

        # Calculate score for this position
        score = 0.0
        
        # Check if it's a known PTM site (highest priority)
        if gene_id in ptm_sites and sequence_position in ptm_sites[gene_id]:  # Need to make sure ptm_sites entries are strings
            score = 1.0
        # Check if it's an acceptor amino acid
        elif amino_acid in ptm_acceptors:
            score = 0.5
        # Otherwise, score is 0.0
        
        scores.append(score)
    
    # Return maximum score in the region
    return max(scores) if scores else 0.0


def writePtmContinueFiles(family_id, tree_directory, species_tag, region_window, ptm_columns, scoring_dict, sequences):
    """
    Write input files for ancestral state reconstruction.
    
    This creates:
    1. A tab-separated file with PTM scores for each gene at each PTM column
    
    Args:
        family_id: Gene family identifier (e.g., "ENSTREE_00001")
        tree_directory: Base directory containing gene family data
        species_tag: Species dataset tag (e.g., "species_n10")
        region_window: Window size around PTM site
        ptm_columns: List of alignment columns containing PTM sites
        scoring_dict: Dictionary mapping {gene_id: {column: score}}
        sequences: Dictionary of aligned sequences
    """
    # Calculate directory structure
    tree_number = int(family_id.split("_")[1])
    supra_folder = f"ENSTREE_{str(int(str(tree_number).zfill(5)[0:2]) + 1).zfill(2)}000"
    
    # Create output directory
    output_dir = Path(tree_directory) / supra_folder / family_id / species_tag / f"region_w{region_window}"
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Base filename for output files
    basename = output_dir / family_id
    
    # Write PTM score file for ancestral reconstruction
    score_file = f"{basename}_ptm_continue_region_w{region_window}.txt"
    
    # Convert 0-based column indices to 1-based for output
    real_column_with_ptm = [str(int(column) + 1) for column in ptm_columns]
    
    with open(score_file, 'w') as score_out:
        # Write header
        score_out.write("Gene\t" + "\t".join(real_column_with_ptm) + "\n")
        
        # Write scores for each gene
        for gene_id in sequences:
            line_parts = [gene_id]
            for column in ptm_columns:
                score = "0.0"
                if gene_id in scoring_dict and column in scoring_dict[gene_id]:
                    score = str(scoring_dict[gene_id][column])
                line_parts.append(score)
            score_out.write("\t".join(line_parts) + "\n")


def get_species_name_fallback(gene_id):
    """Extract species name from Ensembl gene ID"""
    if gene_id.startswith("ENSP0"):
        return "Homo sapiens"
    elif gene_id.startswith("ENSMUSP"):
        return "Mus musculus"
    elif gene_id.startswith("ENSRNOP"):
        return "Rattus norvegicus"
    elif gene_id.startswith("ENSGALP"):
        return "Gallus gallus"
    elif gene_id.startswith("ENSXETP"):
        return "Xenopus tropicalis"
    elif gene_id.startswith("ENSDARP"):
        return "Danio rerio"
    elif gene_id.startswith("ENSDMEP"):
        return "Drosophila melanogaster"
    elif gene_id.startswith("ENSCESP"):
        return "Saccharomyces cerevisiae"
    else:
        return "Unknown species"



def remove_all_gap_columns(input_records, output_file):
    """
    Remove alignment columns that consist entirely of '-' from a FASTA file.

    Args:
        input_records: Fasta file as a list, read in using SeqIO.parse
        output_file: File to save the final fasta in
        returns - output_records: Fasta file as a list, read in using SeqIO.parse
    """
    # Prep sequences
    #records = list(SeqIO.parse(input_fasta, "fasta"))
    ids = [r.id for r in input_records]
    seqs = [str(r.seq) for r in input_records]

    # Check all sequences have same length
    lengths = {len(s) for s in seqs}
    if len(lengths) != 1:
        raise ValueError(f"Sequences have inconsistent lengths: {lengths}")

   # Identify which columns to keep
    keep = [any(c != '-' for c in col) for col in zip(*seqs)]

    # Filter sequences
    filtered_seqs = [''.join(ch for ch, k in zip(seq, keep) if k) for seq in seqs]

    # Build new SeqRecord objects
    output_records = [SeqRecord(Seq(s), id=i, description="") for i, s in zip(ids, filtered_seqs)]

    # Write FASTA output
    with open(output_file, "w") as fh:
        for rec in output_records:
            fh.write(f">{rec.id}\n{str(rec.seq)}\n")

    return output_records