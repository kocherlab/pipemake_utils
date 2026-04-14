import sys

import pandas as pd

# Create a list to store the longest protein ids
longest_protein_ids = []

# Get the gene list filename from the command line arguments
gene_list_filename = sys.argv[1]

# Read in the gene list file
with open(gene_list_filename, 'r') as gene_list_file:
    gene_list = [line.strip() for line in gene_list_file.readlines()]

# Get the sequence filename from the command line arguments
sequence_filename = sys.argv[2]

# Read in the sequence file and check if the id is in the gene list
with open(sequence_filename, 'r') as sequence_file:
    for sequence_line in sequence_file:
        if not sequence_line.startswith('>'):
            continue
        sequence_list = sequence_line[1:].strip().split()

        # Check if the id is in the gene list
        if sequence_list[0] not in gene_list:
            continue

        # using this example - MGEN-H1_v3.3_006397 [protein_id=MGEN-H1_v3.3_006397-R1] [product=kinesin-like protein KIF20A] get the protein_id
        longest_protein_id = ''
        for sequence_attrribute in sequence_list[1:]:
            sequence_attrribute = sequence_attrribute.strip().replace('[', '').replace(']', '')
            if not sequence_attrribute.startswith('protein_id='):
                continue
            longest_protein_id = sequence_attrribute.split('=')[1]

        if not longest_protein_id:
            raise ValueError(f'Could not find protein_id for {sequence_list[0]}')

        longest_protein_ids.append(longest_protein_id)

# Assing the GFF filename from the command line arguments
gff_filename = sys.argv[3]

# Read in the GFF file and filter for the longest protein ids
gff_df = pd.read_csv(gff_filename, sep='\t', header=None, comment='#')

# Filter the GFF dataframe to only include mRNA or transcript features
gff_df = gff_df[gff_df[2].isin(['mRNA', 'transcript'])]

# Create a new column for the feature ID
gff_df['feature_id'] = gff_df[8].str.extract(r'ID=([^;]+)')

# Check if the feature ID is in the longest protein ids list
gff_df = gff_df[gff_df['feature_id'].isin(longest_protein_ids)]

# Adjust the start and end coordinates to be 0-based and half-open
gff_df[3] = gff_df[3] - 1

# Convert the filtered GFF to bed format and write to a new file
bed_df = gff_df[[0, 3, 4, 'feature_id', 6]]
bed_df.to_csv('longest_proteins.bed', sep='\t', header=False, index=False)git 