import argparse

import pandas as pd

from Bio import SeqIO
from collections import defaultdict

def updateParser():
    parser = argparse.ArgumentParser(description="Update an assembly FASTA file")
    parser.add_argument("--assembly", dest = 'input_filename', help = "Assembly input in FASTA format", type = str, required = True)
    parser.add_argument("--output", dest = 'output_filename', help = "Assembly output filename", type = str, required = True)
    parser.add_argument("--blob", dest = 'blob_filename', help = "BlobTools blob file in TSV format", type = str, required = True)
    parser.add_argument("--phylum", dest = 'phylum_to_keep', type=str, default='', help="Phylum to keep in the assembly")

    return vars(parser.parse_args())

def updateAssembly(input_filename, output_filename, blob_filename, phylum_to_keep, file_format="fasta", **kwargs):

    # Open the blob tsv uing the pandas library
    blob_df = pd.read_csv(blob_filename, sep='\t', header=0)

    # Remove any sequences that are not in the specified phylum
    blob_df = blob_df[blob_df['phylum'] == phylum_to_keep]

    # Loop through the sequence lengths and write to the output file
    with open(output_filename, 'w') as output_file:
        with open(input_filename, 'r') as input_file:
            for record in SeqIO.parse(input_file, file_format):
                if record.id in blob_df['# record'].values:
                    SeqIO.write(record, output_file, file_format)

def main():
    # Parse the arguments
    update_args = updateParser()

    # Update the assembly
    updateAssembly(**update_args)

if __name__ == "__main__":
    main()