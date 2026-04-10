import os
import argparse

from Bio import SeqIO

from collections import defaultdict

from pipemake_utils.misc import confirmFile
from pipemake_utils.logger import startLogger, logArgDict

def subsetParser ():

    # Create an argument parser for the subset FASTA script
    subset_parser = argparse.ArgumentParser()
    subset_parser.add_argument('--fasta-file', help = 'File containing the paths to the sample FASTA files', type = str)#, required = True, action = confirmFile())
    subset_parser.add_argument('--loci-str', help = 'Formatted as "locus1,locus2,locus3:start-end". Coordinates are 1-based, half-open', type = str, required = True)
    subset_parser.add_argument('--output-filename', help = 'File to write the subset FASTA file to', type = str, default = 'Subset_FASTA.fa')

    return vars(subset_parser.parse_args())

def subsetFasta(input_filename, subset_loci, subset_filename):

    # Create a dictionary to hold the loci to subet
    subset_loci_dict = defaultdict(set)

    # Loop the loci and assign the locus name and coordinates as (-1, -1) if no coordinates are provided, or as (start, end) if coordinates are provided
    for sample_locus in subset_loci:
        if ':' not in sample_locus:
            subset_loci_dict[sample_locus] = (-1, -1)
        else:
            locus_name, locus_coords = sample_locus.split(':')
            start, end = map(int, locus_coords.split('-'))
            subset_loci_dict[locus_name] = (start, end)

    # Write the subset FASTA file for the sample, including only records with IDs in the set of loci
    with open(subset_filename, 'w') as subset_fasta:
        for record in SeqIO.parse(input_filename, 'fasta'):
            if record.id in subset_loci_dict:
                if subset_loci_dict[record.id] != (-1, -1):
                    start, end = subset_loci_dict[record.id]
                    record.seq = record.seq[start - 1:end]
                SeqIO.write(record, subset_fasta, 'fasta')

def main():

    # Assign the arguments from the CLI
    subset_args = subsetParser()

    # Start the logger
    startLogger()
    logArgDict(subset_args)

    # Parse the loci string into a tuple of loci
    loci_to_subset = [_locus.strip() for _locus in subset_args['loci_str'].split(',')]

    # Subset the FASTA file
    subsetFasta(subset_args['fasta_file'], loci_to_subset, subset_args['output_filename'])
if __name__ == "__main__":
    main()    