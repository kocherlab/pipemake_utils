import os
import argparse

from collections import defaultdict

from pipemake_utils.misc import confirmFile
from pipemake_utils.logger import startLogger, logArgDict
from pipemake_utils.subset_fasta import subsetFasta

def subsetParser ():

    # Create an argument parser for the subset FASTA script
    subset_parser = argparse.ArgumentParser()
    subset_parser.add_argument('--fasta-list', help = 'File containing the paths to the sample FASTA files', type = str, required = True, action = confirmFile())
    subset_parser.add_argument('--loci-file', help = 'File containing the loci for each sample. Formatted as "sample: locus1, locus2, locus3:start-end". Coordinates are 1-based, half-open', type = str, required = True, action = confirmFile())
    subset_parser.add_argument('--output-dir', help = 'Directory to write the subset FASTA files to', type = str, default = 'Subset_FASTAs')

    return vars(subset_parser.parse_args())

def readLociFile(loci_filename):

    # Create a dictionary to hold the loci for each sample
    loci_dict = defaultdict(set)

    # Read the loci file and populate the dictionary with sets of loci for each sample
    with open(loci_filename, 'r') as loci_file:
        for loci_line in loci_file:
            sample_name, loci = loci_line.strip().split(': ')
            loci_dict[sample_name].update(loci.split(', '))

    return loci_dict

def readFastaList(fasta_list_filename):

    # Create a list to hold the paths to the sample FASTA files
    fasta_list = []

    # Read the FASTA list file and populate the list with the paths to the sample FASTA files
    with open(fasta_list_filename, 'r') as fasta_list_file:
        for fasta_list_line in fasta_list_file:

            # Assign the path to the sample FASTA file
            fasta_filename = fasta_list_line.strip()

            # Confirm that the input exists
            if not os.path.exists(fasta_filename):
                raise ValueError(f"{fasta_filename} does not exist")
            
            # Confirm that the input is a file
            if not os.path.isfile(fasta_filename):
                raise ValueError(f"{fasta_filename} is not a file")
            
            # Append the path to the sample FASTA file to the list
            fasta_list.append(fasta_filename)

    return fasta_list

def main():

    # Assign the arguments from the CLI
    subset_args = subsetParser()

    # Start the logger
    startLogger()
    logArgDict(subset_args)

    # Create a dictionary mapping sample names to sets of loci
    loci_dict = readLociFile(subset_args['loci_file'])
  
    # Subset the FASTA files from the input list
    for sample_filename in readFastaList(subset_args['fasta_list']):

        # Assign the sample name by splitting the filename at '_genome'
        sample_name = ''

        # Assign the sample name by finding the first locus in the loci dictionary that is found in the filename
        for locus in loci_dict:
            if locus not in sample_filename:
                continue
            if not sample_name:
                sample_name = locus
            else:
                raise ValueError(f"Multiple samples associated with {sample_filename} in the loci dictionary")

        # Report an error if the sample name is not found in the loci dictionary
        if sample_name not in loci_dict:
            raise ValueError(f"Sample name {sample_name} not found in loci dictionary")
        
        # Create the output directory if it doesn't exist
        os.makedirs(subset_args['output_dir'], exist_ok=True)

        # Assign the output filename for the subset FASTA file
        subset_filename = os.path.join(subset_args['output_dir'], f"{sample_name}.fa")

        # Subset the FASTA file for the sample
        subsetFasta(sample_filename, loci_dict[sample_name], subset_filename)


if __name__ == "__main__":
    main()    