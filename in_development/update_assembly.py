import argparse

from Bio import SeqIO
from collections import defaultdict

from pipemake_utils.misc import confirmFile
from pipemake_utils.logger import startLogger, logArgDict


def updateParser():
    parser = argparse.ArgumentParser(description="Update an assembly FASTA file")
    parser.add_argument("--assembly", dest = 'input_filename', help = "Assembly input in FASTA format", type = str, required = True, action = confirmFile())
    parser.add_argument("--output", dest = 'output_filename', help = "Assembly output filename", type = str, required = True)
    parser.add_argument("--species", type=str, default='', help="Species id")
    remove_headers = parser.add_mutually_exclusive_group(required=True)
    remove_headers.add_argument("--remove-headers", help = "Headers to remove from the assembly", type = str, nargs='*', default=[])
    remove_headers.add_argument("--remove-headers-file", help = "File containing headers to remove", type = str, action = confirmFile())

    return vars(parser.parse_args())

def assignHeadersToRemove(remove_headers = [], remove_headers_file = None):
    if remove_headers_file:
        with open(remove_headers_file, 'r') as f:
            remove_headers.extend([line.strip() for line in f if line.strip()])
    elif remove_headers:
        remove_headers = [header.strip().replace(',', '') for header in remove_headers if header.strip()]
    return remove_headers

def updateAssembly(input_filename, output_filename, remove_headers=[], species = '', file_format="fasta", **kwargs):

    # Create a dictionary to store the sequence lengths
    seq_lengths = defaultdict(int)

    # Read in the sequences from the input file
    with open(input_filename, 'r') as input_file:
        for record in SeqIO.parse(input_file, file_format):
            seq_lengths[record.id] = len(record.seq)

    # Remove specified headers
    if remove_headers:
       for header in remove_headers:
           if header not in seq_lengths:
               raise ValueError(f"Header '{header}' not found in the input file.")
           del seq_lengths[header]

    # Sort the sequences by length, only keeping the id
    seq_lengths = sorted(seq_lengths, key=seq_lengths.get, reverse=True)

    # Create a dictionary to store the ids with chromosome names
    chrom_ids = {}

    # Assign the chromosome names
    for chrom_int, seq_id in enumerate(seq_lengths, 1):
        chrom_id = f"chr_{chrom_int}"
        if species:
            chrom_id = f"{species}_" + chrom_id
        chrom_ids[seq_id] = chrom_id

    # Index the input file
    input_file = SeqIO.index(input_filename, file_format)

    # Loop through the sequence lengths and write to the output file
    with open(output_filename, 'w') as output_file:
        for seq_id in seq_lengths:
            print(f"Processing sequence: {seq_id} with new id: {chrom_ids[seq_id]}")
            record = input_file[seq_id]
            record.id = chrom_ids[seq_id]
            record.description = ''
            record.name = ''
            SeqIO.write(record, output_file, file_format)

def main():
    # Parse the arguments
    update_args = updateParser()

    # Start the logger
    startLogger()

    # Log the arguments
    logArgDict(update_args)

    # Assign headers to remove
    update_args["remove_headers"] = assignHeadersToRemove(update_args["remove_headers"], update_args["remove_headers_file"])

    # Update the assembly
    updateAssembly(**update_args)

if __name__ == "__main__":
    main()
