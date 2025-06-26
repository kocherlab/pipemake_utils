import argparse

from Bio import SeqIO
from collections import defaultdict

def filterMinLength(seq_lengths, min_length):
    filtered_lengths = {k: v for k, v in seq_lengths.items() if v >= min_length}
    return filtered_lengths

def processAssembly(input_filename, output_filename, min_length = None, chrom_count = None, mt_filename = '', no_mt = None, input_format = "fasta", output_format = "fasta", species = ''):

    # Create a dictionary to store the sequence lengths
    seq_lengths = defaultdict(int)

    # Read in the sequences from the input file
    with open(input_filename, 'r') as input_file:
        for record in SeqIO.parse(input_file, input_format):
            seq_lengths[record.id] = len(record.seq)

    # Filter by minimum length if specified
    if min_length is not None:
        seq_lengths = filterMinLength(seq_lengths, min_length)

    # Sort the sequences by length, only keeping the id
    seq_lengths = sorted(seq_lengths, key=seq_lengths.get, reverse=True)

    # Filter by chromosome count if specified
    if chrom_count is not None:
        seq_lengths = list(seq_lengths)[:chrom_count]

    # Create a dictionary to store the ids with chromosome names
    chrom_ids = {}

    # Assign the chromosome names
    for chrom_int, seq_id in enumerate(seq_lengths, 1):
        chrom_id = f"chr_{chrom_int}"
        if species:
            chrom_id = f"{species}_" + chrom_id
        chrom_ids[seq_id] = chrom_id

    # Index the input file
    input_file = SeqIO.index(input_filename, input_format)

    # Loop through the sequence lengths and write to the output file
    with open(output_filename, 'w') as output_file:
        for seq_id in seq_lengths:
            record = input_file[seq_id]
            record.id = chrom_ids[seq_id]
            record.description = ''
            record.name = ''
            SeqIO.write(record, output_file, output_format)

        # Check if there is a mitochondrial sequence to add
        if not no_mt:

            # Open the mitochondrial input file using seqIO
            mt_chr = list(SeqIO.parse(mt_filename, input_format))

            # Confirm a single mitochondrial sequence
            if len(mt_chr) != 1:
                raise ValueError("Mitochondrial input file must contain a single sequence")
            
            # Write the mitochondrial sequence
            SeqIO.write(mt_chr, output_file, output_format)
            
def processAssemblyParser():
    parser = argparse.ArgumentParser(description="Process an assembly file.")
    parser.add_argument("--input", dest="input_filename", type=str, help="Input assembly file", required=True)
    parser.add_argument("--output", dest="output_filename", type=str, help="Output assembly file", default="output.fasta")
    parser.add_argument("--species", type=str, default='', help="Species id")
    parser.add_argument("--min-length", type=int, default=None, help="Minimum length of sequences to keep")
    parser.add_argument("--chrom-count", type=int, default=None, help="Number of longest sequences to keep")
    parser.add_argument("--input-format", type=str, default="fasta", help="Input file format")
    parser.add_argument("--output-format", type=str, default="fasta", help="Output file format")
    mt_parser = parser.add_mutually_exclusive_group(required = True)
    mt_parser.add_argument("--mt-input", dest="mt_filename", type=str, help = "Mitochondrial sequence input file")
    mt_parser.add_argument("--no-mt", help = "Indicates there is no mitochondrial sequence", action = "store_true")
    return vars(parser.parse_args())

if __name__ == "__main__":
    args = processAssemblyParser()
    processAssembly(**args)