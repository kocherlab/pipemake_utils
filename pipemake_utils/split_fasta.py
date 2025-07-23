import os
import logging
import argparse

from pipemake_utils.logger import startLogger, logArgDict

from Bio import SeqIO

def splitArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fasta", type = str, help = "Input FASTA file containing sequences to be chunked.", required = True)
    parser.add_argument("--output-dir", type = str, help = "Output directory to write the chunked sequences.", default = "split_sequences")
    return vars(parser.parse_args())

def splitFasta(input_file, output_dir):
    for record in SeqIO.parse(input_file, "fasta"):
        os.makedirs(output_dir, exist_ok=True)
        with open(os.path.join(output_dir, f"{record.id}.fasta"), 'w') as out_fasta:
            SeqIO.write(record, out_fasta, "fasta")
            logging.info(f"Wrote sequence {record.id} to {os.path.join(output_dir, f'{record.id}.fasta')}")

def main ():

    # Parse the arguments
    split_args = splitArgs()

    # Start logger and log the arguments
    startLogger(f"{split_args['output_dir']}.log")
    logArgDict(split_args)

    splitFasta(split_args['input_fasta'], split_args['output_dir'])

if __name__ == "__main__":
    main()