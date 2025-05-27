import logging
import argparse

from pipemake_utils.logger import startLogger, logArgDict

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

def chuckArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-fasta", type = str, help = "Input FASTA file containing sequences to be chunked.", required = True)
    parser.add_argument("--chunk-size", type = int, help = "Size of each chunk in base pairs.", required = True)
    parser.add_argument("--output-fasta", type = str, help = "Output FASTA file to write the chunked sequences.", default = "chunked_sequences.fasta")
    return vars(parser.parse_args())

def chunkFasta(input_file, output_file, chunk_size):
    
    with open(output_file, 'w') as out_fasta:
        for record in SeqIO.parse(input_file, "fasta"):
            seq = record.seq
            for i in range(0, len(seq), chunk_size):
                chunk_seq = seq[i:i+chunk_size]
                chunk_id = f"{record.id}:{i + 1}-{i + chunk_size}"
                chunked_record = SeqRecord(chunk_seq, id=chunk_id, description="")
                SeqIO.write(chunked_record, out_fasta, "fasta")
                logging.info(f"Wrote chunk {chunk_id} from sequence {record.id}")
def main ():

    # Parse the arguments
    chunk_args = chuckArgs()

    # Start logger and log the arguments
    startLogger("chunk_fasta.log")
    logArgDict(chunk_args)

    chunkFasta(chunk_args['input_fasta'], chunk_args['output_fasta'], chunk_args['chunk_size'])

if __name__ == "__main__":
    main()