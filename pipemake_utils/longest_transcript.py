import re
import argparse

from Bio import SeqIO

from pipemake_utils.logger import startLogger, logArgDict
from pipemake_utils.seqDB import sniffSeqFile, DBFileReader


def seqParser():
    # Create the argument parser
    seq_parser = argparse.ArgumentParser(description="Process NCBI fasta files")

    # Input and output arguments
    db_choices = ["ncbi", "flybase", "pipemake"]
    seq_parser.add_argument(
        "--database",
        help="Source database for the input file",
        required=True,
        choices=db_choices,
    )
    seq_parser.add_argument("--input-filename", help="Input filename", required=True)
    seq_parser.add_argument(
        "--output-prefix", help="Output filename prefix", required=True
    )
    seq_parser.add_argument("--input-format", help="Input format", default="fasta")
    type_choices = ["CDS", "AA"]
    seq_parser.add_argument("--input-type", help="Input type", choices=type_choices)
    seq_parser.add_argument(
        "--output-type", help="Output type", choices=type_choices, default="AA"
    )

    # Optional arguments
    seq_parser.add_argument(
        "--output-primary-id", help="Output primary id", default="protein_id"
    )

    return vars(seq_parser.parse_args())


def returnGeneID(record, primary_id):
    # Check if the gene is the primary id
    if primary_id == "gene":
        return record.id

    # Split by [ and ] to get the attributes
    for _s in re.split(r"\[|\]", record.description):
        # Skip if the string is not an attribute
        if "=" not in _s.strip():
            continue
        if "gene=" in _s:
            return _s.split("=")[1]

    raise ValueError(f"Gene ID not found in record: {record}")


def main():
    # Assign the arguments from the command-line
    seq_args = seqParser()

    # Start the logger
    startLogger(f"{seq_args['output_prefix']}.log")

    # Log the arguments
    logArgDict(seq_args)

    # Check if the input type is provided
    if not seq_args["input_type"]:
        # Create a dictionary to convert the input type
        input_type_cvt = {"dna": "CDS", "protein": "AA"}

        # Sniff the sequence file, if not provided
        seq_args["input_type"] = input_type_cvt[
            sniffSeqFile(seq_args["input_filename"])
        ]

    # Warn if the input cannot create the CDS
    if seq_args["input_type"] == "AA" and seq_args["output_type"] == "CDS":
        raise ValueError("Cannot create DNA from protein sequence")

    # Create a dict to store the longest transcript
    longest_transcript = {}

    # Open the output file
    with open(f"{seq_args['output_prefix']}.fa", "w") as output_file:
        for record in DBFileReader.read(
            seq_args["database"],
            seq_args["input_filename"],
            file_type=seq_args["input_type"],
            file_format=seq_args["input_format"],
            primary_id=seq_args["output_primary_id"],
        ):
            # Assign the gene id
            gene_id = returnGeneID(record, seq_args["output_primary_id"])

            # Check if the gene id is in the longest transcript
            if gene_id not in longest_transcript:
                longest_transcript[gene_id] = record
            elif len(record.seq) > len(longest_transcript[gene_id].seq):
                longest_transcript[gene_id] = record

        # Write the longest transcript to the output file
        for record in longest_transcript.values():
            # Translate the sequence if needed
            if seq_args["input_type"] == "CDS" and seq_args["output_type"] == "AA":
                record.seq = record.seq.translate()

            # Write the record to the output file
            SeqIO.write(record, output_file, "fasta")


if __name__ == "__main__":
    main()
