import argparse

from pipemake_utils.misc import confirmFile
from pipemake_utils.align import mapGapsFromIndex
from pipemake_utils.logger import startLogger, logArgDict


def codonAlignmentParser():
    # Create the parser
    codon_alignment_parser = argparse.ArgumentParser(
        description="Translate a sequence file"
    )

    # Add the arguments
    codon_alignment_parser.add_argument(
        "--aligned-aa",
        help="Input aligned amino acid sequence file",
        required=True,
        action=confirmFile(),
    )
    codon_alignment_parser.add_argument(
        "--unaligned-cds",
        help="Input unaligned CDS sequence file",
        required=True,
        action=confirmFile(),
    )
    codon_alignment_parser.add_argument(
        "--output", help="Output sequence file", required=True
    )
    codon_alignment_parser.add_argument(
        "--file-format", help="Sequence file format", default="fasta"
    )

    # Return the arguments
    return vars(codon_alignment_parser.parse_args())


def main():
    # Parse the arguments
    codon_args = codonAlignmentParser()

    # Start the logger
    startLogger()

    # Log the arguments
    logArgDict(codon_args)

    # Translate the sequence file
    mapGapsFromIndex(
        codon_args["unaligned_cds"],
        codon_args["aligned_aa"],
        codon_args["output"],
        codon_args["file_format"],
        confirm_match=True,
    )


if __name__ == "__main__":
    main()
