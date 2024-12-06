import argparse

from pipemake_utils.misc import confirmFile
from pipemake_utils.align import translateCDS
from pipemake_utils.logger import startLogger, logArgDict


def translateSeqParser():
    translate_parser = argparse.ArgumentParser(description="Translate a sequence file")
    translate_parser.add_argument(
        "--input", help="Input sequence file", required=True, action=confirmFile()
    )
    translate_parser.add_argument(
        "--output", help="Output sequence file", required=True
    )
    translate_parser.add_argument(
        "--file-format", help="Sequence file format", default="fasta"
    )

    return vars(translate_parser.parse_args())


def main():
    # Parse the arguments
    translate_args = translateSeqParser()

    # Start the logger
    startLogger()

    # Log the arguments
    logArgDict(translate_args)

    # Translate the sequence file
    translateCDS(
        translate_args["input"], translate_args["output"], translate_args["file_format"]
    )


if __name__ == "__main__":
    main()
