import os
import argparse

from pipemake_utils.misc import confirmDir
from pipemake_utils.align import confirmCDS
from pipemake_utils.logger import startLogger, logArgDict


def confirmCDSParser():
    confirm_parser = argparse.ArgumentParser(
        description="Confirm the CDS sequences in a file"
    )
    confirm_parser.add_argument(
        "--input-dir",
        help="Input sequence directory",
        required=True,
        action=confirmDir(),
    )
    confirm_parser.add_argument(
        "--output-dir", help="Output sequence directory", required=True
    )
    confirm_parser.add_argument(
        "--file-format", help="Sequence file format", default="fasta"
    )

    return vars(confirm_parser.parse_args())


def main():
    # Parse the arguments
    confirm_args = confirmCDSParser()

    # Start the logger
    startLogger()
    logArgDict(confirm_args)

    # Create the output directory
    if not os.path.exists(confirm_args["output_dir"]):
        os.makedirs(confirm_args["output_dir"])

    # Loop through the files in the input directory
    for input_filename in os.listdir(confirm_args["input_dir"]):
        # Assign the input path
        input_path = os.path.join(confirm_args["input_dir"], input_filename)

        # Assign the output path
        output_path = os.path.join(confirm_args["output_dir"], input_filename)

        # Confirm the CDS sequences in the file
        confirmCDS(input_path, output_path, confirm_args["file_format"])


if __name__ == "__main__":
    main()
