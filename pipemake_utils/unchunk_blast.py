import csv
import argparse

from pipemake_utils.logger import startLogger, logArgDict


def unchunkArgs():
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--chunked-blast",
        help="Input file containing chunked BLAST results.",
        required=True,
    )
    parser.add_argument(
        "--unchunked-blast",
        help="Output file for unchunked BLAST results.",
        default="unchunked_blast.out",
    )
    return vars(parser.parse_args())


def unchunkBlast(input_filename, output_filename):
    with open(input_filename, "r") as input_file, open(
        output_filename, "w"
    ) as output_file:
        input_reader = csv.reader(input_file, delimiter="\t")
        for input_row in input_reader:
            # Process the chunked header
            input_header, input_range = input_row[0].split(":")
            input_start, _ = map(int, input_range.split("-"))

            # Update the input_row with the unchunked header and adjusted positions
            input_row[0] = input_header
            input_row[3] = input_header
            input_row[9] = int(input_row[9]) + input_start - 1
            input_row[10] = int(input_row[10]) + input_start - 1

            # Write the modified row to the output file
            output_file.write("\t".join(map(str, input_row)) + "\n")


def main():
    args = unchunkArgs()

    # Start logger and log the arguments
    startLogger("unchunked_blast.log")
    logArgDict(args)

    unchunkBlast(args["chunked_blast"], args["unchunked_blast"])


if __name__ == "__main__":
    main()
