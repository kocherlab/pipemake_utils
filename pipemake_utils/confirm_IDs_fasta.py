#!/usr/bin/env python3
import logging
import argparse

from Bio import SeqIO

from pipemake_utils.misc import confirmFile
from pipemake_utils.logger import startLogger, logArgDict


def parseIDArgs ():
    parser = argparse.ArgumentParser(description="Find which ID prefixes appear in a FASTA file's headers.")

    parser.add_argument("--fasta-file", help = "Path to the input FASTA file.", required = True, action = confirmFile())
    prefix_group = parser.add_mutually_exclusive_group(required = True)
    prefix_group.add_argument("--ids-file", help = "File containing one ID prefix per line.", action = confirmFile())
    prefix_group.add_argument("--ids", nargs = "*", default = [], help = "One or more ID prefixes given directly on the command line.")
    parser.add_argument("--output", help = "Output file for matched prefixes (default: matched_prefixes.txt)", default = "matched_ids.txt")
    parser.add_argument("--ignore-case", help = "Match prefixes case-insensitively.", action = "store_true")
    parser.add_argument("--first-only", help = "Report only the first matching ID.", action = "store_true")
    parser.add_argument("--prefix-ids", help = "Treat the provided IDs as prefixes to match against the FASTA headers.", action = "store_true")

    return vars(parser.parse_args())

def loadIDs (path):
    ids = []
    with open(path, "r") as id_file:
        for id_line in id_file:
            if not id_line.strip() or id_line.startswith("#"):
                continue
            ids.append(id_line.strip())

    logging.info(f"Loaded {len(ids)} IDs from {path}.")

    return ids

def find_matches(fasta_path, id_list, ignore_case=False, first_only=False, prefix_ids=False, **kwargs):

    # Create a set to store matched
    matches = set()

    logging.info(f"Searching for {len(id_list)} ID(s) in {fasta_path}.")

    for record in SeqIO.parse(fasta_path, "fasta"):

        # If the user wants to ignore case, convert the header to lowercase
        header = record.id.lower() if ignore_case else record.id

        # If the user wants to treat the IDs as prefixes, check if the header starts with the ID
        for to_match in id_list:

            # If the user wants to ignore case, convert the ID to lowercase
            to_match = to_match.lower() if ignore_case else to_match
            
            if prefix_ids and header.startswith(to_match):
                matches.add(header)
                logging.debug(f"Matched prefix: {to_match} in header: {header}")

            elif to_match == header:
                matches.add(header)
                logging.debug(f"Matched ID: {to_match} in header: {header}")

            # Check for the first_only flag to break out of the loop if a match was found
            if len(matches) > 0 and first_only:
                logging.info(f"Matched {len(matches)} ID(s) in {fasta_path}.")
                return list(matches)        

    logging.info(f"Matched {len(matches)} ID(s) in {fasta_path}.")
    return list(matches)

def main():

    # Parse command-line arguments
    id_args = parseIDArgs()

    # Start the logger
    startLogger()

    # Log the arguments
    logArgDict(id_args)

    # Load the IDs from the file or command line arguments
    if id_args["ids_file"]:
        id_list = loadIDs(id_args["ids_file"])
    else:
        id_list = [p.strip() for p in id_args["ids"] if p.strip()]

    # Assign the matches
    matches = find_matches(id_args["fasta_file"], id_list, **id_args)
    
    # Report if no matches were found
    if not matches:
        raise SystemExit("No matches found.")

    logging.info(f"Writing {len(matches)} matches to {id_args['output']}.")

    # Write the results to the output file
    with open(id_args["output"], "w") as matched_output:
        for match in matches:
            matched_output.write(match + "\n")

    logging.info(f"Finished writing matches to {id_args['output']}.")

if __name__ == "__main__":
    main()