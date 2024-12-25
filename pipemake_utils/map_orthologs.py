import os
import argparse

import pandas as pd

from Bio import SeqIO

from misc import confirmDir, confirmFile


def mapParser():
    map_parser = argparse.ArgumentParser(description="Map orthologs onto MSFs")
    map_parser.add_argument(
        "--input-dir",
        help="Directory containing input sequence files",
        type=str,
        required=True,
        action=confirmDir(),
    )
    map_parser.add_argument(
        "--input-species", help="Name of input species column", type=str, required=True
    )
    map_parser.add_argument(
        "--map-file",
        help="File containing the sequences to map",
        type=str,
        required=True,
        action=confirmFile(),
    )
    map_parser.add_argument(
        "--map-species", help="Name of species column to map", type=str, required=True
    )
    map_parser.add_argument(
        "--orthogroups-file",
        help="Orthogroups file from OrthoFinder",
        type=str,
        required=True,
        action=confirmFile(),
    )
    map_parser.add_argument(
        "--out-dir", help="Output directory", type=str, default="Mapped_MSFs"
    )

    return vars(map_parser.parse_args())


def ifSGO(orthogroup):
    # Check if the series contains a NaN value
    if orthogroup.isnull().values.any():
        return False

    # Check if the series contains a comma
    if orthogroup.str.contains(",").any():
        return False

    return True


def main():
    # Read in the arguments
    map_args = mapParser()

    # Read in the orthogroups file
    orthogroups = pd.read_csv(map_args["orthogroups_file"], sep="\t")

    # Drop the Orthogroup column
    orthogroups.drop("Orthogroup", axis=1, inplace=True)

    # Filter the orthogroups to only include single copy orthologs
    orthogroups = orthogroups[orthogroups.apply(ifSGO, axis=1)]

    # Set the input species column as the index
    orthogroups.set_index(map_args["input_species"], inplace=True)

    # Create the orthogroup dictionary
    orthogroup_dict = orthogroups[map_args["map_species"]].to_dict()

    # Index the map fasta file
    map_index = SeqIO.index(map_args["map_file"], "fasta")

    # Loop through the orthogroups
    for input_fasta in os.listdir(map_args["input_dir"]):
        # Assign the input id
        input_id = input_fasta.rsplit(".", 1)[0]

        # Skip if the input id is not in the orthogroup dictionary
        if input_id not in orthogroup_dict:
            continue

        # Create the output directory if it does not exist
        if not os.path.exists(map_args["out_dir"]):
            os.makedirs(map_args["out_dir"])

        # Assign the output path
        output_path = os.path.join(map_args["out_dir"], f"{input_id}.fasta")

        # Open the output file
        with open(output_path, "w") as output_file:
            # Parse the input fasta file
            for record in SeqIO.parse(
                os.path.join(map_args["input_dir"], input_fasta), "fasta"
            ):
                # Write the record to the output file
                SeqIO.write(record, output_file, "fasta")

            # Save the orthogroup record
            mapped_record = map_index[orthogroup_dict[input_id]]
            mapped_record.id = "outgroup"
            mapped_record.description = ""
            SeqIO.write(mapped_record, output_file, "fasta")


if __name__ == "__main__":
    main()
