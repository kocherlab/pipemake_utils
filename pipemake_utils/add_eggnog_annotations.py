#!/usr/bin/env python

import argparse
import gffutils

import pandas as pd

from Bio import SeqIO

from pipemake_utils.misc import confirmFile
from pipemake_utils.seqDB import DBFileReader
from pipemake_utils.logger import startLogger, logArgDict


def eggnogParser():
    eggnog_parser = argparse.ArgumentParser(
        description="Annotate sequence files with gene names from eggNOG output"
    )
    eggnog_parser.add_argument(
        "--eggnog-excel",
        help="eggNOG output file in Excel format",
        type=str,
        required=True,
        action=confirmFile(),
    )
    eggnog_parser.add_argument(
        "--gff", help="GFF3 file to be annotated", type=str, action=confirmFile()
    )
    eggnog_parser.add_argument(
        "--fasta-cds",
        help="CDS fasta file to be annotated",
        type=str,
        action=confirmFile(),
    )
    eggnog_parser.add_argument(
        "--fasta-aa",
        help="AA fasta file to be annotated",
        type=str,
        action=confirmFile(),
    )
    eggnog_parser.add_argument(
        "--output-prefix", help="Output prefix", type=str, default="out"
    )

    return vars(eggnog_parser.parse_args())


def eggnogGeneName(eggnog_filename):
    # Read in the eggNOG output file
    eggnog_dataframe = pd.read_excel(eggnog_filename, comment="#")

    # Remove columns that have None in the query column (comment lines at the end of the file)
    eggnog_dataframe = eggnog_dataframe[eggnog_dataframe["query"].notnull()]

    # Remove columns that have a dash in the Preferred_name column (i.e. no annotation)
    eggnog_dataframe = eggnog_dataframe[
        ~eggnog_dataframe["Preferred_name"].str.contains("-")
    ]

    # Create a dictionary of the eggNOG output file
    return eggnog_dataframe.set_index("query")["Preferred_name"].to_dict()


def createGff3(gff3_filename, output_prefix, eggnog_dict):
    with open(f"{output_prefix}.gff3", "w") as gff_output:
        # Write the GFF3 header
        gff_output.write("##gff-version 3\n")

        # Iterate over the GFF3 file
        for feature in gffutils.DataIterator(gff3_filename):
            # Skip the feature if it is not a gene
            if feature.featuretype == "gene":
                # Confirm a single feature ID was found, otherwise report an error
                if len(feature.attributes["ID"]) > 1:
                    raise ValueError("Multiple feature IDs found for a single feature")

                # Assign the feature ID to a variable
                feature_id = feature.attributes["ID"][0]

                # Assign the preferred name to the Name attribute using get
                feature.attributes["Name"] = eggnog_dict.get(feature_id, [feature_id])

            # Print the feature
            gff_output.write(f"{feature}\n")


def createFasta(fasta_filename, output_prefix, file_type, eggnog_dict):
    with open(f"{output_prefix}.fasta", "w") as fasta_output:
        # Iterate over the fasta file
        for record in DBFileReader.read(
            "pipemake",
            fasta_filename,
            file_type=file_type,
            file_format="fasta",
            primary_id="protein_id",
        ):
            # Assign the gene name to the record
            gene_name = eggnog_dict.get(record.gene_id, record.gene_id)

            # Assign the transcript number
            transcript_number = record.transcript_id.split("-")[1]

            # Add the gene name to th
            record.addAttribute("protein", f"{gene_name} isoform {transcript_number}")

            # Write the record to the output file
            SeqIO.write(record, fasta_output, "fasta")


def main():
    # Parse the arguments
    eggnog_args = eggnogParser()

    # Confirm at least one file to be annotated was provided
    if not any([eggnog_args["gff"], eggnog_args["fasta_cds"], eggnog_args["fasta_aa"]]):
        raise ValueError("At least one file to be annotated must be provided")

    # Start logger and log the arguments
    startLogger(f"{eggnog_args['output_prefix']}.log")
    logArgDict(eggnog_args)

    # Create a dictionary of the eggNOG output file
    eggnog_dict = eggnogGeneName(eggnog_args["eggnog_excel"])

    # Create the GFF3 file, if provided
    if eggnog_args["gff"]:
        createGff3(eggnog_args["gff"], eggnog_args["output_prefix"], eggnog_dict)

    # Create the CDS fasta file, if provided
    if eggnog_args["fasta_cds"]:
        createFasta(
            eggnog_args["fasta_cds"],
            f"{eggnog_args['output_prefix']}_trans",
            "CDS",
            eggnog_dict,
        )

    # Create the AA fasta file, if provided
    if eggnog_args["fasta_aa"]:
        createFasta(
            eggnog_args["fasta_aa"],
            f"{eggnog_args['output_prefix']}_pep",
            "AA",
            eggnog_dict,
        )


if __name__ == "__main__":
    main()
