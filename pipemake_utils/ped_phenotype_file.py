#!/usr/bin/env python

import logging
import argparse

import numpy as np
import pandas as pd

from pipemake_utils.misc import confirmFile
from pipemake_utils.logger import startLogger, logArgDict
from pipemake_utils.model import ModelFile


def argParser():
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="Update plink fam file with model category information"
    )
    parser.add_argument(
        "--fam", help="The plink table", type=str, action=confirmFile(), required=True
    )
    parser.add_argument(
        "--model-file",
        help="The model file",
        type=str,
        action=confirmFile(),
        required=True,
    )
    parser.add_argument(
        "--model-name",
        help="The model to assign from the model file",
        type=str,
        required=True,
    )
    parser.add_argument(
        "--out-prefix", help="The output prefix", type=str, default="out"
    )
    parser.add_argument(
        "--out-format",
        help="The output format",
        type=str,
        choices=["plink2", "basic"],
        default="basic",
    )
    parser.add_argument(
        "--pheno-header", help="The header for the phenotype column", type=str
    )
    parser.add_argument(
        "--numeric",
        help="Report phenotypes numerically (1, 2, etc.)",
        action="store_true",
    )
    parser.add_argument(
        "--binary", help="Report phenotypes in binary (0, 1)", action="store_true"
    )
    parser.add_argument(
        "--phenotype-limit",
        help="The total number of phenotypes allowed",
        type=int,
        default=2,
    )

    # Parse the arguments
    return vars(parser.parse_args())


def mapPhenotype(row, phenotype_map):
    # Assign the individual
    row_ind = row[1]

    # Check if the individual is in the map
    if row_ind not in phenotype_map:
        # Assign the phenotype as missing
        row[5] = np.nan

        logging.warn(f"Individual ({row_ind}) not found in the phenotype map")

    else:
        # Map the model
        row_phenotype = phenotype_map[row_ind]
        row[5] = row_phenotype

        logging.info(f"Mapped individual: {row_ind} to {row_phenotype}")

    return row


def main():
    # Parse the arguments
    map_args = argParser()

    # Start logger and log the arguments
    startLogger(f"{map_args['out_prefix']}.pheno.log")
    logArgDict(map_args)

    # Assign the category
    models = ModelFile.read(map_args["model_file"])
    model_name = models[map_args["model_name"]]

    # Check if the user specified both binary and numeric
    if map_args["binary"] and map_args["numeric"]:
        raise Exception("Cannot specify both binary and numeric")

    # Check if the user specified neither binary with a limit greater than 2
    if not map_args["binary"] and map_args["phenotype_limit"] > 2:
        raise Exception(
            f"Too many phenotypes for binary mapping: {map_args['phenotype_limit']}"
        )

    # Create a dictionary to map the individual to the phenotype
    ind_to_phenotypes = {}

    # Map the individuals to the phenotypes and log the information
    for phenotype_idx, (phenotype_str, inds) in enumerate(model_name.ind_dict.items()):
        # Assign the phenotype id
        phenotype_id = phenotype_str

        # Check if the phenotypes should be reported numerically or in binary
        if map_args["binary"] or map_args["numeric"]:
            # Check if the phenotypes should be reported in binary
            if map_args["binary"]:
                phenotype_id = phenotype_idx
            elif map_args["numeric"]:
                phenotype_id = phenotype_idx + 1
            logging.info(f"Mapping: {phenotype_str} to {phenotype_idx}")

        for ind in inds:
            ind_to_phenotypes[ind] = phenotype_id

    # Read the plink table
    plink_table = pd.read_csv(map_args["fam"], sep=" ", header=None)

    # Check if the plink table only has one column, and if so, read it in again with tab delimiters
    if len(plink_table.columns) == 1:
        plink_table = pd.read_csv(map_args["fam"], sep="\t", header=None)

    # Map the phenotypes
    plink_table = plink_table.apply(
        mapPhenotype, phenotype_map=ind_to_phenotypes, axis=1
    )

    # Check if the phenotypes were mapped
    if len(plink_table[5].unique()) == 1 and np.isnan(plink_table[5].unique()[0]):
        raise Exception("No phenotypes were mapped")

    # Check if the user specified basic format
    if map_args["out_format"] == "basic":
        # Report an error if the user specified a phenotype header
        if map_args["pheno_header"]:
            raise Exception("Cannot specify a phenotype header with basic format")

        # Write the file
        plink_table[5].to_csv(
            f"{map_args['out_prefix']}.pheno.txt", sep=" ", header=False, index=False
        )

    # Check if the user specified plink2 format
    elif map_args["out_format"] == "plink2":
        # Check if the user specified a phenotype header and rename the column
        if map_args["pheno_header"]:
            plink_table = plink_table.rename(columns={5: map_args["pheno_header"]})

        # Rename the column and write the file
        plink_table = plink_table.rename(columns={1: "#IID"})
        plink_table[["#IID", map_args["pheno_header"]]].to_csv(
            f"{map_args['out_prefix']}.pheno.txt",
            sep="\t",
            header=True,
            index=False,
            na_rep="NF",
        )

    else:
        raise Exception(f"Unknown format: {map_args['out_format']}")


if __name__ == "__main__":
    main()
