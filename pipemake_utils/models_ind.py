#!/usr/bin/env python

import os
import logging
import argparse

from pipemake_utils.misc import confirmFile
from pipemake_utils.logger import startLogger, logArgDict
from pipemake_utils.model import ModelFile


def argParser():
    # Create argument parser
    parser = argparse.ArgumentParser(
        description="Create a file of the individuals from a model"
    )
    parser.add_argument(
        "--model-file",
        help="The model file",
        type=str,
        action=confirmFile(),
        required=True,
    )
    parser.add_argument(
        "--models",
        help="The name to assign from the model file",
        type=str,
        nargs="+",
        required=True,
    )
    output = parser.add_mutually_exclusive_group(required=True)
    output.add_argument("--out-prefix", help="The output prefix", type=str)
    output.add_argument(
        "--out-dir", help="The output directory. The output filename ", type=str
    )

    # Parse the arguments
    return vars(parser.parse_args())


def main():
    # Parse the arguments
    map_args = argParser()

    # Start logger and log the arguments
    startLogger(f"{map_args['out_prefix']}.ind.log")
    logArgDict(map_args)

    # Assign the output filename
    if map_args["out_prefix"]:
        out_filename = f"{map_args['out_prefix']}.ind.txt"
    else:
        out_filename = os.path.join(
            map_args["out_dir"], f"{map_args['model_name']}.ind.txt"
        )

    # Assign the model
    models_file = ModelFile.read(map_args["model_file"])

    # Create list to store the inds
    inds_list = []

    # Loop through models
    for model_name in map_args["models"]:
        if model_name not in models_file:
            raise ValueError(f"Model name {model_name} not found in {map_args['model_file']}")
        
        model = models_file[model_name]

        if not inds_list:
            inds_list = sorted(_ind for _inds in model.ind_dict.values() for _ind in _inds)
        
        else:
            chk_list = sorted(_ind for _inds in model.ind_dict.values() for _ind in _inds)

            if inds_list != chk_list:
                raise ValueError("Individuals do not match between models")

    # Print the individuals to a file
    with open(out_filename, "w") as ind_file:
        for ind in inds_list:
            ind_file.write(f"{ind}\n")
            logging.info(f"Added {ind} to {map_args['out_prefix']}.ind.txt")


if __name__ == "__main__":
    main()
