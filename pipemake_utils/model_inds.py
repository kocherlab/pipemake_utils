import argparse

import pandas as pd

from pipemake_utils.misc import *
from pipemake_utils.logger import *
from pipemake_utils.model import readModelFile

def argParser ():

    # Create argument parser
    parser = argparse.ArgumentParser(description = 'Create a file of the individuals from a model')
    parser.add_argument('--model-file', help = 'The model file', type = str, action = confirmFile(), required = True)
    parser.add_argument('--model-name', help = 'The name to assign from the model file', type = str, required = True)
    output = parser.add_mutually_exclusive_group(required = True)
    output.add_argument('--out-prefix', help = 'The output prefix', type = str)
    output.add_argument('--out-dir', help = 'The output directory. The output filename ', type = str)

    # Parse the arguments
    return vars(parser.parse_args())

def main():

    # Parse the arguments
    map_args = argParser()

    # Start logger and log the arguments
    startLogger(f"{map_args['out_prefix']}.ind.log")
    logArgDict(map_args)

    # Assign the model
    models = readModelFile(map_args['model'])
    model = models[map_args['model_name']]

    # Assign the output filename
    if map_args['out_prefix']: out_filename = f"{map_args['out_prefix']}.ind.txt"
    else: out_filename = os.path.join(map_args['out_dir'], f"{map_args['model_name']}.ind.txt")

    # Print the individuals to a file
    with open(out_filename, 'w') as ind_file:
    
        # Loop through the model and write the individuals to a file
        for inds in model.ind_dict.items():
            for ind in inds:
                nd_file.write(f"{ind}\t{ind}\n")
                logging.info(f"Added {ind} to {map_args['out_prefix']}.ind.txt")

if __name__ == '__main__':
	main()
