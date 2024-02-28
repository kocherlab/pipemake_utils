import argparse

import pandas as pd

from kp_tools.misc import *
from kp_tools.logger import *
from kp_tools.model import readModelFile

def argParser ():

    # Create argument parser
    parser = argparse.ArgumentParser(description = 'Create a file of the individuals from a model category')
    parser.add_argument('--model', help = 'The model file', type = str, action = confirmFile(), required = True)
    parser.add_argument('--model-category', help = 'The category to assign from the model file', type = str, required = True)
    parser.add_argument('--out-prefix', help = 'The output prefix', type = str, default = 'out')

    # Parse the arguments
    return vars(parser.parse_args())

def main():

    # Parse the arguments
    map_args = argParser()

    # Start logger and log the arguments
    startLogger(f"{map_args['out_prefix']}.log")
    logArgDict(map_args)

    # Assign the category
    models = readModelFile(map_args['model'])
    model_category = models[map_args['model_category']]

    # Print the individuals to a file
    with open(f"{map_args['out_prefix']}.ind", 'w') as ind_file:
        for inds in model_category.ind_dict.values():
            for ind in inds:
                ind_file.write(f"{ind}\n")
                logging.info(f"Added {ind} to {map_args['out_prefix']}.ind")

if __name__ == '__main__':
	main()
