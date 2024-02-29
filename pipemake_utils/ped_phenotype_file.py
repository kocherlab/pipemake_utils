import argparse

import pandas as pd

from pipemake_utils.misc import *
from pipemake_utils.logger import *
from pipemake_utils.model import readModelFile

def argParser ():

    # Create argument parser
    parser = argparse.ArgumentParser(description = 'Update plink fam file with model category information')
    parser.add_argument('--fam', help = 'The plink table', type = str, action = confirmFile(), required = True)
    parser.add_argument('--model', help = 'The model file', type = str, action = confirmFile(), required = True)
    parser.add_argument('--model-category', help = 'The category to assign from the model file', type = str, required = True)
    parser.add_argument('--categories-limit', help = 'The total number of categories allowed', type = int, default = 2)
    parser.add_argument('--out-prefix', help = 'The output prefix', type = str, default = 'out')

    # Parse the arguments
    return vars(parser.parse_args())

def mapCategory (row, category_map):

    # Assign the individual
    row_ind = row[1]

    # Check if the individual is in the map
    if row_ind not in category_map: raise Exception(f"Individual ({row_ind}) not found in the category map")
    
    # Map the category
    row_category = category_map[row_ind]

    # Update the row
    row[5] = row_category

    logging.info(f"Mapped individual: {row_ind} to {row_category}")

    return row

def main():

    # Parse the arguments
    map_args = argParser()

    # Start logger and log the arguments
    startLogger(f"{map_args['out_prefix']}.pheno.log")
    logArgDict(map_args)

    # Assign the category
    models = readModelFile(map_args['model'])
    model_category = models[map_args['model_category']]

    # Store the categorys
    categories = list(model_category.ind_dict)

    # Create a dictionary to map the individual to the category
    ind_to_categories = {}

    # Map the individuals to the categories and log the information
    for category_idx, (category_str, inds) in enumerate(model_category.ind_dict.items(), 1):
        logging.info(f"Mapping: {category_str} to {category_idx}")
        for ind in inds: 
            ind_to_categories[ind] = category_idx
           
    # Read the plink table
    plink_table = pd.read_csv(map_args['fam'], sep = ' ', header = None)

    # Map the categories
    plink_table = plink_table.apply(mapCategory, category_map = ind_to_categories, axis = 1)

    # Confirm that the categories were mapped
    if len(plink_table[5].unique()) != map_args['categories_limit']:
        raise Exception(f"Categories ({', '.join(categories)}) greater than the limit: {map_args['categories_limit']}")

    # Save the updated plink table
    plink_table[5].to_csv(f"{map_args['out_prefix']}.pheno.txt", sep = ' ', header = False, index = False)

if __name__ == '__main__':
	main()
