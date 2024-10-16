#!/usr/bin/env python

import argparse

import pandas as pd
import scipy.stats as sstats

from pipemake_utils.logger import *
from pipemake_utils.misc import *

def argParser ():

	parser = argparse.ArgumentParser(description = 'Z normalize a column in a file')

	parser.add_argument('--input-file', help = 'Input filename', type = str, action = confirmFile(), required = True)
	parser.add_argument('--input-format', help = 'Input file format', type = str, default = 'tsv', choices = ['tsv', 'csv'])
	parser.add_argument('--normalize-col', help = 'Column to normalize', type = str, required = True)
	parser.add_argument('--out-prefix', help = 'Output prefix', type = str, default = 'out')

	# Add mutually exclusive arguments for output options
	option_parser = parser.add_mutually_exclusive_group()
	option_parser.add_argument('--abs', help = 'Report ABS z-scores', action = 'store_true')
	option_parser.add_argument('--negative-to-zero', help = 'Convert negative z-scores to zero', action = 'store_true')

	return vars(parser.parse_args())

def zNormalize (array_col):
	return sstats.zscore(array_col, nan_policy = 'omit')

def main():

	# Parse the arguments
	normalize_args = argParser()

	# Assign the output filename
	normalize_args['out_filename'] = f"{normalize_args['out_prefix']}.{normalize_args['input_format']}"

	# Start logger and log the arguments
	startLogger(f"{normalize_args['out_prefix']}.normalize.log")
	logArgDict(normalize_args)

	# Assign the input file separator
	if normalize_args['input_format'] == 'tsv': input_separator = '\t'
	elif normalize_args['input_format'] == 'csv': input_separator = ','
	else: raise Exception (f'Unknown input format: {normalize_args["input_format"]}')

	# Read in the table as a pandas dataframe
	normalized_table = pd.read_csv(normalize_args['input_file'], sep = input_separator)

	# Assign the name of the column to be z-normalized
	z_normalized_col = f"Z({normalize_args['normalize_col']})"

	# Z-normalize the column
	normalized_table[z_normalized_col] = zNormalize(normalized_table[normalize_args['normalize_col']])

	# Check if the user wants to convert negative z-scores to zero
	if normalize_args['negative_to_zero']: normalized_table.loc[normalized_table[z_normalized_col] < 0, z_normalized_col] = 0
  	
	# Check if the user wants to report the absolute z-scores
	if normalize_args['abs']: normalized_table[z_normalized_col] = abs(normalized_table[z_normalized_col])

	# Write the z-normalized file to a new file
	normalized_table.to_csv(normalize_args['out_filename'], sep = input_separator, index = False)

if __name__ == '__main__':
	main()