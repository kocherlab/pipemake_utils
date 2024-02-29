import sys
import argparse

import pandas as pd
import numpy as np

from pipemake_utils.misc import *
from pipemake_utils.logger import *

def argParser ():

	parser = argparse.ArgumentParser(description = 'Filter GEMMA output')

	parser.add_argument('--gemma-file', help = 'GEMMA output file', type = str, action = confirmFile(), required = True)
	parser.add_argument('--min-log-pvalue', help = 'Minimum -log10(p-value) to keep', type = float, default = 0.5)
	parser.add_argument('--out-prefix', help = 'Output prefix', type = str, default = 'out')
	parser.add_argument('--pvalue-col', help = 'P-value column to filter', type = str, default = 'p_wald')
	parser.add_argument('--chrom-col', help = 'CRHOM column', type = str, default = 'chr')
	parser.add_argument('--pos-col', help = 'POS column', type = str, default = 'ps')

	return vars(parser.parse_args())

def main():

	# Parse the arguments
	filter_args = argParser()

	# Start logger and log the arguments
	startLogger(f"{filter_args['out_prefix']}.filtered.log")
	logArgDict(filter_args)

	# Open the GEMMA output as a dataframe
	gemma_dataframe = pd.read_csv(filter_args['gemma_file'], sep = '\t')

	# Assign the name of the -log10(p-value) columns
	log_y_col = f"-log10({filter_args['pvalue_col']})"

	# Create the -log10(p-value) columns using the p-values
	gemma_dataframe[log_y_col] = -np.log10(gemma_dataframe[filter_args['pvalue_col']])
	gemma_dataframe[log_y_col] = gemma_dataframe[log_y_col].fillna(0)

	# Store the number of rows before filtering
	pre_filter_values = gemma_dataframe.shape[0]

	logging.info(f"Filtering GEMMA column {log_y_col} to remove p-values less than {filter_args['min_log_pvalue']}")

	# Filter the GEMMA dataset
	gemma_dataframe = gemma_dataframe[(gemma_dataframe[log_y_col] >= filter_args['min_log_pvalue'])]

	logging.info(f"Filtering complete: Removed {pre_filter_values - gemma_dataframe.shape[0]} values ({gemma_dataframe.shape[0]} remaining)")

	# Save the filtered GEMMA dataset	
	gemma_dataframe.to_csv(f"{filter_args['out_prefix']}.filtered.assoc.txt", sep = '\t', index = False)

if __name__ == '__main__':
	main()