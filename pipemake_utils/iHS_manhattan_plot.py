#!/usr/bin/env python

import argparse

import pandas as pd
import matplotlib.pyplot as plt

from pipemake_utils.misc import *
from pipemake_utils.logger import *

def argParser ():

	parser = argparse.ArgumentParser(description = 'Filter normalized iHS output')

	parser.add_argument('--iHS-file', help = 'iHS output file', type = str, action = confirmFile(), required = True)
	parser.add_argument('--out-prefix', help = 'Output prefix', type = str, default = 'out')
	parser.add_argument('--iHS-col', help = 'iHS column name', type = str, default = 'Normalized iHS')
	parser.add_argument('--chrom-col', help = 'CRHOM column name', type = str, default = 'CHROM')
	parser.add_argument('--pos-col', help = 'POS column name', type = str, default = 'POS')
	parser.add_argument('--plot-abs', help = 'Plot the absolute value instead', action='store_true')
	parser.add_argument('--plot-dpi', help = 'Plot DPI', type = int, default = 100)

	return vars(parser.parse_args())

def main():

	# Parse the arguments
	plot_args = argParser()

	# Start logger and log the arguments
	startLogger(f"{plot_args['out_prefix']}.plot.log")
	logArgDict(plot_args)

	# Read in the PBS data
	plot_dataframe = pd.read_csv(plot_args['iHS_file'], sep = '\t', header = None)

	# Assign iHS columns
	plot_dataframe[[plot_args['chrom_col'], plot_args['pos_col']]] = plot_dataframe[0].str.rsplit('_', n = 1, expand = True)
	plot_dataframe[plot_args['pos_col']] = plot_dataframe[plot_args['pos_col']].astype(int)
	plot_dataframe = plot_dataframe.rename(columns={6: plot_args['iHS_col']})
	plot_dataframe = plot_dataframe[[plot_args['chrom_col'], plot_args['pos_col'], plot_args['iHS_col']]]

	# Plot the absolute value, if specified
	if plot_args['plot_abs']: plot_dataframe[plot_args['iHS_col']] = plot_dataframe[plot_args['iHS_col']].abs()

	# Assign the integer chromosome column name
	int_col_name = f"{plot_args['chrom_col']}_int"

	# Sort by CHROM and POS or CHROM and BIN_START
	plot_dataframe[int_col_name] = plot_dataframe[plot_args['chrom_col']].str.split('_').str[-1].astype(int)
	plot_dataframe = plot_dataframe.sort_values(by = [int_col_name, plot_args['pos_col']])

	# Get data length
	pbs_data_len = len(plot_dataframe.index)
	plot_dataframe['ORDER'] = range(pbs_data_len)

	# Group the dataframe
	grouped_plot_dataframe = plot_dataframe.groupby([int_col_name])

	# Create the plot
	plt.rcParams["figure.figsize"] = (24, 6)
	fig = plt.figure()
	ax = fig.add_subplot(111)
	colors = ['royalblue','navy']
	bg_colors = ['white','whitesmoke']
	x_labels = []
	x_labels_pos = []

	for num, (name, group) in enumerate(grouped_plot_dataframe):
		if group.empty: continue

		group.plot(kind='scatter', x='ORDER', y=plot_args['iHS_col'], color=colors[num % len(colors)], ax=ax, s = 1, zorder = 2)
		plt.axvspan(group['ORDER'].iloc[0], group['ORDER'].iloc[-1], facecolor = bg_colors[num % len(bg_colors)], zorder=1)

		x_labels.append(group[plot_args['chrom_col']].iloc[0])
		x_labels_pos.append((group['ORDER'].iloc[-1] - (group['ORDER'].iloc[-1] - group['ORDER'].iloc[0])/2))

	ax.set_xticks(x_labels_pos)
	ax.set_xticklabels(x_labels)
	ax.set_xlim([0, pbs_data_len])
	ax.set_ylim([plot_dataframe[plot_args['iHS_col']].min() * 1.5, plot_dataframe[plot_args['iHS_col']].max() * 1.5])
	ax.set_xlabel('Chromosome')

	plt.xticks(rotation=90, fontsize=10)
	plt.tight_layout()
	plt.savefig(f"{plot_args['out_prefix']}.manhattan.png", dpi = plot_args['plot_dpi'])	

if __name__ == '__main__':
	main()