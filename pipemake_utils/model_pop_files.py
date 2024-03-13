import argparse

from pipemake_utils.misc import *
from pipemake_utils.logger import *
from pipemake_utils.model import readModelFile


def argParser ():

	# Create argument parser
	parser = argparse.ArgumentParser(description = 'Create population files from a model file')
	parser.add_argument('--model-file', help = 'Model filename', type = str, action = confirmFile(), required = True)
	parser.add_argument('--model-name', help = 'Model name', type = str, required = True)
	parser.add_argument('--out-dir', help = 'Output directory', type = str, required = True)

	# Parse the arguments
	return vars(parser.parse_args())

def main():

	# Parse the arguments
	model_args = argParser()

	# Read the model file
	models = readModelFile(model_args['model_file'])

	# Create the population files
	models[model_args['model_name']].create_pop_files(file_ext = 'pop', file_path = model_args['out_dir'])

if __name__ == '__main__':
	main()