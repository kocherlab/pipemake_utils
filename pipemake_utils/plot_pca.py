import argparse

import seaborn as sns
import pandas as pd

from pipemake_utils.misc import *
from pipemake_utils.logger import *
from pipemake_utils.model import readModelFile

def argParser ():

	# Create argument parser
	parser = argparse.ArgumentParser(description = 'Create PCA plot using PLINK PCA files and a model file')
	parser.add_argument('--pca-dir', help = 'The PCA directory', type = str, action = confirmDir(), required = True)
	parser.add_argument('--model-file', help = 'The model file', type = str, action = confirmFile(), required = True)
	parser.add_argument('--model-name', help = 'The name to assign from the model file', type = str, required = True)
	parser.add_argument('--out-prefix', help = 'The output prefix', type = str, default = 'out')
	
	# Parse the arguments
	return vars(parser.parse_args())

class PCA ():
	def __init__ (self, eigenvectors_dataframe, eigenvalues_dataframe, **kwargs):

		# Confirm the eigenvectors dataframe has the expected columns
		if not set(['IID', 'PC1', 'PC2']).issubset(set(eigenvectors_dataframe.columns)):
			raise Exception(f'Error parsing eigenvectors columns: {eigenvectors_dataframe.columns}')

		# Confirm the eigenvalues dataframe has the expected column
		if not set(['PC1', 'PC2']).issubset(set(eigenvalues_dataframe.columns)):
			raise Exception(f'Error parsing eigenvalues columns: {eigenvalues_dataframe.columns}')

		if 'Eigenvalues' not in eigenvalues_dataframe.index:
			raise Exception(f'Error parsing eigenvalues index: {eigenvalues_dataframe.index}')

		# Store the dataframes
		self._eigenvectors_dataframe = eigenvectors_dataframe
		self._eigenvalues_dataframe = eigenvalues_dataframe

		logging.info(f'Successfully created PCA object from PLINK directory')

	@classmethod
	def fromPlinkDir (cls, directory, **kwargs):

		# Create strings to store PCA filenames
		eigenvectors_filename = ''
		eigenvalues_filename = ''
			
		# Look for the files
		for pca_file in os.listdir(directory):
			if pca_file.endswith('.eigenvec'):
				if eigenvectors_filename: raise Exception(f'PLINK eigenvectors already assigned')
				eigenvectors_filename = os.path.join(directory, pca_file)
			elif pca_file.endswith('.eigenval'):
				if eigenvalues_filename: raise Exception(f'PLINK eigenvalues already assigned')
				eigenvalues_filename = os.path.join(directory, pca_file)

		# Confirm the files were found
		if not eigenvectors_filename or not eigenvalues_filename: raise Exception(f'Unable to load PLINK PCA')

		# Create dataframes of the files
		eigenvectors_dataframe = pd.read_csv(eigenvectors_filename, sep = '\t')
		eigenvalues_dataframe = pd.read_csv(eigenvalues_filename, sep = '\t', header = None, names = ['Eigenvalues'])
		eigenvalues_dataframe.index = [f'PC{_i + 1}' for _i in eigenvalues_dataframe.index]

		return cls(eigenvectors_dataframe, eigenvalues_dataframe.T, **kwargs)

	def plotUsingModel (self, model_file, model_name, out_prefix, **kwargs):

		# Assign the Model
		models = readModelFile(model_file)
		model = models[model_name]

		# Store the ind to category
		ind_to_category = {}
		for category, inds in model.ind_dict.items():
			for ind in inds: ind_to_category[ind] = category.capitalize()

		# Update the dataframe with the category
		self._eigenvectors_dataframe[model] = self._eigenvectors_dataframe['IID'].map(ind_to_category)

		logging.info(f'Successfully assigned model to PCA dataframe')

		# Create the PCA plot
		pca_plot = sns.scatterplot(data = self._eigenvectors_dataframe, x = 'PC1', y = 'PC2', hue = model, alpha = 0.75, linewidth = 0.0)

		# Update the axis and plot titles
		pca_plot.set(xlabel = f"PC1 ({self._eigenvalues_dataframe.loc['Eigenvalues']['PC1']}%)", 
					 ylabel = f"PC2 ({self._eigenvalues_dataframe.loc['Eigenvalues']['PC2']}%)", 
					 title = f'{model} PCA')
		
		# Save the figure
		pca_fig = pca_plot.get_figure()
		pca_fig.savefig(f'{out_prefix}.pdf')

		logging.info(f'Successfully created PCA plot: {out_prefix}.pdf')

def main():

	# Parse the arguments
	pca_args = argParser()

	# Start logger and log the arguments
	startLogger(f"{pca_args['out_prefix']}.plot.log")
	logArgDict(pca_args)

	# Create the PCA object using the PLINK directory
	pca = PCA.fromPlinkDir(pca_args['pca_dir'])

	# Plot the PCA using the model file and model name
	pca.plotUsingModel(**pca_args)

if __name__ == '__main__':
	main()
	