import argparse
import pandas as pd

def pveParser():
    pve_parser = argparse.ArgumentParser(description="Calculate PVE for each SNP.")
    pve_parser.add_argument("--input", dest = 'input_filename', help="Input filename", type=str, required=True)
    pve_parser.add_argument("--output", dest = 'output_filename', help="Output filename", type=str, required=True)
    pve_parser.add_argument("--gwas-method", help="GWAS method used", type=str, choices=['GEMMA'], default='GEMMA')
    return vars(pve_parser.parse_args())


def main():

    # Parse arguments
    pve_args = pveParser()

    # Convert input file to dataframe
    input_dataframe = pd.read_csv(pve_args['input_filename'], sep="\t")

    # Check if the GWAS method is GEMMA
    if pve_args['gwas_method'] == 'GEMMA':

        # Calculate PVE for GEMMA
        input_dataframe['pve'] = (2 * (input_dataframe['beta']**2 * input_dataframe['af'] * (1 - input_dataframe['af']))) / (2 * input_dataframe['beta'] * input_dataframe['af'] * (1 - input_dataframe['af']) + input_dataframe['se']**2 * 2 * 1000 * input_dataframe['af'] * (1 - input_dataframe['af']))

    # Round PVE values to 7 decimal places
    input_dataframe['pve'] = input_dataframe['pve'].round(7)

    # Save the results to output file
    input_dataframe.to_csv(pve_args['output_filename'], sep="\t", index=False)


if __name__ == "__main__":
    main()