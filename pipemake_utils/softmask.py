#!/usr/bin/env python

import argparse

from Bio import SeqIO
from Bio.Seq import Seq

def maskParser ():
	
	# Create the parser
	mask_parser = argparse.ArgumentParser(description = 'Soft mask a genome using a hard-masked genome')

	# Add the arguments
	mask_parser.add_argument('--input-fasta', help = 'Genome fasta', type = str, required = True)
	mask_parser.add_argument('--hard-masked-fasta', help = 'Hard masked genome fasta', type = str, required = True)
	mask_parser.add_argument('--output-fasta', help = 'The output prefix for the soft masked fasta', type = str, default = 'out.fa.masked')
	return vars(mask_parser.parse_args())

def main ():

	# Assign the arguments
	mask_args = maskParser()

	# Creae an index of the hard-masked genome
	hard_masked_dict = SeqIO.index(mask_args['hard_masked_fasta'], 'fasta')

	# Create the soft-masked fasta
	with open(mask_args['output_fasta'], 'w') as soft_masked_fasta:
		for unmasked_record in SeqIO.parse(mask_args['input_fasta'], 'fasta'):
			hard_masked_record = hard_masked_dict[unmasked_record.id]

			if len(unmasked_record) != len(hard_masked_record):
				raise Exception(f'Unable to match unmasked and hard-masked sequence: {unmasked_record.id}')

			# Create a str of the nucleotides to quickly soft mask
			soft_masked_seq = ''

			for i, (unmasked, hard_masked) in enumerate(zip(str(unmasked_record.seq), str(hard_masked_record.seq))):
				if unmasked.upper() == hard_masked.upper(): soft_masked_seq += unmasked
				else: soft_masked_seq += unmasked.lower()

			# Update the sequence
			unmasked_record.seq = Seq(''.join(soft_masked_seq))

			# Write to the output fasta
			SeqIO.write(unmasked_record, soft_masked_fasta, 'fasta')

if __name__ == '__main__':
	main()