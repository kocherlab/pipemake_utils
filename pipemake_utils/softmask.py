#!/opt/conda/envs/bio_python/bin/python

import sys

from Bio import SeqIO
from Bio.Seq import Seq

# Assign the genome files
unmasked_genome = sys.argv[1]
hard_masked_genome = sys.argv[2]
soft_masked_genome = sys.argv[3]

# Creae an index of the hard-masked genome
hard_masked_dict = SeqIO.index(hard_masked_genome, 'fasta')

with open(soft_masked_genome, 'w') as soft_masked_fasta:
	for unmasked_record in SeqIO.parse(unmasked_genome, 'fasta'):
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