#!/usr/bin/env python

import os
import sys
import argparse
import gffutils

from Bio import SeqIO
from itertools import chain
from collections import defaultdict

def processParser ():

	def confirmFile ():
		'''Custom action to confirm file exists'''
		class customAction(argparse.Action):
			def __call__(self, parser, args, value, option_string=None):
				if not os.path.isfile(value):
					raise IOError(f'Unable to find file: {value}')
				setattr(args, self.dest, value)
		return customAction

	# Create the parser
	process_parser = argparse.ArgumentParser(description="process-braker: Confirm and process BRAKER3 output")

	process_parser.add_argument("--fasta-aa", help = "BRAKER3 amino acid fasta filename", type = str, action = confirmFile(), required = True)
	process_parser.add_argument("--fasta-cds", help = "BRAKER3 CDS fasta filename", type = str, action = confirmFile(), required = True)
	process_parser.add_argument("--gff", help = "BRAKER3 gff3 filename", type = str, action = confirmFile(), required = True)
	process_parser.add_argument("--species", help = "Species", type = str, required = True)
	process_parser.add_argument("--assembly-version", help = "Assembly version", type = str, required = True)
	process_parser.add_argument("--annotation-version", help = "Annotation version", type = str, required = True)
	process_parser.add_argument("--out-dir", help = "Output directory", type = str, required = True)

	return vars(process_parser.parse_args())

def createSeqList (filename, seq_format = 'fasta'):

	# Create dict to IDs
	id_dict = defaultdict(list)

	# Read the fasta and store gene and transcript IDs
	for transcript_record in SeqIO.parse(filename, seq_format):
		gene_id, _ = transcript_record.id.split('.')
		id_dict[gene_id].append(transcript_record.id)

	return id_dict

def checkGFF (gff_filename, id_dict):

	# Store DB in memory
	gff_db = gffutils.create_db(gff_filename, ':memory:')

	# Create set of expected features
	expected_features = set(['gene', 'mRNA', 'CDS', 'exon', 'intron', 'start_codon', 'stop_codon'])

	# Create set of the features within the gff
	gff_features = set(gff_db.featuretypes())
	if expected_features != gff_features:
		raise Exception (f'Unable to match GFF features: {gff_features}\n Expected: {expected_features}')
	
	# Loop the ID dict confirm they match the GFF
	for gene_id, transcript_list in id_dict.items():
		
		# Load the gene, if possible
		try: gff_gene = gff_db[gene_id]
		except: raise Exception (f'Unable to locate ID in gff: {gene_id}')

		# Create list of the transcripts of the gene
		gff_mrna = [gff_mrna['ID'] for gff_mrna in gff_db.children(gff_gene, featuretype='mRNA')]
		gff_mrna = set(chain.from_iterable(gff_mrna))

		# Confirm the mRNA IDs match
		if set(transcript_list) != gff_mrna:
			raise Exception (f'Unable to match mRNA IDs between: {set(transcript_list)} and {gff_mrna}')

def createNewIDCvtDict (id_list, annotations_prefix):

	# Create dict to IDs
	new_id_cvt_dict = defaultdict(list)

	# Assign the number of characters to justify
	id_char_len = len(str(len(id_list)))

	# Loop the ID list
	for id_counter, id_str in enumerate(id_list, 1):

		# Create the new ID
		new_id = f"{annotations_prefix}_{str(id_counter).rjust(id_char_len, '0')}"

		# Add new ID to the dict
		new_id_cvt_dict[id_str] = new_id

	return new_id_cvt_dict

def updateSeqFile (filename, id_cvt_dict, out_filename, seq_format = 'fasta'):

	# Open the sequence output filename 
	with open(out_filename, 'w') as out_file:
   
		# Read the fasta and store gene and transcript IDs
		for transcript_record in SeqIO.parse(filename, seq_format):
			
			# Assign the gene ID
			gene_id, _ = transcript_record.id.split('.')

			# Update the record
			transcript_record.description = ''
			transcript_record.id = transcript_record.id.replace(f"{gene_id}.", f"{id_cvt_dict[gene_id]}-")

			# Add the updated record to the output file
			SeqIO.write(transcript_record, out_file, seq_format)

def updateGFF (filename, id_cvt_dict, out_filename):

	def attributesDict (attributes_str):
		
		attributes_dict = {}

		for attribute_str in attributes_str.split(';'):
			if not attribute_str: continue
			attribute_key, attribute_value = attribute_str.split('=')
			attributes_dict[attribute_key] = attribute_value

		return attributes_dict

	# Open the gff output filename 
	out_file = open(out_filename, 'w')

	# Store the old and new gene id
	old_gene_id = None
	new_gene_id = None

	with open(filename, 'r') as gff_input:
		for feature_line in gff_input:

			# Write the header/comments to the out file
			if feature_line.startswith('#'): 
				out_file.write(feature_line)
				continue

			# Assign the basic information of the current feature line
			_, _, feature_type, _, _, _, _, _, feature_attributes_str = feature_line.strip().split('\t')
			
			# Assign the gene ID and update the feature
			if feature_type == 'gene':

				# Assign the gene attributes
				gene_attributes_dict = attributesDict(feature_attributes_str)

				# Confirm the ID was found
				if 'ID' not in gene_attributes_dict:
					raise Exception (f'Unable to assign gene ID from: {feature_attributes_str}')
				
				# Confirm the ID is within the conversion dict
				if gene_attributes_dict['ID'] not in id_cvt_dict:
					raise Exception (f"Unable to convert gene ID: {gene_attributes_dict['ID']}")
							   
				# Update the old and new gene IDs
				old_gene_id = gene_attributes_dict['ID']
				new_gene_id = id_cvt_dict[old_gene_id]
				
				# Update the feature line
				feature_line = feature_line.replace(old_gene_id, new_gene_id)

			else:

				# Update the feature line
				feature_line = feature_line.replace(f'{old_gene_id}.', f'{new_gene_id}-')
				feature_line = feature_line.replace(old_gene_id, new_gene_id)

			# Write the updated feature line
			out_file.write(feature_line)

	# Close the output file
	out_file.close()

def main():

	# Assign the process args
	process_args = processParser()

	# Create an ID list from the amino acid sequences
	aa_id_dict = createSeqList(process_args['fasta_aa'])

	# Create an ID list from the CDS
	cds_id_dict = createSeqList(process_args['fasta_cds'])

	# Confirm the amino acid and CDS IDs are the same
	if aa_id_dict != cds_id_dict: 
		raise Exception (f'Unable to match IDs between: {prot_filename} and {transcript_filename}')

	# Confirm the GFF matches the sequence files
	checkGFF(process_args['gff'], aa_id_dict)

	# Assign the annotations prefix
	annotations_prefix = f"{process_args['species']}_{process_args['assembly_version']}.{process_args['annotation_version']}"

	# Create the new ID conversion dict
	new_id_cvt_dict = createNewIDCvtDict(aa_id_dict, annotations_prefix)

	# Create the output directory, if it does not alraady exits
	if not os.path.exists(process_args['out_dir']): os.makedirs(process_args['out_dir'])

	# Create the output filename prefix
	out_prefix = os.path.join(process_args['out_dir'], f"{process_args['species']}_OGS_{process_args['assembly_version']}.{process_args['annotation_version']}")

	# Update the sequence files
	updateSeqFile(process_args['fasta_aa'], new_id_cvt_dict, f'{out_prefix}_pep.fa')
	updateSeqFile(process_args['fasta_cds'], new_id_cvt_dict, f'{out_prefix}_trans.fa')

	# Update the GFF
	updateGFF(process_args['gff'], new_id_cvt_dict, f'{out_prefix}.gff3')

if __name__ == '__main__':
	main()