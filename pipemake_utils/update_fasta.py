
import csv
import argparse

from Bio import SeqIO

from pipemake_utils.misc import *
from pipemake_utils.logger import *

def gffAttributesToDict (attribute_str):
    attribute_dict = {}
    for attribute in attribute_str.split(';'):
        if attribute.strip():
            key_value = attribute.split('=')
            if len(key_value) == 2:
                key, value = key_value
                attribute_dict[key.strip()] = value.strip()
            else:
                pass
    return attribute_dict


def generateGffDict(gff_filename, attributes, **kwargs):

    # Initialize a dictionary to hold mRNA attributes for the GFF file
    gff_mRNA_attribute_dict = {}

    # Add the exception attribute to the arguments
    attributes.append('exception')

    # Open the GFF file for reading
    with open(gff_filename, 'r') as gff_file:
        gff_reader = csv.reader(gff_file, delimiter='\t')
        for row in gff_reader:

            # Check if the row is a comment or header
            if row[0].startswith('#') or row[2] != 'mRNA':
                continue

            # Parse the attributes column into a dictionary
            mRNA_attribute_dict = gffAttributesToDict(row[8])

            if 'ID' not in mRNA_attribute_dict:
                raise ValueError(f"ID attribute not found in GFF attributes: {row[8]}")

            # Extract the ID and other attributes
            mRNA_id = mRNA_attribute_dict['ID']

            # Filter attributes based on the provided list
            mRNA_attribute_dict = {_attr: mRNA_attribute_dict[_attr] for _attr in attributes if _attr in mRNA_attribute_dict}

            # Add the mRNA ID and its attributes to the dictionary
            gff_mRNA_attribute_dict[mRNA_id] = mRNA_attribute_dict

    return gff_mRNA_attribute_dict

def checkDNAForInternalStopCodon(sequence):
    stop_codons = ['TAA', 'TAG', 'TGA']

    for codon in range(0, len(str(sequence)), 3):
        if sequence[codon: codon + 3] in stop_codons:
            return True

def updateFasta(gff_attribute_dict, fasta_filename, out_filename, seq_type, **kwargs):
    
    with open(f"{out_filename}", 'w') as fasta_out_file, open(fasta_filename, 'r') as fasta_in_file:

        # Iterate thhe fasta file using biopython SeqIO
        for record in SeqIO.parse(fasta_in_file, 'fasta'):

            if record.id not in gff_attribute_dict:
                raise ValueError(f"Sequence ID {record.id} not found in GFF attributes dictionary.")

            attributes = gff_attribute_dict[record.id]            
           
            # Update the record ID with the attributes from the GFF dictionary
            record.id +=  ' ' + ' '.join([f"[{key}={value}]" for key, value in attributes.items()])
            record.description = ""  # Clear the description to avoid duplication
                
            # Check for internal stop codons within a protein sequence
            if seq_type == 'protein' and '.' in record.seq:
                
                logging.info(f"Internal stop codon found in sequence {record.id} flagging and replacing with 'X'.")

                # Replace the internal stop codon with 'X'
                record.seq = record.seq.replace('.', 'X')

                # Add a warning to the record ID
                record.id += ' [warn=internal_stop_codon]'

            # Check for internal stop codons within a nucleotide sequence
            elif seq_type == 'nucleotide' and checkDNAForInternalStopCodon(record.seq):

                logging.info(f"Internal stop codon found in nucleotide sequence {record.id} flagging.")

                # Add a warning to the record ID
                record.id += ' [warn=internal_stop_codon]'

            # Write the modified record to the output FASTA file
            SeqIO.write(record, fasta_out_file, 'fasta')

def processNCBIAnnotationsParser():
    parser = argparse.ArgumentParser(description='Process NCBI GTF annotations andd generate GFF and GTF files with modified IDs.')
    parser.add_argument('--gff-file', dest = 'gff_filename', help = 'GTF output from the NCBI annotation pipeline', type = str, required = True, action = confirmFile())
    parser.add_argument('--fasta-file', dest = 'fasta_filename', help = 'FASTA file with sequences', type = str, required = True, action = confirmFile())
    parser.add_argument('--attributes', help = 'One or more GFF attributes to extract from the GTF file', type = str, nargs = '+', default = ['gene', 'product'])
    parser.add_argument('--seq-type', choices=['nucleotide', 'protein'], help = 'Type of sequences in the FASTA file', type = str, required = True)
    parser.add_argument('--out-file', dest = 'out_filename', help = 'Output prefix for FASTA', type = str, default = 'out.fa')
    
    return vars(parser.parse_args())


def main():

    process_ncbi_args = processNCBIAnnotationsParser()

    # Start logger and log the arguments
    startLogger(f"{process_ncbi_args['out_filename']}.log")
    logArgDict(process_ncbi_args)

    # Generate the GFF dictionary with mRNA attributes
    gff_attribute_dict = generateGffDict(**process_ncbi_args)

    # Update the FASTA file with the attributes from the GFF dictionary
    updateFasta(gff_attribute_dict, **process_ncbi_args)

if __name__ == '__main__':
    main()