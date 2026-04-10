import csv
import argparse

from Bio import SeqIO   
from collections import defaultdict

def updateParser():
    parser = argparse.ArgumentParser(description="Extract the longest transcript per gene from a GFF and FASTA file")
    parser.add_argument("--gff", dest='gff_filename', help="Input GFF file with transcript annotations", type=str, required=True)
    parser.add_argument("--fasta", dest='fasta_filename', help="Input FASTA file with transcript sequences", type=str, required=True)
    parser.add_argument("--output", dest='output_filename', help="Output FASTA file with longest transcripts per gene", type=str, required=True)
    
    return vars(parser.parse_args())


def longest_transcript(gff_filename, fasta_filename, output_filename):

    # Parse the GFF file to map transcripts to genes
    transcript_to_gene = defaultdict(str)

    # Read the GFF file and build the mapping
    with open(gff_filename, 'r') as gff_file:
        gff_reader = csv.reader(gff_file, delimiter='\t')
        for gff_row in gff_reader:
            if gff_row[0].startswith('#'):
                continue

            if 'ID=transcript' in gff_row[8]:
                gene_id = None
                transcript_id = None
                for attr in gff_row[8].split(';'):
                    attr_id, attr_value = attr.split('=')
                    if attr_id == 'ID':
                        transcript_id = attr_value
                    elif attr_id == 'Parent':
                        gene_id = attr_value
                
                if not transcript_id or not gene_id:
                    raise ValueError(f"Could not parse transcript or gene ID from GFF line: {gff_row}")

                transcript_to_gene[transcript_id] = gene_id

    gene_longest_transcript = defaultdict(str)
    gene_longest_length = defaultdict(int)

    # Read the FASTA file and determine the longest transcript for each gene
    for record in SeqIO.parse(fasta_filename, "fasta"):

        # Ensure the transcript ID exists in the mapping
        if record.id not in transcript_to_gene:
            raise ValueError(f"Transcript ID {record.id} not found in GFF mapping.")

        # Get the gene ID and sequence length
        gene_id = transcript_to_gene[record.id]
        seq_length = len(record.seq)

        # Update the longest transcript for the gene if necessary
        if gene_id not in gene_longest_length:
            gene_longest_length[gene_id] = seq_length
            gene_longest_transcript[gene_id] = record
        elif seq_length > gene_longest_length[gene_id]:
            gene_longest_length[gene_id] = seq_length
            gene_longest_transcript[gene_id] = record

    # Write the longest transcripts to the output FASTA file
    with open(output_filename, 'w') as output_file:
        for record in SeqIO.parse(fasta_filename, "fasta"):
            gene_id = transcript_to_gene[record.id]
            if gene_longest_transcript[gene_id].id == record.id:
                SeqIO.write(record, output_file, "fasta")


def main():
    
    # Parse the arguments
    args = updateParser()

    # Run the longest transcript extraction
    longest_transcript(**args)

if __name__ == "__main__":
    main()