import csv
import argparse

from pipemake_utils.misc import *
from pipemake_utils.logger import *

def gffAttributesToDict (gff_attributes):

    gff_attributes_dict = {}

    for attribute in gff_attributes.split(';'):
        if not attribute:
            continue
        key, value = attribute.split('=')
        if key in gff_attributes_dict:
            if isinstance(gff_attributes_dict[key], list):
                gff_attributes_dict[key].append(value)
            else:
                gff_attributes_dict[key] = [gff_attributes_dict[key], value]
        else:
            gff_attributes_dict[key] = value

    return gff_attributes_dict

def gffAttributesToStr (gff_attributes, order = []):

    gff_attribute_str = ''

    if order:

        for ordered_gff_attribute in order:
            if ordered_gff_attribute not in gff_attributes:
                continue
            gff_attribute_str += f'{ordered_gff_attribute}={gff_attributes[ordered_gff_attribute]};'
            del gff_attributes[ordered_gff_attribute]

    gff_attribute_str += ';'.join([f'{key}={value}' for key, value in gff_attributes.items()])

    if gff_attribute_str.endswith(';'):
        gff_attribute_str = gff_attribute_str[:-1]

    return gff_attribute_str

def gtfAttributesToStr (gtf_attributes, order = []):

    gtf_attribute_str = ''

    if order:

        for ordered_gtf_attribute in order:
            if ordered_gtf_attribute not in gtf_attributes:
                continue
            gtf_attribute_str += f'{ordered_gtf_attribute} "{gtf_attributes[ordered_gtf_attribute]}"; '
            del gtf_attributes[ordered_gtf_attribute]

    gtf_attribute_str += '; '.join([f'{key} "{value}"' for key, value in gtf_attributes.items()])

    if not gtf_attribute_str.endswith(';'):
        gtf_attribute_str += ';'
    
    return gtf_attribute_str

def gffGeneAttributeStr (gff_gene_attributes, prefix = 'MT'):

    gff_gene_attributes['Name'] = f"{prefix}-{gff_gene_attributes['Name']}"
    gff_gene_attributes['ID'] = f"{prefix}-{gff_gene_attributes['ID']}"

    return gffAttributesToStr(gff_gene_attributes, order = ['ID', 'Parent', 'Name', 'product'])

def gtfGeneAttributeStr (gff_gene_attributes, prefix = 'MT'):

    gtf_gene_attributes = {}

    gtf_gene_attributes['gene_id'] = f"{prefix}-{gff_gene_attributes['ID']}"
    gtf_gene_attributes['gbkey'] = 'Gene'
    gtf_gene_attributes['transcript_id'] = ''

    return gtfAttributesToStr(gtf_gene_attributes, order = ['gene_id', 'transcript_id', 'gbkey'])

def gffTranscriptAttributeStr (gff_transcript_attributes, prefix = 'MT'):

    if 'Name' in gff_transcript_attributes:
        gff_transcript_attributes['Name'] = f"{prefix}-{gff_transcript_attributes['Name']}"
    gff_transcript_attributes['ID'] = f"{prefix}-{gff_transcript_attributes['ID']}"
    gff_transcript_attributes['Parent'] = f"{prefix}-{gff_transcript_attributes['Parent']}"

    return gffAttributesToStr(gff_transcript_attributes, order = ['ID', 'Parent', 'Name', 'product'])

def gffExonAttributeStr (gff_exon_attributes):

    gff_exon_attributes['Parent'] = f"MT-{gff_exon_attributes['Parent']}"
    
    return gffAttributesToStr(gff_exon_attributes, order = ['ID', 'Parent'])

def gffCDSAttributeStr (gff_cds_attributes, prefix = 'MT'):

    gff_cds_attributes['ID'] = f"{prefix}-{gff_cds_attributes['ID']}"
    gff_cds_attributes['Parent'] = f"{prefix}-{gff_cds_attributes['Parent']}"
    gff_cds_attributes['Name'] = f"{prefix}-{gff_cds_attributes['Name']}"

    return gffAttributesToStr(gff_cds_attributes, order = ['ID', 'Parent', 'Name'])

def processGFF (gff_filename, out_prefix, species_tag):

    gff_output = open(f'{out_prefix}.gff', 'w')

    # Open the GFF file for reading
    with open(gff_filename, 'r') as gff_file:
        gff_reader = csv.reader(gff_file, delimiter='\t')

        for row in gff_reader:
    
            # Stop processing if we reach the FASTA section
            if row[0].startswith('##FASTA'):
                break

            # Check if the row is a comment or header
            if row[0].startswith('#'):
                continue

            row[0] = f'{species_tag}_chr_MT'
            row[1] = 'mitohifi'

            # Parse the attributes column into a dictionary
            gff_attributes = gffAttributesToDict(row[8])

            # Create a the new GFF line with the modified attributes
            if row[2] in ['gene']:
                gff_attribute_str = gffGeneAttributeStr(gff_attributes)
            elif row[2] in ['mRNA', 'tRNA', 'rRNA']:
                gff_attribute_str = gffTranscriptAttributeStr(gff_attributes)
            elif row[2] in ['exon']:
                gff_attribute_str = gffExonAttributeStr(gff_attributes)
            elif row[2] in ['CDS']:
                gff_attribute_str = gffCDSAttributeStr(gff_attributes)
            elif row[2] in ['DNA']:
                continue
            else:
                raise ValueError(f"Unknown feature type: {row[2]}")
            
            # Write the GFF line to the output file
            gff_output.write('\t'.join(row[:8] + [gff_attribute_str]) + '\n')

    gff_output.close()

def processNCBIAnnotationsParser():
    parser = argparse.ArgumentParser(description='Process mitohifi GFF annotations andd generate GFF and GTF files with modified IDs.')
    parser.add_argument('--gff-file', dest = 'gff_filename', help = 'GFF output from the NCBI annotation pipeline', type = str, required = True, action = confirmFile())
    parser.add_argument('--species-tag', help = 'Species tag to replace in IDs', type = str, required = True)
    parser.add_argument('--out-prefix', help = 'Output prefix for GFF and GTF files', type = str, default = 'out')
    
    return vars(parser.parse_args())

def main():

    process_ncbi_args = processNCBIAnnotationsParser()

    # Start logger and log the arguments
    startLogger(f"{process_ncbi_args['out_prefix']}.cvt.log")
    logArgDict(process_ncbi_args)

    processGFF(**process_ncbi_args)

if __name__ == '__main__':
    main()