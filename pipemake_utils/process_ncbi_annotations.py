import re
import csv
import argparse

from pipemake_utils.misc import *
from pipemake_utils.logger import *

def processID (id_attribute_dict, species_tag):

    for id_str in id_attribute_dict:
        
        id_attribute_dict[id_str] = id_attribute_dict[id_str].replace('%2C', ',')
        id_attribute_dict[id_str] = id_attribute_dict[id_str].replace('egapxtmp_', f'{species_tag}_')

    return id_attribute_dict

def processTranscript (transcript_attribute_dict, species_tag):

    for transcript_str in list(transcript_attribute_dict):

        transcript_attribute_dict[transcript_str] = transcript_attribute_dict[transcript_str].replace('gnl|WGS:ZZZZ|egapxtmp_', f'{species_tag}_')

        if bool(re.search(r'orig_(protein|transcript)_id', transcript_str)):
            del transcript_attribute_dict[transcript_str]

    return transcript_attribute_dict

def processExon (exon_attribute_dict, species_tag):

    for exon_str in list(exon_attribute_dict):

        exon_attribute_dict[exon_str] = exon_attribute_dict[exon_str].replace('gnl|WGS:ZZZZ|egapxtmp_', f'{species_tag}_')

        if bool(re.search(r'orig_(protein|transcript)_id', exon_str)):
            del exon_attribute_dict[exon_str]

    return exon_attribute_dict

def processCDS (cds_attribute_dict, species_tag):

    for cds_str in list(cds_attribute_dict):

        cds_attribute_dict[cds_str] = cds_attribute_dict[cds_str].replace('gnl|WGS:ZZZZ|egapxtmp_', f'{species_tag}_')
        cds_attribute_dict[cds_str] = cds_attribute_dict[cds_str].replace('WGS:ZZZZ:egapxtmp_', f'{species_tag}_')

        if bool(re.search(r'orig_(protein|transcript)_id', cds_str)):
            del cds_attribute_dict[cds_str]

    return cds_attribute_dict

def processCodon (codon_attribute_dict, species_tag):
    
    for codon_str in list(codon_attribute_dict):

        codon_attribute_dict[codon_str] = codon_attribute_dict[codon_str].replace('gnl|WGS:ZZZZ|egapxtmp_', f'{species_tag}_')
        codon_attribute_dict[codon_str] = codon_attribute_dict[codon_str].replace('WGS:ZZZZ:egapxtmp_', f'{species_tag}_')

        if bool(re.search(r'orig_(protein|transcript)_id', codon_str)):
            del codon_attribute_dict[codon_str]

    return codon_attribute_dict

def gtfAttributesToDict (attribute_str):

    attribute_dict = {}
    for attribute in attribute_str.split(';'):
        if attribute.strip() == '':
            continue
        key, value = attribute.strip().split(' ', 1)
        value = value.strip('"')
        attribute_dict[key] = value

    return attribute_dict

def gtfAttributesToGFF (gtf_attributes):

    gff_attributes = {}
    
    for key, value in gtf_attributes.items():
        if key == 'db_xref':
            gff_attributes['Dbxref'] = value
        else:
            gff_attributes[key] = value

    return gff_attributes

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

def gffGeneAttributeStr (gff_gene_attributes):

    for gene_attribute in list(gff_gene_attributes):

        if gene_attribute == 'gene':
             gff_gene_attributes['Name'] = gff_gene_attributes[gene_attribute]
        elif gene_attribute == 'gene_id':
            gff_gene_attributes['ID'] = gff_gene_attributes[gene_attribute]
            del gff_gene_attributes[gene_attribute]
        elif gene_attribute == 'transcript_id':
            del gff_gene_attributes[gene_attribute]

    if 'Name' not in gff_gene_attributes:
        gff_gene_attributes['Name'] = gff_gene_attributes['ID']
    
    return gffAttributesToStr(gff_gene_attributes, order = ['ID', 'Name', 'description'])

def gffTranscriptAttributeStr (gff_transcript_attributes):

    for transcript_attribute in list(gff_transcript_attributes):

        if transcript_attribute == 'gene_id':
            gff_transcript_attributes['Parent'] = gff_transcript_attributes[transcript_attribute]
            del gff_transcript_attributes[transcript_attribute]
        elif transcript_attribute == 'transcript_id':
            gff_transcript_attributes['ID'] = gff_transcript_attributes[transcript_attribute]
            gff_transcript_attributes['Name'] = gff_transcript_attributes[transcript_attribute]
            del gff_transcript_attributes[transcript_attribute]

    return gffAttributesToStr(gff_transcript_attributes, order = ['ID', 'Parent', 'Name', 'product'])

def gffExonAttributeStr (gff_exon_attributes):

    if 'transcript_id' not in gff_exon_attributes and 'number' not in gff_exon_attributes:
        raise ValueError("Exon attributes must contain 'transcript_id' and 'number' keys.")
    
    gff_exon_attributes['ID'] = gff_exon_attributes['transcript_id'] + ':exon-' + gff_exon_attributes['exon_number']

    for exon_attribute in list(gff_exon_attributes):

        if exon_attribute == 'transcript_id':
            gff_exon_attributes['Parent'] = gff_exon_attributes[exon_attribute]
            del gff_exon_attributes[exon_attribute]
        elif exon_attribute == 'exon_number':
            gff_exon_attributes['number'] = gff_exon_attributes[exon_attribute]
            del gff_exon_attributes[exon_attribute]
        elif exon_attribute == 'gene_id':
            del gff_exon_attributes[exon_attribute]

    return gffAttributesToStr(gff_exon_attributes, order = ['ID', 'Parent', 'number'])

def gffCDSAttributeStr (gff_cds_attributes):

    if 'transcript_id' not in gff_cds_attributes and 'number' not in gff_cds_attributes:
        raise ValueError("Exon attributes must contain 'transcript_id' and 'number' keys.")
    
    gff_cds_attributes['ID'] = gff_cds_attributes['transcript_id'] + ':cds-' + gff_cds_attributes['exon_number']

    for exon_attribute in list(gff_cds_attributes):

        if exon_attribute == 'transcript_id':
            gff_cds_attributes['Parent'] = gff_cds_attributes[exon_attribute]
            del gff_cds_attributes[exon_attribute]
        elif exon_attribute == 'exon_number':
            gff_cds_attributes['number'] = gff_cds_attributes[exon_attribute]
            del gff_cds_attributes[exon_attribute]
        elif exon_attribute == 'protein_id':
            del gff_cds_attributes[exon_attribute]

    return gffAttributesToStr(gff_cds_attributes, order = ['ID', 'Parent', 'number'])

def processGTF (gtf_filename, out_prefix, species_tag):

    gff_output = open(f'{out_prefix}.gff', 'w')
    gff_output.write('##gff-version 3\n')

    gtf_output = open(f'{out_prefix}.gtf', 'w')
    gtf_output.write('#gtf-version 2.2\n')

    # Open the GFF file for reading
    with open(gtf_filename, 'r') as gtf_file:
        gtf_reader = csv.reader(gtf_file, delimiter='\t')

        for row in gtf_reader:

            # Check if the row is a comment or header
            if row[0].startswith('#'):
                continue

            # Parse the attributes column into a dictionary
            gtf_attributes = gtfAttributesToDict(row[8])

            # Process attributes based on feature type
            if row[2] in ['gene']:
                pass
            elif row[2] in ['transcript']:
                gtf_attributes = processTranscript(gtf_attributes, species_tag)
            elif row[2] in ['exon']:
                gtf_attributes = processExon(gtf_attributes, species_tag)
            elif row[2] in ['CDS']:
                gtf_attributes = processCDS(gtf_attributes, species_tag)
            elif row[2] in ['start_codon', 'stop_codon']:
                gtf_attributes = processCodon(gtf_attributes, species_tag)
            else:
                raise ValueError(f"Unknown feature type: {row[2]}")
            
            # Process basic IDs for all feature types
            gtf_attributes = processID(gtf_attributes, species_tag)
            
            # Create the gtf attribute string
            gtf_attribute_str = '; '.join([f'{key} "{value}"' for key, value in gtf_attributes.items()]) + ';'

            # Create a the new GTF line with the modified attributes
            gtf_output.write('\t'.join(row[:8] + [gtf_attribute_str]) + '\n')

            # Convert GTF attributes to GFF attributes
            gff_attributes = gtfAttributesToGFF(gtf_attributes)

            # Create a the new GFF line with the modified attributes
            if row[2] in ['gene']:
                gff_attribute_str = gffGeneAttributeStr(gff_attributes)
            elif row[2] in ['transcript']:
                gff_attribute_str = gffTranscriptAttributeStr(gff_attributes)
                if 'mRNA' == gff_attributes['gbkey']:
                    row[2] = 'mRNA'
            elif row[2] in ['exon']:
                gff_attribute_str = gffExonAttributeStr(gff_attributes)
            elif row[2] in ['CDS']:
                gff_attribute_str = gffCDSAttributeStr(gff_attributes)
            elif row[2] in ['start_codon', 'stop_codon']:
                continue
            else:
                raise ValueError(f"Unknown feature type: {row[2]}")
            
            # Write the GFF line to the output file
            gff_output.write('\t'.join(row[:8] + [gff_attribute_str]) + '\n')

    gff_output.close()
    gtf_output.close()

def processNCBIAnnotationsParser():
    parser = argparse.ArgumentParser(description='Process NCBI GTF annotations andd generate GFF and GTF files with modified IDs.')
    parser.add_argument('--gtf-file', dest = 'gtf_filename', help = 'GTF output from the NCBI annotation pipeline', type = str, required = True, action = confirmFile())
    parser.add_argument('--species-tag', help = 'Species tag to replace in IDs', type = str, required = True)
    parser.add_argument('--out-prefix', help = 'Output prefix for GFF and GTF files', type = str, default = 'out')
    
    return vars(parser.parse_args())

def main():

    process_ncbi_args = processNCBIAnnotationsParser()

    # Start logger and log the arguments
    startLogger(f"{process_ncbi_args['out_prefix']}.cvt.log")
    logArgDict(process_ncbi_args)

    processGTF(**process_ncbi_args)

if __name__ == '__main__':
    main()