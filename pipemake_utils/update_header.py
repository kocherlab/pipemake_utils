import logging
import argparse
import gffutils

from Bio import SeqIO

from pipemake_utils.misc import *
from pipemake_utils.logger import *


def strAttributeValue(feature, attribute):
    values = feature.attributes.get(attribute)
    return values[0] if values else None


def getTranscript(gff_db, transcript_id):

    try:
        transcript = gff_db[transcript_id]
    except gffutils.exceptions.FeatureNotFoundError:
        transcript = None

    if transcript is None or transcript.featuretype not in TRANSCRIPT_TYPES:
        raise ValueError(f"Sequence ID {transcript_id} not found among the GFF transcript features ({', '.join(sorted(TRANSCRIPT_TYPES))})")

    return transcript


def buildHeader(gff_db, transcript, id_attr, attributes):

    # The kocherlab ID must be on the transcript line itself (never inherited from the gene)
    new_id = strAttributeValue(transcript, id_attr)
    if not new_id:
        raise ValueError(f"Transcript {transcript.id} has no '{id_attr}' attribute in the GFF")

    # Get the gene of the transcript, if present
    gene = next(gff_db.parents(transcript, featuretype=tuple(GENE_TYPES)), None)

    # Keep the current ID in the header, unless it is already the header's ID
    fields = [("ID", transcript.id)] if id_attr != "ID" else []
    for attribute in attributes:
        if attribute == id_attr:
            continue

        # Assign the attribute value from the transcript if possible
        value = strAttributeValue(transcript, attribute)

        # If the transcript does not have the attribute, assign it from the gene if possible
        if not value and gene is not None:
            value = strAttributeValue(gene, attribute)

        # If the attribute is present, add it to the list of fields
        if value:
            fields.append((attribute, value))

    # Return the new ID and the description (the bracketed fields)
    return new_id, " ".join(f"[{key}={value}]" for key, value in fields)


def updateHeader(gff_file, fasta_file, out_file, id_attr, attributes, **kwargs):

    # Create an in-memory database of the GFF
    gff_db = gffutils.create_db(gff_file, ":memory:", merge_strategy="create_unique", disable_infer_genes=True, disable_infer_transcripts=True)

    header_count = 0
    with open(fasta_file, "r") as fasta_in_file, open(out_file, "w") as fasta_out_file:

        # Iterate the fasta file using biopython SeqIO
        for record in SeqIO.parse(fasta_in_file, "fasta"):

            # gffread uses the transcript ID as the record ID
            transcript = getTranscript(gff_db, record.id)

            # Replace the record ID and description
            record.id, record.description = buildHeader(gff_db, transcript, id_attr, attributes)
            record.name = record.id

            # Write the modified record to the output FASTA file
            SeqIO.write(record, fasta_out_file, "fasta")

            # Increment the header count
            header_count += 1

    logging.info(f"Rewrote {header_count} headers.")


def updateHeaderParser():
    parser = argparse.ArgumentParser(
        description="Rewrite the headers of a gffread FASTA using the GFF: >new_id [ID=...] [attribute=value] ..."
    )
    parser.add_argument("--gff-file", help="GFF3 used to create the FASTA", type=str, required=True, action=confirmFile())
    parser.add_argument("--fasta-file", help="FASTA created by gffread (headers are the transcript IDs)", type=str, required=True, action=confirmFile())
    parser.add_argument("--out-file", help="Output FASTA", type=str, required=True)
    parser.add_argument("--id-attr", help="Transcript attribute holding the ID the header starts with (default: the current ID)", type=str, default="ID")
    parser.add_argument(
        "--attributes",
        help="Attributes to add to the header, in order. Each is taken from the transcript, else its gene",
        type=str,
        nargs="+",
        default=["Name", "description"],
    )

    return vars(parser.parse_args())


def main():

    update_header_args = updateHeaderParser()

    # Start logger and log the arguments
    startLogger()
    logArgDict(update_header_args)

    # Rewrite the headers
    updateHeader(**update_header_args)


if __name__ == "__main__":
    main()
