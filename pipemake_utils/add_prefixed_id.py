import re
import sys
import logging
import argparse

from pipemake_utils.misc import *
from pipemake_utils.logger import *


def gffAttributesToDict(attribute_str):
    attribute_dict = {}
    for attribute in attribute_str.split(";"):
        if "=" in attribute:
            key, value = attribute.split("=", 1)
            attribute_dict[key] = value
    return attribute_dict


def gtfAttributesToDict(attribute_str):
    attribute_dict = {}

    # Match a key followed by a quoted value (group 2) or an unquoted value (group 3)
    for match in re.finditer(r'\s*([^\s"]+)\s+(?:"([^"]*)"|([^;\s]+))\s*;?', attribute_str):
        attribute_dict[match.group(1)] = match.group(2) if match.group(2) is not None else match.group(3)
    return attribute_dict


def getFeatureID(attribute_dict, id_key):

    # Use the gene_id / transcript_id attribute when present
    if id_key in attribute_dict:
        return attribute_dict[id_key]

    # Otherwise fall back to the GFF3 ID, removing the "gene:" / "transcript:" prefix
    raw_id = attribute_dict.get("ID", "")
    return raw_id.split(":", 1)[1] if ":" in raw_id else raw_id


def appendAttribute(attribute_str, annotation_format, key, value):
    attribute_str = attribute_str.rstrip().rstrip(";")
    if annotation_format == "gtf":
        return f'{attribute_str}; {key} "{value}";'
    return f"{attribute_str};{key}={value}"


def addPrefixedID(
    in_filename,
    out_filename,
    annotation_format,
    prefix,
    sep,
    attr_name,
    gene_types,
    transcript_types,
    **kwargs,
):

    # Select the attribute parser for the annotation format
    attributesToDict = gtfAttributesToDict if annotation_format == "gtf" else gffAttributesToDict

    # Count the number of features updated
    feature_counts = {"gene": 0, "transcript": 0, "missing_id": 0}

    with open(in_filename, "r") as annotation_file, open(out_filename, "w") as out_file:
        for line in annotation_file:

            # Pass comments and blank lines through
            if line.startswith("#") or not line.strip():
                out_file.write(line)
                continue

            # Split the line into columns
            annotation_columns = line.rstrip().split("\t")

            # Check that the line has at least 9 columns (GFF3/GTf format)
            if len(annotation_columns) < 9:
                out_file.write(line)
                continue

            # Only gene-level and transcript-level features are updated
            feature_type = annotation_columns[2]
            if feature_type in gene_types:
                feature_level, id_key = "gene", "gene_id"
            elif feature_type in transcript_types:
                feature_level, id_key = "transcript", "transcript_id"
            else:
                out_file.write(line)
                continue

            # Convert the attribute string to a dictionary for easier manipulation
            attribute_dict = attributesToDict(annotation_columns[8])

            # Get the current ID, skipping features without an ID
            current_id = getFeatureID(attribute_dict, id_key)

            # Raise an error if there is no ID for the feature. This may indicate a problem, incomplete, or a non-standard GFF3/GTF file.
            if not current_id:
                raise Exception(f"Feature {feature_type} on line {line.strip()} has no ID; cannot add prefixed ID")

            # Raise an error if the attribute already exists
            if attr_name in attribute_dict:
                raise Exception(f"Attribute '{attr_name}' already exists on a {feature_type} line; choose another --attr-name")

            # Append the new attribute to the attribute string and write the modified line to the output file
            annotation_columns[8] = appendAttribute(annotation_columns[8], annotation_format, attr_name, f"{prefix}{sep}{current_id}")

            # Write the modified line to the output file
            out_file.write("\t".join(annotation_columns) + "\n")

            # Increment the count of features updated
            feature_counts[feature_level] += 1

    logging.info(
        f"Added {attr_name} to {feature_counts['gene']} gene-level and {feature_counts['transcript']} "
        f"transcript-level features ({feature_counts['missing_id']} skipped, no ID)."
    )


def addPrefixedIDParser():
    parser = argparse.ArgumentParser(
        description="Add a prefixed-ID attribute (<prefix><sep><current ID>) to the gene-level and transcript-level features of a GFF3 or GTF."
    )
    parser.add_argument("--in-file", dest="in_filename", help="Input GFF3 or GTF (uncompressed)", type=str, required=True, action=confirmFile())
    parser.add_argument("--out-file", dest="out_filename", help="Output filename", type=str, required=True)
    parser.add_argument("--format", dest="annotation_format", help="Input/output format", choices=["gff3", "gtf"], type=str, required=True)
    parser.add_argument("--prefix", help="Prefix prepended to the current ID", type=str, required=True)
    parser.add_argument("--sep", help="Separator placed between the prefix and the ID", type=str, default="_")
    parser.add_argument("--attr-name", help="Name of the attribute to add", type=str, default="kocherlab_id")
    parser.add_argument("--gene-types", help="Gene-level feature types", type=str, nargs="+", default=GENE_TYPES)
    parser.add_argument("--transcript-types", help="Transcript-level feature types", type=str, nargs="+", default=TRANSCRIPT_TYPES)

    return vars(parser.parse_args())


def main():

    add_prefixed_id_args = addPrefixedIDParser()

    # Start logger and log the arguments
    startLogger()
    logArgDict(add_prefixed_id_args)

    # Add the prefixed ID to the annotation file
    addPrefixedID(**add_prefixed_id_args)


if __name__ == "__main__":
    main()
