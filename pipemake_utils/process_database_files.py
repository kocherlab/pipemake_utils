import re
import logging
import argparse

from Bio import SeqIO
from collections import defaultdict


def sniffSeqFile(input_filename, file_format="fasta", limit=10):
    def sniff(seq):
        seq_types = {
            "dna": re.compile("^[acgtn]*$", re.I),
            "protein": re.compile("^[acdefghiklmnpqrstvwy\*]*$", re.I),
        }
        seq_matches = [_t for _t, _s in seq_types.items() if _s.search(str(seq))]

        # Check if format errors
        if not seq_matches:
            raise ValueError("Could not determine the sequence type")
        if len(seq_matches) > 1:
            raise ValueError("Multiple sequence types found")

        return seq_matches[0]

    # Create dict to store sniffs
    sniff_dict = defaultdict(int)

    # Parse the fasta file
    for sniffs, record in enumerate(SeqIO.parse(input_filename, file_format)):
        if sniffs >= limit:
            break
        sniff_dict[sniff(record.seq)] += 1

    # Check if no sequences found
    if not sniff_dict:
        raise ValueError("No sequences found in the file")

    return sorted(sniff_dict.items(), key=lambda item: item[1], reverse=True)[0][0]


class DBFileReader:
    def __init__(
        self,
        database,
        input_filename,
        file_type,
        file_format="fasta",
        primary_id="gene",
        limit_attributes=["gene", "protein", "protein_id"],
    ):
        self.database = database.lower()
        self.input_filename = input_filename
        self.type = file_type.lower()
        self.file_format = file_format
        self.primary_id = primary_id
        self.limit_attributes = limit_attributes

    def __iter__(self):
        if self.type == "dna":
            return self._readDNA()
        elif self.type == "protein":
            return self._readAAs()

    @classmethod
    def read(cls, database, input_filename, file_type, file_format="fasta"):
        return cls(database, input_filename, file_type, file_format)

    @classmethod
    def readNCBI(cls, input_filename, file_type, file_format="fasta"):
        return cls("ncbi", input_filename, file_type, file_format)

    @classmethod
    def readFlybase(cls, input_filename, file_type, file_format="fasta"):
        return cls("flybase", input_filename, file_type, file_format)

    def _readDNA(self):
        # Parse the fasta file
        for record in SeqIO.parse(self.input_filename, self.file_format):
            # Check if record is divisible by 3
            if len(record.seq) % 3 != 0:
                logging.warning(f"{record.id} is not divisible by 3")
                continue

            # Get the record attributes
            record_attributes = self._parseAttributes(
                record.description, self.limit_attributes
            )

            # Update the record id
            record.id = record_attributes[self.primary_id]

            # Update the record description
            record.description = self._attributeStr(
                record_attributes, skip_attributes=[self.primary_id]
            )

            yield record

    def _readAAs(self):
        # Parse the fasta file
        for record in SeqIO.parse(self.input_filename, self.file_format):
            # Get the record attributes
            record_attributes = self._parseAttributes(
                record.description, self.limit_attributes
            )

            # Update the record id
            record.id = record_attributes[self.primary_id]

            # Update the record description
            record.description = self._attributeStr(
                record_attributes, skip_attributes=[self.primary_id]
            )

            yield record

    def _parseAttributes(self, *args, **kwargs):
        # Parse NCBI database attributes
        if self.database == "ncbi":
            return self._parseNCBIAttributes(*args, **kwargs)
        elif self.database == "flybase":
            return self._parseFlybaseAttributes(*args, **kwargs)
        elif self.database == "pipemake":
            return self._parsePipemakeAttributes(*args, **kwargs)

    @staticmethod
    def _parseNCBIAttributes(record_description, limit_attributes=[]):
        record_attributes = {}

        # Split by [ and ] to get the attributes
        for _s in re.split(r"\[|\]", record_description):
            # Skip if the string is not an attribute
            if "=" not in _s.strip():
                continue

            # Split by = to get the key and value of the attribute
            attribute_dict = _s.strip().split("=")

            # Skip if the attribute is not in the limit_attributes
            if limit_attributes and attribute_dict[0] not in limit_attributes:
                continue

            # Update the record_attributes dictionary
            record_attributes[attribute_dict[0]] = attribute_dict[1]

        return record_attributes

    @staticmethod
    def _parsePipemakeAttributes(record_description, limit_attributes=[]):
        record_attributes = {}

        # Assign the gene and protein_id attributes
        record_attributes["gene"] = record_description.split("-")[0]
        record_attributes["protein_id"] = record_description

        return record_attributes

    @staticmethod
    def _parseFlybaseAttributes(record_description, limit_attributes=[]):
        record_attributes = {}

        # Split by [ and ] to get the attributes
        for _s in re.split(r";|\ ", record_description):
            # Skip if the string is not an attribute
            if "=" not in _s.strip():
                continue

            # Split by = to get the key and value of the attribute
            attribute_dict = _s.strip().split("=")

            if "parent" == attribute_dict[0]:
                attribute_dict[0] = "gene"
                attribute_dict[1] = attribute_dict[1].split(",")[0]
                if "FBgn" not in attribute_dict[1]:
                    raise ValueError(
                        f"parent attribute does not contain FBgn: {attribute_dict[1]}"
                    )

            elif "name" == attribute_dict[0]:
                attribute_dict[0] = "protein"

            elif "dbxref" == attribute_dict[0]:
                attribute_dict[0] = "protein_id"
                for _dbxref in attribute_dict[1].split(","):
                    if "FBpp" in _dbxref.split(":")[1]:
                        attribute_dict[1] = _dbxref.split(":")[1]
                        break
                if "FBpp" not in attribute_dict[1]:
                    raise ValueError(
                        f"dbxref attribute does not contain FBpp: {attribute_dict[1]}"
                    )

            # Skip if the attribute is not in the limit_attributes
            if limit_attributes and attribute_dict[0] not in limit_attributes:
                continue

            # Update the record_attributes dictionary
            record_attributes[attribute_dict[0]] = attribute_dict[1]

        return record_attributes

    @staticmethod
    def _attributeStr(record_attributes, skip_attributes=[]):
        attribute_str = ""

        # Loop through the record attributes
        for _k, _v in record_attributes.items():
            # Add space if the attribute_str is not empty
            if attribute_str:
                attribute_str += " "

            # Skip if the attributes if in skip_attributes
            if skip_attributes and _k in skip_attributes:
                continue

            # Update the attribute_str
            attribute_str += "[%s=%s]" % (_k, _v)

        return attribute_str


def databaseParser():
    db_parser = argparse.ArgumentParser(description="Process NCBI fasta files")

    # Input and output arguments
    db_parser.add_argument(
        "--database", help="Source database for the input file", required=True
    )
    db_parser.add_argument("--input-filename", help="Input filename", required=True)
    db_parser.add_argument("--output-prefix", help="Output prefix", required=True)

    # Optional arguments
    db_parser.add_argument(
        "--no-cds", help="Don't output CDS sequences", action="store_true"
    )
    db_parser.add_argument("--input-type", help="Input type (e.g. DNA, AA)")
    db_parser.add_argument("--input-format", help="Input format", default="fasta")
    db_parser.add_argument("--output-format", help="Output format", default="fasta")
    db_parser.add_argument(
        "--output-primary-id", help="Output primary id", default="gene"
    )
    db_parser.add_argument(
        "--log-filename", help="Log filename", default="process_ncbi_cds.log"
    )

    return vars(db_parser.parse_args())


def main():
    # Assign the arguments from the command-line
    db_args = databaseParser()

    # Set the logging configuration
    logging.basicConfig(
        filename=db_args["log_filename"],
        level=logging.INFO,
        format="%(asctime)s %(levelname)s %(message)s",
        datefmt="%Y-%m-%d %H:%M:%S",
    )

    # Check if the input type is provided
    if not db_args["input_type"]:
        # Sniff the sequence file, if not provided
        db_args["input_type"] = sniffSeqFile(db_args["input_filename"])

    # Warn if the input cannot create the CDS
    if db_args["input_type"] == "protein" and not db_args["no_cds"]:
        logging.warning("Protein sequences cannot be converted to CDS")
        db_args["no_cds"] = True

    # Open files to store the processed CDS and protein sequences
    protein_file = open(f"{db_args['output_prefix']}_aa.fa", "w")
    if not db_args["no_cds"]:
        cds_file = open(f"{db_args['output_prefix']}_cds.fa", "w")

    for record in DBFileReader.read(
        db_args["database"],
        db_args["input_filename"],
        file_type=db_args["input_type"],
        file_format=db_args["input_format"],
    ):
        # Write the record to the cds file, if not no_cds
        if not db_args["no_cds"]:
            SeqIO.write(record, cds_file, db_args["output_format"])

        # Translate the DNA record
        if db_args["input_type"] == "dna":
            record.seq = record.seq.translate()

        # Write the record to the protein file
        SeqIO.write(record, protein_file, db_args["output_format"])

    # Close the files
    protein_file.close()
    if not db_args["no_cds"]:
        cds_file.close()


if __name__ == "__main__":
    main()
