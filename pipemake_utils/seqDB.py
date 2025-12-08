import re
import logging

from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

from collections import defaultdict


def sniffSeqFile(input_filename, file_format="fasta", limit=10):
    def sniff(seq):
        seq_types = {
            "dna": re.compile("^[acgtn\-]*$", re.I),
            "protein": re.compile("^[acdefghiklmnpqrstvwy\*\-]*$", re.I),
        }
        seq_matches = [_t for _t, _s in seq_types.items() if _s.search(str(seq))]

        # Check if format errors
        if not seq_matches:
            raise ValueError("Could not determine the sequence type")
        if len(seq_matches) == 2:
            seq_matches = ["dna"]

        return seq_matches[0]

    logging.info(f"Sniffing the sequence file: {input_filename}")

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

    # Assign the sequence type
    seq_type = max(sniff_dict, key=sniff_dict.get)

    return seq_type


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
        # Parse the fasta file
        for record in SeqIO.parse(self.input_filename, self.file_format):
            # Check if record is divisible by 3
            if self.type == "cds" and len(record.seq) % 3 != 0:
                logging.warning(f"{record.id} is not divisible by 3")
                continue

            # Get the record attributes
            record_attributes = self._parseAttributes(
                record.description, self.limit_attributes
            )

            # Update the record id
            record.id = record_attributes[self.primary_id]
            record.name = record_attributes[self.primary_id]

            # Update the record description
            record.description = self._attributeStr(
                record_attributes, skip_attributes=[self.primary_id]
            )

            yield DBRecord(record, self.primary_id)

    @classmethod
    def read(cls, database, input_filename, file_type, file_format="fasta", **kwargs):
        return cls(database, input_filename, file_type, file_format, **kwargs)

    @classmethod
    def readNCBI(cls, input_filename, file_type, file_format="fasta", **kwargs):
        return cls("ncbi", input_filename, file_type, file_format, **kwargs)

    @classmethod
    def readFlybase(cls, input_filename, file_type, file_format="fasta", **kwargs):
        return cls("flybase", input_filename, file_type, file_format, **kwargs)

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
    def _parsePipemakeAttributes(record_description, limit_attributes=[]):
        record_attributes = {}

        # Assign the gene and protein_id attributes
        record_attributes["gene"] = record_description.rsplit("-", 1)[0]
        record_attributes["protein_id"] = record_description

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


class DBRecord(SeqRecord):
    def __init__(self, record, primary_id):
        super().__init__(
            record.seq,
            id=record.id,
            name=record.name,
            description=record.description,
        )

        # Assign the primary id from the database
        self.primary_id = primary_id
        self._gene_id = None
        self._protein_id = None

    @property
    def gene_id(self):
        # Check if the gene id is already assigned
        if self._gene_id:
            return self._gene_id

        # Assign the gene id
        self._gene_id = self._returnGeneID()
        return self._gene_id

    @property
    def protein_id(self):
        # Check if the gene id is already assigned
        if self._protein_id:
            return self._protein_id

        # Assign the gene id
        self._protein_id = self._returnTranscriptID()
        return self._protein_id

    @property
    def transcript_id(self):
        # Check if the gene id is already assigned
        if self._protein_id:
            return self._protein_id

        # Assign the gene id
        self._protein_id = self._returnTranscriptID()
        return self._protein_id

    def _returnGeneID(self):
        # Check if the gene is the primary id
        if self.primary_id == "gene":
            return self.id

        # Split by [ and ] to get the attributes
        for _s in re.split(r"\[|\]", self.description):
            # Skip if the string is not an attribute
            if "=" not in _s.strip():
                continue
            if "gene=" in _s:
                return _s.split("=")[1]

    def _returnTranscriptID(self):
        # Check if the gene is the primary id
        if self.primary_id == "protein_id":
            return self.id

        # Split by [ and ] to get the attributes
        for _s in re.split(r"\[|\]", self.description):
            # Skip if the string is not an attribute
            if "=" not in _s.strip():
                continue
            if "protein_id=" in _s:
                return _s.split("=")[1]

    def addAttribute(self, key, value):
        self.description += f" [{key}={value}]"
