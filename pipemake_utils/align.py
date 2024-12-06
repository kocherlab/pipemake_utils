import os
import logging

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def confirmCDS(input_filename, output_filename, file_format="fasta"):
    # Create a list to store the records
    records = []

    # Loop through the records
    for record in SeqIO.parse(input_filename, file_format):
        # Check if the sequence is divisible by 3
        if len(record.seq) % 3 != 0:
            logging.warning(
                f"{input_filename}: {record.id} is not divisible by 3. File will not be written."
            )
            records = []
            break

        # Translate the sequence
        translated_seq = record.seq.translate()

        # Check for internal stop codons
        if "*" in translated_seq and not translated_seq.endswith("*"):
            logging.warning(
                f"{input_filename}: {record.id} has internal stop codons. File will not be written."
            )
            records = []
            break

        # Append the record to the list
        records.append(record)

    # Check if there are records
    if not records:
        return

    # Write the records to the output file
    with open(output_filename, "w") as out_handle:
        SeqIO.write(records, out_handle, file_format)

    logging.info(f"Confirmed the CDS sequences in the file: {input_filename}")


def translateCDS(input_filename, output_filename, file_format="fasta"):
    with open(output_filename, "w") as out_handle:
        for record in SeqIO.parse(input_filename, file_format):
            record.seq = record.seq.translate(to_stop=True)
            SeqIO.write(record, out_handle, "fasta")

    logging.info(f"Translated the alignment file: {output_filename}")


def mapGapsFromIndex(
    index_filename,
    msa_filename,
    output_filename,
    seq_format="fasta",
    confirm_match=True,
):
    # Create the sequence index
    seq_index = SeqIO.index(index_filename, seq_format)

    # Create the output file
    output_file = open(output_filename, "w")

    with open(msa_filename) as msa_file:
        for msa_record in SeqIO.parse(msa_file, seq_format):
            # Get the DNA sequence for the record
            dna_record = seq_index[msa_record.id]

            # Map the DNA sequence
            aa_count = 0
            mapped_seq = ""
            for amino_acid in msa_record.seq:
                # Map the gap
                if amino_acid == "-":
                    mapped_seq += "---"

                # Map the codon
                else:
                    # Assign and map the codon
                    codon = dna_record.seq[aa_count * 3 : (aa_count * 3) + 3]
                    mapped_seq += codon
                    aa_count += 1

                    # Confirm the codon, raise exception if not possible
                    if confirm_match and amino_acid != codon.translate():
                        output_file.close()
                        os.remove(output_filename)
                        raise Exception(
                            f"Codon does not match to amino acid in {msa_record.id}"
                        )

            # Create the mapped record
            mapped_record = SeqRecord(
                Seq(mapped_seq), id=msa_record.id, description=msa_record.description
            )

            # Check if the mapped record is the same as the original
            if mapped_record.seq.translate() != msa_record.seq:
                output_file.close()
                os.remove(output_filename)
                raise Exception(
                    f"Mapped sequence does not match original sequence in {msa_record.id}"
                )

            # Save the mapped record
            SeqIO.write(mapped_record, output_file, seq_format)

    # Close the output file
    output_file.close()
