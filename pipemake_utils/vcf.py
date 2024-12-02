import os
import gffutils


from collections import defaultdict
from cyvcf2 import VCF
from Bio import SeqIO


class DBVCF:
    def __init__(self, vcf_file):
        self.vcf = vcf_file
        self.samples = self.vcf.samples
        self._ploidy = -1

        # Determine the ploidy of the VCF file
        for vcf_record in self.vcf:
            self._ploidy = vcf_record.ploidy
            break

        # Check if the ploidy was determined
        if self._ploidy == -1:
            raise ValueError("Could not determine the ploidy of the VCF file")

    @classmethod
    def openVCF(cls, vcf_file):
        return cls(VCF(vcf_file))

    def transcipts(self, genome_filename="", gff_filename="", format="fasta"):
        # Index the reference genome
        seq_index = SeqIO.index(genome_filename, format=format)

        # Check if a GFFutils database already exists, if so don't create it
        if os.path.exists(f"{gff_filename}.db"):
            gff_db = gffutils.FeatureDB(f"{gff_filename}.db")
        else:
            # Open the GFF file
            gff_db = gffutils.create_db(f"{gff_filename}", f"{gff_filename}.db")

        # Loop through the transcripts
        for gff_transcript in gff_db.features_of_type(["mRNA", "transcript"]):
            # Create a list of CDS coordinates
            cds_coords_list = []

            # Store the CDS coordinates
            for cds in gff_db.children(gff_transcript, featuretype="CDS"):
                cds_coords_list.append([cds.chrom, cds.start, cds.end])

            # Sort the CDS coordinates
            exon_coords_list = sorted(cds_coords_list, key=lambda x: x[1])

            # Reverse the exon coordinates if the transcript is on the negative strand
            if gff_transcript.strand == "-":
                exon_coords_list = exon_coords_list[::-1]

            # Create a dictionary to store the CDS sequences
            cds_seqs = defaultdict(str)

            # Loop through the exons of the transcript
            for exon_coords in exon_coords_list:
                # Create a dictionary to store the exon sequences
                exon_seqs = defaultdict(list)

                # Get the reference sequence
                ref_seq = seq_index[exon_coords[0]][exon_coords[1] - 1 : exon_coords[2]]

                # Convert the reference sequence to a list
                ref_seq = list(str(ref_seq.seq).upper())

                # Store the reference sequence
                exon_seqs["ref"] = ref_seq

                # Loop through the samples, and store the reference sequence
                for sample in self.samples:
                    for i in range(1, self._ploidy + 1):
                        exon_seqs[f"{sample}_{i}"] = ref_seq.copy()

                # Loop through the variants
                for variant in self.vcf(
                    f"{exon_coords[0]}:{exon_coords[1]}-{exon_coords[2]}"
                ):
                    # Assign the sequence position
                    sequence_pos = variant.POS - exon_coords[1]

                    # Check if the reference sequence matches the VCF file
                    if variant.REF != exon_seqs["ref"][sequence_pos]:
                        raise ValueError("Reference sequence does not match VCF file")

                    # Remove the phase information from the genotype
                    variant_array = variant.genotype.array()[:, :-1]

                    # Update the sequences with their respective genotypes
                    for sample_pos, genotype in enumerate(variant_array):
                        for genotype_pos in range(len(genotype)):
                            exon_seqs[f"{self.samples[sample_pos]}_{genotype_pos + 1}"][
                                sequence_pos
                            ] = self.replaceGenotype(variant, genotype[genotype_pos])

                # Concatenate the sequences
                for sample, seq_list in exon_seqs.items():
                    if gff_transcript.strand == "-":
                        cds_seqs[sample] += self.reverseComplement(seq_list)
                    else:
                        cds_seqs[sample] += "".join(seq_list)

            # for sample, seq in cds_seqs.items():
            #    print(f">{sample}")
            # if sample == 'ref':
            #    print(seq == cmp.upper())
            #    #print()

    @staticmethod
    def reverseComplement(seq_list):
        return "".join(seq_list)[::-1].translate(str.maketrans("ACGT", "TGCA"))

    @staticmethod
    def replaceGenotype(variant, genotype_pos):
        if genotype_pos == 0:
            return variant.REF
        elif genotype_pos > 0:
            return variant.ALT[genotype_pos - 1]
        else:
            return "N"


vcf = DBVCF.openVCF("LBAL_v3_filtered_10.vcf.gz")
vcf.transcipts(
    gff_filename="LBAL_OGS_v3.1.gff3", genome_filename="LBAL_genome_v3.fasta.gz"
)
