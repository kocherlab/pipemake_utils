import os
import gffutils


from collections import defaultdict
from cyvcf2 import VCF
from Bio import SeqIO

"""
class VCF:
    def __init__(self, vcf_file):
        self.vcf = VCF(vcf_file)
        self.samples = self.vcf.samples
"""


# Check if a GFFutils database already exists, if so don't create it
if not os.path.exists("LBAL_OGS_v3.1.gff3.db"):
    gff_db = gffutils.FeatureDB("LBAL_OGS_v3.1.gff3.db")
else:
    # Open the GFF file
    gff_db = gffutils.create_db(
        "LBAL_OGS_v3.1.gff3", "LBAL_OGS_v3.1.gff3.db"
    )  #':memory:')


# Open the VCF file
vcf_data = VCF("LALB_w_LBAL_v3_f2_merge_ref0.vcf.gz")

# Open the FASTA file
for record in SeqIO.parse(open("LBAL_OGS_v3.1_trans.fa"), "fasta"):
    if len(list(gff_db.children(record.id, featuretype="CDS"))) == 1:
        continue
    if gff_db[record.id].strand != "+":
        continue
    cmp = record.seq
    break

# Index the reference genome
seq_index = SeqIO.index("LBAL_genome_v3.fasta", format="fasta")

# Create a list of CDS coordinates
cds_coords_list = []

# Store the CDS coordinates
for cds in gff_db.children(record.id, featuretype="CDS"):
    cds_coords_list.append([cds.chrom, cds.start, cds.end])

# Sort the CDS coordinates
exon_coords_list = sorted(cds_coords_list, key=lambda x: x[1])

if gff_db[record.id].strand == "-":
    exon_coords_list = exon_coords_list[::-1]

# Create a dictionary to store the CDS sequences
cds_seqs = defaultdict(str)  # type: dict[str, str]

for exon_coords in exon_coords_list:
    exon_seqs = defaultdict(list)

    ref_seq = seq_index[exon_coords[0]][exon_coords[1] - 1 : exon_coords[2]]

    if gff_db[record.id].strand == "-":
        ref_seq = ref_seq.reverse_complement()

    ref_seq = list(str(ref_seq.seq).upper())

    exon_seqs["ref"] = ref_seq

    for sample in vcf_data.samples:
        exon_seqs[f"{sample}_1"] = ref_seq.copy()
        exon_seqs[f"{sample}_2"] = ref_seq.copy()

    for variant in vcf_data(f"{exon_coords[0]}:{exon_coords[1]}-{exon_coords[2]}"):

        def replaceGenotype(gt):
            if gt == 0:
                return variant.REF
            elif gt > 0:
                return variant.ALT[gt - 1]
            else:
                return "N"

        variant_array = variant.genotype.array()[:, :-1]
        for sample_pos, genotype in enumerate(variant_array):
            exon_seqs[f"{vcf_data.samples[sample_pos]}_1"][
                (variant.POS - exon_coords[1])
            ] = replaceGenotype(genotype[0])
            exon_seqs[f"{vcf_data.samples[sample_pos]}_2"][
                (variant.POS - exon_coords[1])
            ] = replaceGenotype(genotype[1])

    for sample, seq in exon_seqs.items():
        cds_seqs[sample] += "".join(seq)

for sample, seq in cds_seqs.items():
    print(f">{sample}\n{seq}")
    break

print(record.format("fasta"))

"""
    for cds in sorted(cds_list, key=lambda x: x[1]):

        cds_seq = seq_index[cds[0]][cds[1] - 1:cds[2]]

        
        if gff_db[record.id].strand == "-":
            cds_seq = cds_seq.reverse_complement()

        transcript_list.append(str(cds_seq.seq).upper())

    if gff_db[record.id].strand == "-":
        transcript_list = transcript_list[::-1]
"""


"""
#print(gff_db[record.id].asinterval)
transcript_list = []

cds_list = []
for cds in gff_db.children(record.id, featuretype='CDS'):
    cds_list.append([cds.chrom, cds.start - 1, cds.end])
for cds in sorted(cds_list, key=lambda x: x[1]):

    cds_seq = seq_index[cds[0]][cds[1]:cds[2]]

    
    if gff_db[record.id].strand == "-":
        cds_seq = cds_seq.reverse_complement()

    transcript_list.append(str(cds_seq.seq).upper())

if gff_db[record.id].strand == "-":
    transcript_list = transcript_list[::-1]

transcript_seq = ''.join(transcript_list)

if cmp.upper() != transcript_seq:
    print(cmp)
    print(transcript_seq)
    print('Error')
sys.exit()
"""

# for variant in test_vcf('20:14000-18000'):


# print(variant.genotypes)
# print(variant)

# print((test_vcf.samples))
