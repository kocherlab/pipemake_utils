import argparse

from Bio import SeqIO
from Bio.Data import CodonTable

from collections import defaultdict

def countArgs():
    """
    Parse command line arguments for counting codons in a sequence file.
    
    Returns:
        argparse.Namespace: Parsed command line arguments.
    """
    
    parser = argparse.ArgumentParser(description="Count codons in a nucleotide sequence file.")
    parser.add_argument('--input', help = "Path to the input sequence file in FASTA format.", type = str, required = True)
    parser.add_argument('--out-prefix', help = "Prefix for the output file where codon counts will be saved", type = str, default = "out")
    
    return parser.parse_args()


def countStats(sequence_file):

    standard_table = CodonTable.unambiguous_dna_by_name["Standard"]

    # Initialize a dictionary to hold codon counts
    codon_counts = defaultdict(lambda: defaultdict(int))

    # Create ints to hold the total characters and gc characters
    total_chars = 0
    gc_chars = 0
    
    for record in SeqIO.parse(sequence_file, "fasta"):
        seq = record.seq
        # Ensure the sequence length is a multiple of 3
        if len(seq) % 3 != 0:
            raise ValueError("Sequence length is not a multiple of 3.")
        
        # Iterate through the sequence in steps of 3 to get codons
        for i in range(0, len(seq), 3):
            codon = str(seq[i:i+3])
            amino_acid = standard_table.forward_table.get(codon, "Stop")
            codon_counts[amino_acid][codon] += 1

        # Count total characters and GC characters
        total_chars += len(seq)
        gc_chars += seq.count('G') + seq.count('C')

    return codon_counts, total_chars, gc_chars

import sys

codon_counts, total_chars, gc_chars = countStats(sys.argv[1])

total_codons = sum(sum(codons.values()) for codons in codon_counts.values())

# import the standard codon table


print('#BIMPv4.2 Codon Bias')
print(f"#GC Content: {gc_chars / total_chars:.2%}")
print("Codon\tAmino Acid\tFrac\tFreq/1000\tNumber")
for amino_acid, codons in codon_counts.items():
    aa_total = sum(codons.values())
    for codon, count in codons.items():
        codon_frac = count / aa_total
        codon_freq_1000 = (count / total_codons) * 1000
        print(f"{codon}\t{amino_acid}\t{codon_frac:.4f}\t{codon_freq_1000:.2f}\t{count}")