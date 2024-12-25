import os
import sys
import logging

import numpy as np
import scipy.stats as stats

from Bio import AlignIO
from Bio.Data import CodonTable
from itertools import product
from collections import defaultdict


def missingData(codons):
    for codon in codons:
        if "-" in codon:
            return True
        if "N" in codon:
            return True
    return False


def outToCodons(iupac_codon):
    # Define the IUPAC-to-nucleotide mapping
    iupac_dict = {
        "A": ["A"],
        "C": ["C"],
        "G": ["G"],
        "T": ["T"],
        "R": ["A", "G"],
        "Y": ["C", "T"],
        "S": ["G", "C"],
        "W": ["A", "T"],
        "K": ["G", "T"],
        "M": ["A", "C"],
    }

    # Create a list of nucleotides for each IUPAC code
    nucleotide_list = [iupac_dict[nucl] for nucl in iupac_codon]

    # Return the product of the nucleotide lists
    return ["".join(codon_list) for codon_list in product(*nucleotide_list)]


def isFourFoldDegenerate(codons):
    # Define the four-fold degenerate codons
    ff_degenerates = {"AC", "CC", "CG", "CT", "GC", "GG", "GT", "TC"}

    codons_set = set(str(codon.seq) for codon in codons[:, 0:2])

    if len(codons_set) != 1:
        return False

    if list(codons_set)[0] in ff_degenerates:
        return True
    return False


def mutateCodon(codon, site, pos):
    codon_list = list(codon)
    codon_list[pos] = site
    return "".join(codon_list)


def isPolymorphic(site_set):
    if len(site_set) > 1:
        return True
    return False


def isSynonymous(ref, sample_alleles, pos, out_list):
    codon_table = CodonTable.unambiguous_dna_by_id[1]

    # Assign the outgroup amino acids
    out_aa_list = [codon_table.forward_table[codon] for codon in out_list]

    # Keep track if the codon is synonymous
    synonymous_list = []

    # Loop through the sample alleles
    for sample_allele in sample_alleles:
        sample_aa = codon_table.forward_table[mutateCodon(ref, sample_allele, pos)]
        for out_aa in out_aa_list:
            synonymous_list.append(out_aa == sample_aa)

    return all(synonymous_list)


def isStopCodon(ref, sample_alleles, pos):
    codon_table = CodonTable.unambiguous_dna_by_id[1]

    # Loop through the sample alleles
    for sample_allele in sample_alleles:
        if mutateCodon(ref, sample_allele, pos) in codon_table.stop_codons:
            return True

    return False


def mk(input_file, format="fasta"):
    # Assign the base filename
    base_filename = os.path.basename(input_file)

    # Define the codon table
    codon_table = CodonTable.unambiguous_dna_by_id[1]

    # Read in the alignment file
    alignment = AlignIO.read(input_file, format)

    # Set the alignment length
    alignment_length = alignment.get_alignment_length()

    # Check if the alignment length is divisible by 3
    if alignment_length % 3 != 0:
        logging.warning("Alignment length is not divisible by 3.")
        return {"ps": np.nan, "pn": np.nan, "ds": np.nan, "dn": np.nan}

    # Confirm the first sequence is the reference sequence
    if alignment[0].id not in ["ref", "reference"]:
        raise ValueError(
            "The first sequence in the alignment is not the reference sequence."
        )

    # Confirm the last sequence is the outgroup sequence
    if alignment[-1].id not in ["out", "outgroup"]:
        raise ValueError(
            "The last sequence in the alignment is not the outgroup sequence."
        )

    # Store the different counting categories
    mk_counter = defaultdict(int, {"ps": 0, "pn": 0, "ds": 0, "dn": 0})

    # Loop through the alignment by 3s
    for i in range(0, alignment_length, 3):
        # Check for missing data

        if missingData(alignment[:, i : i + 3]):
            continue

        # Assign the reference and outgroup sequences
        ref = str(alignment[0, i : i + 3].seq)
        out_list = outToCodons(str(alignment[-1, i : i + 3].seq))

        # Remove the reference codon from the outgroup codons, if an alternative is present
        if len(out_list) > 1:
            out_list = [codon for codon in out_list if codon != ref]

        # Check if the reference or outgroup sequences are stop codons
        if ref in codon_table.stop_codons or set(out_list) & set(
            list(codon_table.stop_codons)
        ):
            # Skip the stop codon, if terminal
            if i + 3 >= alignment_length:
                break

            # Log the stop codon
            if ref in codon_table.stop_codons:
                logging.warning(
                    f"Reference stop codon found at position {i} in sequence: {base_filename}"
                )
            else:
                logging.warning(
                    f"Outgroup stop codon found at position {i} in sequence: {base_filename}"
                )
            return {"ps": np.nan, "pn": np.nan, "ds": np.nan, "dn": np.nan}

        # Loop the codon positions
        for j in range(3):
            # Check if the site is monomorphic
            if len(set(nucl for nucl in alignment[:, i + j])) == 1:
                continue

            # Count the number of variants within the codon
            sample_alleles = set(nucl for nucl in alignment[1:-1, i + j])

            # Check if more than one variant is found
            if len(sample_alleles) > 2:
                raise ValueError(
                    f"More than one variant found at position: {i + j + 1}"
                )

            # Check if the site is a stop codon
            if isStopCodon(ref, sample_alleles, j):
                logging.warning(
                    f"Sample stop codon found at position {i} in sequence: {base_filename}"
                )
                return {"ps": np.nan, "pn": np.nan, "ds": np.nan, "dn": np.nan}

            # Check if the site is polymorphic
            is_polymorphic = isPolymorphic(sample_alleles)

            # Check if the codons are synonymous
            is_synonymous = isSynonymous(ref, sample_alleles, j, out_list)

            # Increment the counting categories
            if is_polymorphic and is_synonymous:
                mk_counter["ps"] += 1
            elif is_polymorphic and not is_synonymous:
                mk_counter["pn"] += 1
            elif not is_polymorphic and is_synonymous:
                mk_counter["ds"] += 1
            elif not is_polymorphic and not is_synonymous:
                mk_counter["dn"] += 1

    return mk_counter


def fisherExactTest(pn=0, ps=0, dn=0, ds=0):
    if any(np.isnan([ds, dn, ps, pn])):
        return np.nan, np.nan
    return stats.fisher_exact([[ps, pn], [ds, dn]])


def calcAlpha(pn=0, ps=0, dn=0, ds=0):
    if 0 in [ds, dn, ps, pn]:
        return np.nan
    elif any(np.isnan([ds, dn, ps, pn])):
        return np.nan
    return 1 - (
        (mk_values["ds"] * mk_values["pn"]) / (mk_values["dn"] * mk_values["ps"])
    )


def calcNeutralityIndex(pn=0, ps=0, dn=0, ds=0):
    if 0 in [ds, dn, ps, pn]:
        return np.nan
    elif any(np.isnan([ds, dn, ps, pn])):
        return np.nan
    return (mk_values["pn"] / mk_values["ps"]) / (mk_values["dn"] / mk_values["ds"])


# Assign the sequence directory
seq_dir = sys.argv[1]

with open("mk_results.csv", "w") as mk_output:
    mk_output.write("ID,PS,PN,DS,DN,FET,FET_PVALUE,ALPHA,NEUTRALITY_INDEX\n")

    # Loop the sequence files in the directory
    for seq_file in os.listdir(seq_dir):
        # Define the sequence file path
        seq_path = os.path.join(seq_dir, seq_file)

        # Assign the sequence id
        seq_id = seq_file.rsplit(".", 1)[0]

        # Calculate the MK test values
        mk_values = mk(seq_path)

        # Calculate the Fisher's exact test
        fet_oddsratio, fet_pvalue = fisherExactTest(**mk_values)

        # Calculate the alpha value
        alpha = calcAlpha(**mk_values)

        # Calculate the neutrality index
        neutrality_index = calcNeutralityIndex(**mk_values)

        # Write the results to the output file
        mk_output.write(
            f'{seq_id},{mk_values["ps"]},{mk_values["pn"]},{mk_values["ds"]},{mk_values["dn"]},{fet_oddsratio},{fet_pvalue},{alpha},{neutrality_index}\n'
        )
