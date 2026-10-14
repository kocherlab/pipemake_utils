#!/usr/bin/env python3

import os
import sys
import json
import logging
import argparse

from ete4 import Tree
from ete4.treeview import TreeStyle
from collections import defaultdict



class GeneTree:
    def __init__(self, tree_file, newick_str = '', file_type = 'newick', invalid_action = 'warning'):

        self.tree_file = os.path.basename(tree_file)
        self._file_type = file_type
        self._origin_tree = Tree(newick_str, parser=1) if newick_str else Tree(tree_file, parser=1)
        self._labelled_tree = None

        if invalid_action not in ['ignore', 'warning', 'error']:
            raise ValueError(f"Invalid default_action: {invalid_action}. Must be one of 'ignore', 'warning', or 'error'.")
        self._invalid_action = invalid_action

        self._gene_to_species_name = {}
        self._valid_tree = True
        self._species_duplicates = defaultdict(list)
        self._species_uniques = {}
        self._species_tree_leaf_names = []

    @classmethod
    def fromNewick(cls, tree_file, invalid_action = 'warning'):

        return cls(tree_file, invalid_action=invalid_action)

    @classmethod
    def fromHyphy(cls, json_file, index = '0', invalid_action = 'warning'):

        # Read the JSON file and extract the Newick tree
        with open(json_file) as f:
            trees = json.load(f)["input"]["trees"]
                    
            if index not in trees:
                sys.exit(f"Tree index {index!r} not found; available: {list(trees)}")
            
            newick_str = trees[index]
            
            # Some HyPhy versions store the Newick string directly, others wrap it in a dict.
            if isinstance(newick_str, dict):
                newick_str = newick_str.get("string") or newick_str.get("newick") or next(iter(newick_str.values()))
    
            # Strip whitespace and ensure the tree ends with a semicolon.
            newick_str = newick_str.strip()
    
            # Ensure the tree ends with a semicolon.
            if not newick_str.endswith(";"):
                newick_str += ";"

        return cls(json_file, newick_str = newick_str, file_type = 'hyphy', invalid_action=invalid_action)

    def isValid (self):          
        return self._valid_tree

    def mapSpeciesTree (self, species_tree):

        # Store the leaf names of the species tree as a set
        self._species_tree_leaf_names = set(species_tree.leaf_names())

        # Store the species to gene names mapping
        species_to_gene_names = defaultdict(set)

        # Map the gene tree leaves to the species tree leaves
        for leaf in self._origin_tree.leaf_names():
            leaf_assigned = False
            for species_name in self._species_tree_leaf_names:
                if leaf.startswith(species_name):
                    species_to_gene_names[species_name].add(leaf)
                    self._gene_to_species_name[leaf] = species_name
                    leaf_assigned = True
                    break

            # If no matching species name is found for the leaf, log an error
            if not leaf_assigned:
                logging.error(f"Leaf '{leaf}' in tree '{self.tree_file}' does not match any species in the species tree.")
                self._valid_tree = False

        # Assign the duplicate status
        for species, gene_names in species_to_gene_names.items():
            if len(gene_names) == 1:
                self._species_uniques[species] = list(gene_names)[0]
            elif gene_names == set(self._origin_tree.common_ancestor(gene_names).leaf_names()):
                self._species_duplicates[species] = [gene_names, True]
            else:
                if self._invalid_action == 'error':
                    raise ValueError(f"Non-sister duplicates found for species ({species}) in {self.tree_file}")
                elif self._invalid_action == 'warning':
                    logging.warning(f"Non-sister duplicates found for species ({species}) in {self.tree_file}")
                    self._valid_tree = False
                elif self._invalid_action == 'ignore':
                    logging.warning(f"Non-sister duplicates found for species ({species}) in {self.tree_file}")
                
                self._species_duplicates[species] = [gene_names, False]

    def label (self, node_label_dict, analysis_label_dict, label_descendants = False, analysis_label_symbols = '{}'):

        # Create a copy of the tree to label
        labelled_tree = self._origin_tree.copy()

        # Create a flag to indicate whether the tree has been labelled
        tree_labelled = False

        # Loop through the node labels in the node label dictionary
        for node_label, groups in node_label_dict.items():

            if node_label not in analysis_label_dict:
                logging.warning(f"Node label {node_label} not found among analysis. Skipping for {self.tree_file}")
                continue

            # Assign the corresponding label from analysis_label_dict to node_label
            mapped_node_label = analysis_label_dict[node_label]

            # Creatge a flag to indicate whether the current node label has been applied
            label_applied = False
            
            for node in labelled_tree.traverse('postorder'):
                if node.is_leaf:
                    continue

                node_species_list = [self._gene_to_species_name[_lf] for _lf in node.leaf_names()]
                
                # Create bool to indicate whether non-sister duplicates are found for the current node
                non_sister_duplicates = False

                # Loop through the node_species_list to check for duplicate species and assign the status of the duplicate species
                for node_species in node_species_list:

                    # Check if the node_species is a duplicate species
                    if node_species not in self._species_duplicates:
                        continue

                    # Assign the status of the duplicate species
                    is_sister = self._species_duplicates[node_species][1]

                    # Set the non_sister_duplicates flag to True if the duplicate species is not a sister duplicate
                    if not is_sister:
                        non_sister_duplicates = True
                        break

                # If non-sister duplicates are found, skip the current node
                if non_sister_duplicates:
                    logging.warning(f"Node label {node_label} cannot be labelled. Non-sister duplicates found for node {node.name}. Skipping for {self.tree_file}.")
                    continue

                # Check that each group in the label_dict is a subset of the node species
                if any(set(group).isdisjoint(node_species_list) for group in groups):
                    continue

                # If we reach this point, it means that the node_species contains at least one member of each group
                group_members = [item for group in groups for item in group]

                if not set(node_species_list).issubset(group_members):
                    continue

                # node_label and continue to the next node label
                node.name = f"{analysis_label_symbols[0]}{mapped_node_label}{analysis_label_symbols[1]}"
                tree_labelled = True
                label_applied = True

                # If the--label-descendants option is not used, break the loop and continue to the next node label
                if not label_descendants:
                    break

                # If the --label-descendants option is used, label all descendants of the current node with the same label
                for child in node.descendants():
                    if child.is_leaf:
                        child.name += f"{analysis_label_symbols[0]}{mapped_node_label}{analysis_label_symbols[1]}"
                    else:
                        child.name = f"{analysis_label_symbols[0]}{mapped_node_label}{analysis_label_symbols[1]}"

                # If the --label-descendants option is used, we can assign the labelled tree to self._labelled_tree and return early
                self._labelled_tree = labelled_tree
                return
            
            if not label_applied:
                logging.warning(f"Node label {node_label} cannot be labelled. Skipping for {self.tree_file}.")

        # Label the non-duplicate leaves with the corresponding label fr analysis_label_dictom
        for leaf in labelled_tree:
            leaf_species = self._gene_to_species_name[leaf.name]
            if leaf_species not in self._species_uniques:
                continue
            if leaf_species not in analysis_label_dict:
                continue
            mapped_node_label = analysis_label_dict[leaf_species]
            leaf.name += f"{analysis_label_symbols[0]}{mapped_node_label}{analysis_label_symbols[1]}"
            tree_labelled = True
            label_applied = True

        # Label the duplicate leaves with the corresponding label from analysis_label_dict
        for species, species_duplicate in self._species_duplicates.items():
            duplicate_leaves = species_duplicate[0]
            if species not in analysis_label_dict:
                continue
            duplicate_node = labelled_tree.common_ancestor(duplicate_leaves)
            if duplicate_leaves != set(duplicate_node.leaf_names()):
                logging.warning(f"Node label {species} cannot be labelled. Non-sister duplicates found. Skipping for {self.tree_file}.")
                continue
            mapped_node_label = analysis_label_dict[species]
            duplicate_node.name = f"{analysis_label_symbols[0]}{mapped_node_label}{analysis_label_symbols[1]}"
            tree_labelled = True
            label_applied = True

        if tree_labelled:
            self._labelled_tree = labelled_tree
        else:
            self._labelled_tree = None

    def writeLabelledTree (self, output_dir, analysis_name, output_naming_style, parser=1):

        # Confirm the file type is newick, or raise an error
        if self._file_type != 'newick':
            raise ValueError(f"Cannot write labelled tree for {self.tree_file}. File type is {self._file_type}, expected 'newick'.")

        # Confirm that there is a labelled tree to write
        if self._labelled_tree is not None:

            # Update the output directory based if necessary
            if output_naming_style == 'subdirectory':
                output_dir = os.path.join(output_dir, analysis_name)
            
            # Create the output directory
            os.makedirs(output_dir, exist_ok=True)

            # Determine the output filename based on the output naming style
            if output_naming_style == 'suffix':
                output_filename_wo_ext = f"{os.path.splitext(self.tree_file)[0]}.{analysis_name}"
            else:
                output_filename_wo_ext = f"{os.path.splitext(self.tree_file)[0]}"

            # Write the labelled tree to the output directory
            self._labelled_tree.write(outfile = os.path.join(output_dir, f"{output_filename_wo_ext}.tre"), parser=parser)

            # Open file to write the labelled tree as an ASCII string
            with open(os.path.join(output_dir, f"{output_filename_wo_ext}.fig"), 'w') as ascii_file:
              ascii_file.write(self._labelled_tree.to_str(props=['name']).replace('⊗', '─'))
        else:
            logging.error(f"No labelled tree to write for {self.tree_file}")
    
    def printTree (self):
        
        print(self._origin_tree.to_str(props=["name"]))

def parseABSREL():
    absrel_parser = argparse.ArgumentParser(description="Extract the tree from an aBSREL JSON.")
    absrel_parser.add_argument('--input-json', help="Path to the aBSREL JSON file", required=True)
 
    return vars(absrel_parser.parse_args())

 
def main():

    # Parse command line arguments
    absrel_args = parseABSREL()

    test= GeneTree.fromHyphy(absrel_args['input_json'])
    test.printTree()

    

    '''# Read the JSON file
    with open(absrel_args['input_json']) as absrel_file:
        absrel_data = json.load(absrel_file)
 
    gene_newick = get_tree(absrel_data)

    gene_tree = ete4.Tree(gene_newick, parser=1)
    
    # Create a dictionary to hold the node names and their corresponding leaf names
    gene_node_to_leaves = defaultdict(list)

    for node in gene_tree.traverse("postorder"):
        if not node.is_leaf:
            gene_node_to_leaves[node.name].extend(leaf for leaf in node.leaf_names())

    print("\nNode to Leaves Mapping:")
    for node_name, leaves in gene_node_to_leaves.items():
        print(f"{node_name}: {leaves}")'''
 
if __name__ == "__main__":
    main()