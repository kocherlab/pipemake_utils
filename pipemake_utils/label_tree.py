import os
import logging
import argparse

from ete4 import Tree
from collections import defaultdict

from pipemake_utils.misc import confirmDir, confirmFile
from pipemake_utils.logger import startLogger, logArgDict

def returnDuplicateSpecies (tree):
    species_list = [leaf.split('-')[0] for leaf in tree.leaf_names()]
    dup_species = set([s for s in species_list if species_list.count(s) > 1])
    return dup_species

def yieldDuplicateLeafNames (tree):

    # Iterate through the duplicate species and yield the species name and the set of duplicate leaves
    for species in returnDuplicateSpecies(tree):
        duplicate_leaves = set()
        for leaf in tree.leaf_names():
            if leaf.startswith(species + '-'):
                duplicate_leaves.add(leaf)
        yield species, duplicate_leaves

def yieldNonDuplicateLeaves (tree):
    # Assign the duplicate species
    dup_species = returnDuplicateSpecies(tree)

    # Iterate through the leaves of the tree
    for leaf in tree:
        species = leaf.name.split('-')[0]
        if species not in dup_species:
            yield leaf

def labelTreeNodes (tree, node_label_dict, analysis_label_dict, analysis_label_symbols = '{}'):

    # Create a copy of the tree to label
    labelled_tree = tree.copy()

    # Create a flag to indicate whether the tree has been labelled
    tree_labelled = False

    for node_label, groups in node_label_dict.items():

        if node_label not in analysis_label_dict:
            print(f"Warning: node_label {node_label} not found in analysis_label_dict. Skipping this label.")
            continue

        # Assign the corresponding label from analysis_label_dict to node_label
        mapped_node_label = analysis_label_dict[node_label]

        # Creatge a flag to indicate whether the current node label has been applied
        label_applied = False
        
        for node in labelled_tree.traverse('postorder'):
            if not node.is_leaf:
                node_species = [_lf.split('-')[0] for _lf in node.leaf_names()]

                # Check that each group in the label_dict is a subset of the node species
                if any(set(group).isdisjoint(node_species) for group in groups):
                    continue

                # If we reach this point, it means that the node_species contains at least one member of each group
                group_members = [item for group in groups for item in group]
                
                if not set(node_species).issubset(group_members):
                    continue

                # node_label and continue to the next node label
                node.name = f"{analysis_label_symbols[0]}{mapped_node_label}{analysis_label_symbols[1]}"
                tree_labelled = True
                label_applied = True
                break
        
        if not label_applied:
            print(f"Warning: node label {node_label} could not be applied to any node in the tree.")

    # Label the non-duplicate leaves with the corresponding label fr analysis_label_dictom
    for leaf in yieldNonDuplicateLeaves(labelled_tree):
        leaf_species = leaf.name.split('-')[0]
        if leaf_species not in analysis_label_dict:
            continue
        mapped_node_label = analysis_label_dict[leaf_species]
        leaf.name += f"{analysis_label_symbols[0]}{mapped_node_label}{analysis_label_symbols[1]}"
        tree_labelled = True
        label_applied = True

    # Label the duplicate leaves with the corresponding label from analysis_label_dict
    for species, duplicate_leaves in yieldDuplicateLeafNames(labelled_tree):
        duplicate_node = labelled_tree.common_ancestor(duplicate_leaves)
        mapped_node_label = analysis_label_dict[species]
        duplicate_node.name = f"{analysis_label_symbols[0]}{mapped_node_label}{analysis_label_symbols[1]}"
        tree_labelled = True
        label_applied = True

    if tree_labelled:
        return labelled_tree
    else:
        return Tree()

def confirmTree (tree):
   
    # Iterate through the leaves of the tree
    for species, duplicate_leaves in yieldDuplicateLeafNames(tree):
    
        # Check if the duplcates are in the same clade
        if duplicate_leaves != set(tree.common_ancestor(duplicate_leaves).leaf_names()):
            print(f"Error: Non-sister duplicates found for species {species}")
            return False
    
    return True

def createLabelDict (species_tree, species_label_symbols='[]'):
    label_dict = defaultdict(list)

    logging.info("Creating label dictionary from species tree")

    for node in species_tree.traverse():
        if not node.is_leaf and node.name:
            node_children = []
            for child in node.children:
                if child.is_leaf:
                    node_children.append([child.name])
                else:
                    node_children.append(list(child.leaf_names()))
            
            # Remove the first and last character of the node name if they are in species_label_symbols
            if node.name[0] in species_label_symbols and node.name[-1] in species_label_symbols:
                node.name = node.name[1:-1]
            elif node.name[0] in species_label_symbols or node.name[-1] in species_label_symbols:
                logging.warning(f"Node name '{node.name}' has mismatched label symbols.")
            
            label_dict[node.name] = node_children

            logging.info(f"Added label '{node.name}' with species: {node_children}")
    
    logging.info(f"Label dictionary created with {len(label_dict)} labels.")

    return dict(sorted(label_dict.items(), key=lambda x: sum(len(sub) for sub in x[1])))

def analysisLabelFile (label_filename):

    logging.info(f"Reading analysis label file from {label_filename}")
    
    # Create a dictionary to hold the analysis labels
    analysis_labels = defaultdict(lambda: defaultdict(str))

    # Read the analysis label file and create a dictionary mapping node labels to analysis labels
    with open(label_filename, 'r') as label_file:
        for label_str in label_file:

            # Skip empty lines
            if not label_str.strip():
                continue

            # Raise an error if the line does not contain a ':' character
            if ':' not in label_str:
                raise IOError(f"Warning: Line '{label_str}' does not contain a ':' character. Please check the format of the label file.")
            
            # Assign the analysis name, and data
            analysis_name, analysis_data = label_str.split()

            logging.info(f"Processing analysis label: {analysis_name.strip()}")

            # Process the analysis label data
            for analysis_label_data in analysis_data.split(';'):
                analysis_label, analysis_nodes_str = analysis_label_data.split(':')
                for analysis_node in analysis_nodes_str.split(','):
                    analysis_labels[analysis_name.strip()][analysis_node.strip()]= analysis_label.strip()

                    logging.info(f"Mapping node '{analysis_node.strip()}' to label '{analysis_label.strip()}'")

    logging.info(f"Analysis label dictionary created with {len(analysis_labels)} analyses.")

    return analysis_labels

def treeParser ():
    tree_parser = argparse.ArgumentParser(description='Label gene trees based on a species tree and an analysis label file')
    tree_parser.add_argument('--species-tree', help = 'Newick string of the species tree with node labels', type=str, required = True, action = confirmFile())
    tree_parser.add_argument('--analysis-label-file', help = 'String of the analysis label file in the format "ANALYSIS_LABEL:NODE1,NODE2;ANALYSIS_LABEL2:NODE3,NODE4"', type=str, required = True, action = confirmFile())
    tree_parser.add_argument('--gene-tree-dir', help = 'Directory containing the gene trees to be labelled', type=str, required = True, action = confirmDir())
    tree_parser.add_argument('--output-dir', help = 'Output directory for the labelled gene trees', type=str, default = 'labelled_trees')
    return vars(tree_parser.parse_args())

def main ():
    # Assign the arguments from the command-line
    tree_args = treeParser()

    # Start the logger and log the arguments
    startLogger(f'{tree_args["output_dir"]}.log')
    logArgDict(tree_args)

    # Read in the species tree
    species_tree = Tree(tree_args['species_tree'], parser=1)
    logging.info(f"Species tree read from {tree_args['species_tree']}")

    # Create a dictionary mapping node labels to the corresponding groups of species
    node_label_dict = createLabelDict(species_tree)

    # Create the analysis label dictionary from the analysis label file
    analysis_labels = analysisLabelFile(tree_args['analysis_label_file'])
    
    # Loop the analysis labels and label the gene trees for each analysis
    for analysis_name, analysis_label_dict in analysis_labels.items():

        logging.info(f"Labelling gene trees for analysis '{analysis_name}'")

        # Loop the gene trees in the specified directory
        for gene_tree_file in os.listdir(tree_args['gene_tree_dir']):

            # Skip files that do not end with '.tre'
            if not gene_tree_file.endswith('.tre'):
                continue

            # Read in the gene tree
            gene_tree = Tree(os.path.join(tree_args['gene_tree_dir'], gene_tree_file), parser=1)

            # Confirm that the gene tree is valid and can be labelled
            if not confirmTree(gene_tree):
                logging.error(f"Tree '{gene_tree_file}' failed the confirmation check. Skipping this tree.")
                continue

            # Label the gene tree
            labelled_gene_tree = labelTreeNodes(gene_tree, node_label_dict, analysis_label_dict)

            # Create the output directory for the labelled gene trees if it does not exist
            output_dir = os.path.join(tree_args['output_dir'], analysis_name)
            os.makedirs(output_dir, exist_ok=True)

            # Write the labelled gene tree to the output directory
            labelled_gene_tree.write(outfile=os.path.join(output_dir, gene_tree_file), parser=1)

            logging.info(f"Labelled tree for '{gene_tree_file}'")

if __name__ == '__main__':
    main()