from ete4 import Tree

from collections import defaultdict

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
                print(f"Warning: Node name '{node.name}' has mismatched label symbols.")
            
            label_dict[node.name] = node_children

    return dict(sorted(label_dict.items(), key=lambda x: sum(len(sub) for sub in x[1])))

def analysisLabelFile (label_filename):
    analysis_labels = {}
    #with open(label_filename, 'r') as label_file:
    for label_str in label_filename.splitlines():
        if not label_str.strip():
            continue
        if ':' not in label_str:
            raise IOError(f"Warning: Line '{label_str}' does not contain a ':' character. Please check the format of the label file.")
        for analysis_label_data in label_str.split(';'):
            analysis_label, analysis_nodes_str = analysis_label_data.split(':')
            for analysis_node in analysis_nodes_str.split(','):
                analysis_labels[analysis_node.strip()] = analysis_label.strip()
    return analysis_labels

species_str = '((((1,2)[A],3)[B],4)[C],(5,6)[D])[F];'

species_tree = Tree(species_str, parser=1)

node_label_dict = createLabelDict(species_tree)

label_file = 'TEST:A,B,C,1;REF:D\n'
analysis_label_dict = analysisLabelFile(label_file)


def testTree (tree_str):
    tree = Tree(tree_str, parser=1)
    if not confirmTree(tree):
        return
    labeled_tree = labelTreeNodes(tree, node_label_dict, analysis_label_dict)
    print(labeled_tree.to_str(props=['name'], compact=True))

example_str = '((((1-1,2-1),3-1),4-1),(5-1,6-1));'
testTree(example_str)

example_str = '((((1-1,1-2),3-1),4-1),(5-1,6-1));'
testTree(example_str)

#example_str = '((((1-1,8-1),3-1),4-1),(5-1,6-1));'
#testTree(example_str)

#example_str = '((((1-1,2-1),1-2),4-1),(5-1,6-1));'
#testTree(example_str)
