from ete4 import Tree

from collections import defaultdict

def labelTreeNodes (tree, label_dict):
    for label, groups in label_dict.items():
        for node in tree.traverse('postorder'):
            if not node.is_leaf:
                node_species = [_lf.split('-')[0] for _lf in node.leaf_names()]

                # Check that each group in the label_dict is a subset of the node_species
                if any(set(group).isdisjoint(node_species) for group in groups):
                    continue

                # If we reach this point, it means that the node_species contains at least one member of each group
                group_members = [item for group in groups for item in group]
                
                if not set(node_species).issubset(group_members):
                    continue

                # Label and continue to the next label
                node.name = label
                break

    return tree

def confirmTree (tree):
    # Get a list of the duplicate species in the tree
    species_list = [leaf.split('-')[0] for leaf in tree.leaf_names()]
    dup_species = set([s for s in species_list if species_list.count(s) > 1])
    
    # Iterate through the leaves of the tree
    for species in dup_species:
        duplicate_nodes = set()
        for leaf in tree.leaf_names():
            if leaf.startswith(species + '-'):
                duplicate_nodes.add(leaf)
        
        # Check if the duplcates are in the same clade
        if duplicate_nodes != set(tree.common_ancestor(duplicate_nodes).leaf_names()):
            print(f"Error: Non-sister duplicates found for species {species}")
            return False
    
    return True

def createLabelDict (species_tree):
    label_dict = defaultdict(list)

    for node in species_tree.traverse():
        if not node.is_leaf and node.name:
            node_children = []
            for child in node.children:
                if child.is_leaf:
                    node_children.append([child.name])
                else:
                    node_children.append(list(child.leaf_names()))
            label_dict[node.name] = node_children

    return dict(sorted(label_dict.items(), key=lambda x: sum(len(sub) for sub in x[1])))

species_str = '((((1,2)[A],3)[B],4)[C],(5,6)[D])[F];'

species_tree = Tree(species_str, parser=1)

label_dict = createLabelDict(species_tree)

def testTree (tree_str):
    tree = Tree(tree_str, parser=1)
    if not confirmTree(tree):
        return
    labeled_tree = labelTreeNodes(tree, label_dict)
    print(labeled_tree.to_str(props=['name'], compact=True))

example_str = '((((1-1,2-1),3-1),4-1),(5-1,6-1));'
testTree(example_str)

example_str = '((((1-1,1-2),3-1),4-1),(5-1,6-1));'
testTree(example_str)

example_str = '((((1-1,8-1),3-1),4-1),(5-1,6-1));'
testTree(example_str)

example_str = '((((1-1,2-1),1-2),4-1),(5-1,6-1));'
testTree(example_str)
