import numpy as np
from copy import deepcopy

def sort_chromosomes(chromosome_list):
    """
    Sorts a list of unordered chromosome names
    :param chromosome_list: list of unordered characters denoting chromosomes '1', '2', ..., 'X', 'Y'
    """
    # Replace X and Y with 23 and 24
    sorted_chromosome_list = np.array(chromosome_list)
    sorted_chromosome_list[np.where(sorted_chromosome_list == "X")[0]] = 23
    sorted_chromosome_list[np.where(sorted_chromosome_list == "Y")[0]] = 24

    # Convert everything to integer
    sorted_chromosome_list = sorted_chromosome_list.astype(int)

    # Sort
    sorted_chromosome_list = np.sort(sorted_chromosome_list)

    # Convert back to string
    sorted_chromosome_list = sorted_chromosome_list.astype(str)
    sorted_chromosome_list[np.where(sorted_chromosome_list == "23")[0]] = "X"
    sorted_chromosome_list[np.where(sorted_chromosome_list == "24")[0]] = "Y"

    return sorted_chromosome_list

def create_fusion_tree(learned_tree, region_neutral_states):
    fusion_tree = deepcopy(learned_tree)

    # Adds 1 node below each node of the learned tree with an added neutral genome
    new_event_dict = dict()
    for i, state in enumerate(region_neutral_states):
        new_event_dict[str(i)]=str(int(state))

    fusion_tree.read_tree_str(learned_tree.tree_str)

    for node in list(fusion_tree.node_dict):
        new_node_id = str(int(node) + 1000)
        fusion_tree.node_dict[new_node_id] = dict(parent_id=node, event_dict=new_event_dict)

    fusion_tree.update_tree_str()

    return fusion_tree
