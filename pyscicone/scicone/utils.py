from copy import deepcopy

def create_fusion_tree(learned_tree, region_neutral_states):
    fusion_tree = deepcopy(learned_tree)

    # Adds 1 node below each node of the learned tree with an added neutral genome
    new_event_dict = dict()
    for i, state in enumerate(region_neutral_states):
        new_event_dict[str(i)]=str(int(state))

    for node in list(fusion_tree.node_dict):
        new_node_id = str(int(node) + 1000)
        fusion_tree.node_dict[new_node_id] = dict(parent_id=node, event_dict=new_event_dict)

    fusion_tree.update_tree_str()

    return fusion_tree
