from ComputeDistance import *
from copy import *
from mdmst import *

class Tree:
    def __init__(self, name, out_edge, in_edge):
        self.name = name  # name of node
        self.out_edge = out_edge  # dictionary of node name - distance pair
        self.in_edge = in_edge  # dictionary of node name - distance pair

    def get_out_degree(self):
        return len(self.out_edge)

    def get_in_degree(self):
        return len(self.in_edge)


def create_tree(nodes, node_name_list,root):
    tree_node_dict = {}
    for node in node_name_list:
        temp_out_edge = {}
        temp_in_edge = {}
        for other_node in node_name_list:
            if not node == other_node:
                temp_out_edge[other_node] = dist(nodes[node], nodes[other_node])
                temp_in_edge[other_node] = dist(nodes[other_node], nodes[node])
        tree_node_dict[node] = temp_out_edge
    return tree_node_dict
