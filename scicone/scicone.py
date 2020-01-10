import os
import subproces
from snakemake.workflow import Workflow, Rules
import snakemake.workflow
from snakemake import shell
from snakemake.logging import setup_logger
import numpy as np
import pandas as pd
import graphviz

# See https://gist.github.com/marcelm/7e432275dfa5e762b4f8 on how to run a
# Snakemake workflow from Python

class Tree(object):
    def __init__(self):
        self.traces\ = None
        self.best_scores = None
        self.

    def read_tree_from_file(self):
        pass

    def learn_tree(self, data):
        pass

    def plot_tree(self, show_genes=True):
        pass


class SCICoNE(object):
    """
    This class  provides an interface to interact with the outputs from the C++
    program that infers a copy number phylogeny from single-cell WGS data.

    It exposes functions to visualize breakpoint detection results, annotate
    the bins, and plot the resulting trees. Its attributes can be further
    processed using scgenpy, a generic package for pre, post-processing and
    visualization of single-cell copy number data.
    """
    def __init__(self, tree_output_path, bp_output_path=None, persistence=False):
        """Create a SCICoNE object.

        tree_results_path : type
            Description of parameter `tree_results_path`.
        bp_results_path : type
            Description of parameter `bp_results_path`.
        persistence : boolean
            Wether to delete output files from C++ after loading them into the class.
        """
        self.tree_results_path = tree_results_path
        self.bp_results_path = bp_results_path
        self.persistence = persistence

        self.n_cells = 0
        self.n_nodes = 0
        self.best_tree = None
        self.tree_list = []


    def simulate_data(self, n_cells=500, n_nodes=50, n_regions=50, n_reads=10000, ploidy=2):
        {params.sim_bin} --n_regions {wildcards.regions} --n_reads {wildcards.reads} --n_iters {params.n_iters} --n_cells {params.n_cells} --n_bins {params.n_bins} --n_nodes \
        {params.n_nodes} --verbosity 0 --ploidy 2 --postfix {wildcards.rep_id};
        try:
            cmd_output = subprocess.run([params.binary, f"--d_matrix_file={input.d_matrix_file}", f"--n_bins={n_bins}",\
                f"--n_cells={n_cells}", f"--window_size={params.window_size}", f"--postfix={params.posfix}",\
                f"--verbosity={params.verbosity}", f"--threshold={params.threshold}", f"--bp_limit={params.bp_limit}"])
        except subprocess.SubprocessError as e:
            print("Status : FAIL", e.returncode, e.output, e.stdout, e.stderr)

    def detect_breakpoints(self):
        pass

    def learn_tree(self, n_reps=10):
        # run tree bin using the complete snakemake workflow

        # if multiple reps:
        # read in all traces
        # pick best tree

        # load best tree
        self.tree = _read_tree_from_file()

        if not self.persistence:
            # delete files
            pass

    def plot_bps(self, cluster_cells=True):
        pass
