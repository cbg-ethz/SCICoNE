import os, shutil
import subprocess, re
from snakemake.workflow import Workflow, Rules
import snakemake.workflow
from snakemake import shell
from snakemake.logging import setup_logger
import numpy as np
import pandas as pd
import graphviz
from graphviz import Source

# See https://gist.github.com/marcelm/7e432275dfa5e762b4f8 on how to run a
# Snakemake workflow from Python

class Tree(object):
    def __init__(self):
        self.traces = None
        self.best_scores = None
        self.graphviz_str = ""

    def read_from_file(self, file, output_path=None):
        """
            reads the file containing a tree and converts it to graphviz format
            :param file: path to the tree file.
            :param output_path: (optional) path to the output file.
        """
        with open(file) as f:
            list_tree_file = list(f)

        graphviz_header = [
            "digraph {",
            'node [style=filled,color="#D4C0D6",fontsize=20,margin=0,shape=oval]'
            'edge [arrowhead=none, color="#602A86"]',
        ]

        graphviz_labels = []
        graphviz_links = []

        graphviz_labels.append("0[label=<<font point-size='30'> Neutral </font>>]")  # root

        for line in list_tree_file:
            if line.startswith("node 0:"):
                continue
            elif line.startswith("node"):
                comma_splits = line.split(",")

                comma_first = re.split(" |:", comma_splits[0])
                node_id = comma_first[1]
                p_id = comma_first[4]
                comma_rest = comma_splits[1:]
                comma_rest[0] = comma_rest[0].lstrip("[")
                comma_rest[-1] = comma_rest[-1].rstrip("]\n")
                merged_labels = []
                [k_begin, previous_v] = (int(x) for x in comma_rest[0].split(":"))
                k_end = k_begin
                for term in comma_rest[1:]:  # events vector
                    [k, v] = (int(x) for x in term.split(":"))
                    if k == k_end + 1 and v == previous_v:
                        k_end = k  # update the end
                    else:
                        if k_begin == k_end:
                            merged_labels.append(f"{previous_v:+}R{k_begin}")
                        else:
                            merged_labels.append(f"{previous_v:+}R{k_begin}:{k_end}")
                        k_begin = k_end = k
                    previous_v = v
                # print the last one
                if k_begin == k_end:
                    merged_labels.append(f"{previous_v:+}R{k_begin}")
                else:
                    merged_labels.append(f"{previous_v:+}R{k_begin}:{k_end}")

                str_merged_labels = (
                    "<i> </i>"
                    + " ".join(
                        f"{x}<i> </i><br/>" if i % 10 == 0 and i > 0 else str(x)
                        for i, x in enumerate(merged_labels)
                    )
                    + "<i> </i>"
                )

                newline = "<br/>"
                endline = "<i> </i>"

                # If there are newlines followed by endlines, switch
                str_merged_labels = str_merged_labels.replace(
                    newline + endline, endline + newline
                )

                # Remove whatever comes after the last endline position
                new_end_pos = (
                    [m.start() for m in re.finditer(endline, str_merged_labels)][-1]
                    + len(endline)
                    - 1
                )
                if len(str_merged_labels) > new_end_pos + 1:
                    str_merged_labels = str_merged_labels[: new_end_pos + 1]

                # Replace multiple endlines with one
                while endline * 2 in str_merged_labels:
                    str_merged_labels = str_merged_labels.replace(endline * 2, endline)

                graphviz_labels.append(
                    f"{node_id}[label=<{str_merged_labels}>]"
                )  # use < > to allow HTML
                graphviz_links.append(f"{p_id} -> {node_id}")

        self.graphviz_str = graphviz_header + graphviz_labels + graphviz_links + ["}"]

        if output_path is not None:
            with open(output_path, "w") as file:
                for line in self.graphviz_str:
                    file.write(f"{line}\n")

        self.graphviz_str = '\n'.join(self.graphviz_str)

    def learn_tree(self, data):
        pass

    def plot_tree(self, show_genes=True):
        if self.graphviz_str != "":
            s = Source(self.graphviz_str)
            return s
        else:
            raise Exception("No graphviz string available. Please read from file or learn a new tree from data.")

class SCICoNE(object):
    """
    This class  provides an interface to interact with the outputs from the C++
    program that infers a copy number phylogeny from single-cell WGS data.

    It exposes functions to visualize breakpoint detection results, annotate
    the bins, and plot the resulting trees. Its attributes can be further
    processed using scgenpy, a generic package for pre, post-processing and
    visualization of single-cell copy number data.
    """
    def __init__(self, binary_path, output_path, persistence=False):
        """Create a SCICoNE object.

        binary_path : type
            Path to SCICoNE binaries.
        output_path : type
            Path to SCICoNE output files.
        persistence : boolean
            Wether to delete output files from C++ after loading them into the class.
        """
        self.binary_path = binary_path
        self.simulation_binary = os.path.join(self.binary_path, 'simulation')
        self.bp_binary = os.path.join(self.binary_path, 'breakpoint_detection')
        self.inference_binary = os.path.join(self.binary_path, 'inference')
        self.score_binary = os.path.join(self.binary_path, 'score')

        self.output_path = output_path
        self.persistence = persistence
        self.postfix = 'PYSCICONETEMP'

        self.n_cells = 0
        self.n_nodes = 0
        self.best_tree = None
        self.tree_list = []

    def simulate_data(self, n_cells=200, n_nodes=5, n_bins=1000, n_regions=40, n_reads=10000, nu=1.0, min_reg_size=10, max_regions_per_node=1, ploidy=2, verbosity=0):
        output_path = os.path.join(self.output_path, f"{self.postfix}_simulation")

        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        os.mkdir(output_path)
        cwd = os.getcwd()

        done = False
        while not done:
            try:
                cmd_output = subprocess.run([self.simulation_binary, f"--n_cells={n_cells}", f"--n_nodes={n_nodes}",\
                    f"--n_regions={n_regions}", f"--n_bins={n_bins}", f"--n_reads={n_reads}", f"--nu={nu}",\
                    f"--min_reg_size={min_reg_size}", f"--max_regions_per_node={max_regions_per_node}",\
                    f"--ploidy={ploidy}", f"--verbosity={verbosity}", f"--postfix={self.postfix}"])
            except subprocess.SubprocessError as e:
                print("SubprocessError: ", e.returncode, e.output, e.stdout, e.stderr)

            # The binary generates 4 files: *_d_mat.csv, *_ground_truth.csv, *_region_sizes.txt and *_tree.txt.
            # We read each one of these files into np.array, np.array, np.array and Tree structures, respectively,
            # and output a dictionary with keys "d_mat", "ground_truth", "region_sizes" and "tree".
            d_mat_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_d_mat.csv"
            ground_truth_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_ground_truth.csv"
            region_sizes_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_region_sizes.txt"
            tree_file = f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_tree.txt"

            try:
                # Move to output directory
                cwd = os.getcwd()
                os.rename(os.path.join(cwd, d_mat_file), os.path.join(output_path, d_mat_file))
                os.rename(os.path.join(cwd, ground_truth_file), os.path.join(output_path, ground_truth_file))
                os.rename(os.path.join(cwd, region_sizes_file), os.path.join(output_path, region_sizes_file))
                os.rename(os.path.join(cwd, tree_file), os.path.join(output_path, tree_file))
                done = True
            except OSError as e:
                pass

        d_mat_file = os.path.join(output_path, f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_d_mat.csv")
        ground_truth_file = os.path.join(output_path, f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_ground_truth.csv")
        region_sizes_file = os.path.join(output_path, f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_region_sizes.txt")
        tree_file = os.path.join(output_path, f"{n_nodes}nodes_{n_regions}regions_{n_reads}reads_{self.postfix}_tree.txt")

        try:
            output = dict()
            output['d_mat'] = np.loadtxt(d_mat_file, delimiter=',')
            output['ground_truth'] = np.loadtxt(ground_truth_file, delimiter=',')
            output['region_sizes'] = np.loadtxt(region_sizes_file, delimiter=',')

            tree = Tree()
            tree.read_from_file(tree_file)

            output['tree'] = tree

        except OSError as e:
            print("OSError: ", e.output, e.stdout, e.stderr)

        # Now that we've read all outputs into memory, we delete the temporary files if persistence==False
        if not self.persistence:
            shutil.rmtree(output_path)

        return output

    def detect_breakpoints(self, data, window_size=30, threshold=3.0, bp_limit=300):
        n_cells = data.shape[0]
        n_bins = data.shape[1]
        verbosity = 1 # > 0 to generate all files

        try:
            # Write the data to a file to be read by the binary
            data_file = f"{self.postfix}_bp_detection.txt"
            np.savetxt(data_file, data, delimiter=',')

            cmd_output = subprocess.run([self.bp_binary, f"--d_matrix_file={data_file}", f"--n_bins={n_bins}",\
                f"--n_cells={n_cells}", f"--window_size={window_size}", f"--threshold={threshold}", \
                f"--bp_limit={bp_limit}", f"--verbosity={verbosity}", f"--postfix={self.postfix}"])

            # Delete the data file
            os.remove(data_file)
        except subprocess.SubprocessError as e:
            print("SubprocessError: ", e.returncode, e.output, e.stdout, e.stderr)

        output = dict()

        try:
            cwd = os.getcwd()
            for fn in os.listdir(cwd):
                if self.postfix in fn:
                    key = fn.split('_', 1)[1].split('.')[0]
                    output[key] = np.loadtxt(fn, delimiter=',')
                    os.remove(fn)
        except OSError as e:
            print("OSError: ", e.output, e.stdout, e.stderr)

        return output

    def learn_tree(self, n_reps=10):
        # run tree bin using the complete snakemake workflow

        # if multiple reps:
        # read in all traces
        # pick best tree

        # load best tree

        if not self.persistence:
            # delete files
            pass

    def plot_bps(self, cluster_cells=True):
        pass
