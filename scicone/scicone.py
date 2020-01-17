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

def generate_bp_roc(inferred_bps_dict, true_bps_indicator, threshold_coeffs, window_size, correct_close=True, add_dummy_bps=True):
    bps = pd.DataFrame(inferred_bps_dict['all_bps_comparison'])
    bps.columns = ['idx','log_sp','range']
    bps.sort_values('idx')['idx'].tolist()
    bps.index = bps['idx']

    bps['ranking'] = bps['log_sp'] / bps['range']
    bps = bps.dropna()
    # threshold_coeffs = sorted(bps['ranking'].values)

    # correcting for the bps 1-2 bins nearby
    if correct_close:
        for index, row in bps.iterrows():
            idx_val = bps.loc[index, 'idx']
            for gt in np.where(true_bps_indicator)[0]:
                if (abs(idx_val - gt) <=np.ceil(window_size/2) and idx_val != gt):
                    # print('correcting ' + str(idx_val) + '->' + str(gt))
                    bps.loc[index,'idx'] = gt

    # Add remaining bins to bps to make sure all bins leads to TPR and FPR == 1
    if add_dummy_bps:
        new_bps_indicator = np.ones(true_bps_indicator.shape[0])
        new_bps_indicator[inferred_bps_dict['all_bps_comparison'][:,0].astype(int)] += 1
        dummy_bps = np.where(new_bps_indicator==1)[0]

        dummy_bps_df = pd.DataFrame({'idx':dummy_bps, 'log_sp':np.zeros(dummy_bps.shape[0]), 'range':np.zeros(dummy_bps.shape[0]), 'ranking':np.zeros(dummy_bps.shape[0])})
        dummy_bps_df.index = dummy_bps_df['idx']
        dummy_bps_df['ranking'] = 1e-6

        bps = pd.concat([bps, dummy_bps_df])

    tpr_values = []
    fpr_values = []
    roc_values = []
    bps_indicators = []
    for thr in threshold_coeffs:
        inferred_bps = []
        for index, row in bps.iterrows():
            if row['ranking'] > thr:
                inferred_bps.append(row['idx'])
            else:
                break

        inferred_bps_indicator = np.zeros(true_bps_indicator.shape[0])
        inferred_bps_indicator[np.array(inferred_bps).astype(int)] = 1.
        inferred_bps_indicator = np.array(inferred_bps_indicator)

        tp = np.count_nonzero(inferred_bps_indicator[np.where(true_bps_indicator==1)[0]]==1)
        fp = np.count_nonzero(inferred_bps_indicator[np.where(true_bps_indicator==0)[0]]==1)
        tn = np.count_nonzero(inferred_bps_indicator[np.where(true_bps_indicator==0)[0]]==0)
        fn = np.count_nonzero(inferred_bps_indicator[np.where(true_bps_indicator==1)[0]]==0)

        tpr = tp / (tp + fn)
        fpr = fp / (fp + tn)

        tpr_values.append(tpr)
        fpr_values.append(fpr)
        bps_indicators.append(inferred_bps_indicator)

    roc_curve = dict()
    roc_curve['tpr'] = tpr_values
    roc_curve['fpr'] = fpr_values
    roc_curve['bps'] = bps_indicators
    roc_curve['auc'] = auc(fpr_values, tpr_values)

    return roc_curve

import pandas as pd
import numpy as np
from collections import Counter
from sklearn.metrics import roc_curve, auc
from scipy.cluster.hierarchy import ward, leaves_list
from scipy.spatial.distance import pdist
import seaborn as sns
import matplotlib.pyplot as plt

new_sci = SCICoNE('/cluster/work/bewi/members/pedrof/sc-dna/build/',
                  '/cluster/work/bewi/members/pedrof/sc-dna/notebooks/', persistence=False)
old_sci = SCICoNE('/cluster/work/bewi/ngs/projects/tumorProfiler/code/dna-pipeline/sc-dna/bin/',
                  '/cluster/work/bewi/members/pedrof/sc-dna/notebooks/', persistence=False)

data = new_sci.simulate_data(n_cells=200, n_nodes=10, n_bins=1000, n_regions=40, n_reads=10000, nu=10.0, max_regions_per_node=2, min_reg_size=20)
sim_tree = data['tree']
sim_tree.plot_tree()

cell_genotypes = pd.DataFrame(data['ground_truth'])
cell_bps = cell_genotypes.diff(periods=1, axis=1)
cell_bps = cell_bps.fillna(value=0.0)
cell_bps[cell_bps != 0] = 1 # replace the non-zeroes by 1
grouped_cell_bps = cell_bps.sum(axis=0)
ground_truth = grouped_cell_bps[grouped_cell_bps > 0]
ground_truth = ground_truth.index.tolist()

bps_indicator = np.zeros(grouped_cell_bps.shape[0])
bps_indicator[grouped_cell_bps>=1] = 1
print(f"Number of true breakpoints: {np.sum(bps_indicator==1)}")

counts = data['d_mat']
Z = ward(pdist(counts))
hclust_index = leaves_list(Z)
counts = counts[hclust_index]
cell_genotypes = cell_genotypes.to_numpy()
cell_genotypes = cell_genotypes[hclust_index]

fig = plt.figure(figsize=(16, 4))
cmap = sns.diverging_palette(220, 10, as_cmap=True)
plt.pcolor(counts, cmap=cmap)
plt.colorbar()
ax = plt.gca()
plt.title('Raw counts with breakpoints')
plt.ylabel('cell')
ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), colors='black', linestyles='dashed', linewidths=2)
plt.show()

fig = plt.figure(figsize=(16, 4))
cmap = sns.diverging_palette(220, 10, as_cmap=True)
plt.pcolor(cell_genotypes, cmap=cmap)
plt.colorbar()
ax = plt.gca()
plt.title('True copy numbers with breakpoints')
plt.xlabel('bin')
plt.ylabel('cell')
ax.vlines(np.where(bps_indicator)[0], *ax.get_ylim(), colors='black', linestyles='dashed', linewidths=2)
plt.show()

# True number of clusters
np.unique(cell_genotypes, return_counts=False, axis=0).shape[0]

# Global thresholds to try
threshold_coeffs = np.linspace(0.0,20.0, 100) # 100 thresholds btw 1 and 16
n_bins = cell_genotypes.shape[1]

# New method
new_bps = new_sci.detect_breakpoints(data['d_mat'], window_size=10, threshold=0.1)
new_roc_curve = generate_bp_roc(new_bps, bps_indicator, threshold_coeffs, 10, correct_close=True, add_dummy_bps=True)

fig = plt.figure(figsize=(16, 12))
ax = plt.subplot(3, 1, 1)
plt.pcolor(counts, cmap=cmap)
# plt.colorbar()
ax = plt.gca()
plt.title('True copy numbers with inferred breakpoints')
plt.xlabel('bin')
plt.ylabel('cell')
ax.vlines(np.where(new_roc_curve['bps'][0]), *ax.get_ylim(), colors='black', linestyles='dashed', linewidths=2)

ax = plt.subplot(3, 1, 2, sharex=ax)
cmap = sns.diverging_palette(220, 10, as_cmap=True)
plt.pcolor(cell_genotypes, cmap=cmap)
# plt.colorbar()
ax = plt.gca()
plt.title('True copy numbers with inferred breakpoints')
plt.xlabel('bin')
plt.ylabel('cell')
ax.vlines(np.where(new_roc_curve['bps'][0]), *ax.get_ylim(), colors='black', linestyles='dashed', linewidths=2)

plt.subplot(3, 1, 3, sharex=ax)
plt.plot(np.log(new_bps['sp_vec']))
ax = plt.gca()
ax.vlines(np.where(new_roc_curve['bps'][0]), *ax.get_ylim(), colors='black', linestyles='dashed', linewidths=2)
plt.show()

# Old method
old_bps = old_sci.detect_breakpoints(data['d_mat'], window_size=10, threshold=0.1)
old_roc_curve = generate_bp_roc(old_bps, bps_indicator, threshold_coeffs, 10, correct_close=True, add_dummy_bps=True)

fig = plt.figure(figsize=(16, 12))
ax = plt.subplot(3, 1, 1)
plt.pcolor(counts, cmap=cmap)
# plt.colorbar()
ax = plt.gca()
plt.title('True copy numbers with inferred breakpoints')
plt.xlabel('bin')
plt.ylabel('cell')
ax.vlines(np.where(old_roc_curve['bps'][0]), *ax.get_ylim(), colors='black', linestyles='dashed', linewidths=2)

ax = plt.subplot(3, 1, 2, sharex=ax)
cmap = sns.diverging_palette(220, 10, as_cmap=True)
plt.pcolor(cell_genotypes, cmap=cmap)
# plt.colorbar()
ax = plt.gca()
plt.title('True copy numbers with inferred breakpoints')
plt.xlabel('bin')
plt.ylabel('cell')
ax.vlines(np.where(old_roc_curve['bps'][0]), *ax.get_ylim(), colors='black', linestyles='dashed', linewidths=2)

plt.subplot(3, 1, 3, sharex=ax)
plt.plot(np.log(old_bps['sp_vec']))
ax = plt.gca()
ax.vlines(np.where(old_roc_curve['bps'][0]), *ax.get_ylim(), colors='black', linestyles='dashed', linewidths=2)
plt.show()


# Plot the ROC of each method
plt.figure(figsize=(8,8))
plt.plot(new_roc_curve['fpr'], new_roc_curve['tpr'], color="darkorange", label='New ROC curve (area = %0.2f)' % auc(new_roc_curve['fpr'], new_roc_curve['tpr']), marker='.', alpha=0.4)
plt.plot(old_roc_curve['fpr'], old_roc_curve['tpr'], color="navy", label='Old ROC curve (area = %0.2f)' % auc(old_roc_curve['fpr'], old_roc_curve['tpr']), marker='.', alpha=0.4)
plt.plot([0, 1], [0, 1], color='gray', linestyle='--')
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title('Receiver operating characteristic')
plt.legend(loc="lower right")
plt.show()

fig = plt.figure(figsize=(20,4))
plt.plot(np.log(old_bps['sp_vec']))
plt.plot(np.log(new_bps['sp_vec']))
plt.xticks(ticks=np.linspace(0, 1000, 500).astype(int), labels='')
plt.show()

plt.plot(threshold_coeffs, new_roc_curve['tpr'], marker='.', label='New')
plt.plot(threshold_coeffs, old_roc_curve['tpr'], marker='.', label='Old')
plt.legend()
plt.xlabel('threshold')
plt.title('TPR')

plt.plot(threshold_coeffs, new_roc_curve['fpr'], marker='.', label='New')
plt.plot(threshold_coeffs, old_roc_curve['fpr'], marker='.', label='Old')
plt.legend()
plt.xlabel('threshold')
plt.title('FPR')

# Multiple runs
n_reps = 10
n_cells = 200
n_nodes = 10
n_bins = 1000
n_regions = 40
n_reads = 10000
nu = 1.0
window_size = 10
threshold_coeffs = np.linspace(0.0, 20.0, 100)
roc_curve_lists = dict()
roc_curve_lists['new'] = []
roc_curve_lists['old'] = []
for rep in range(n_reps):
    print(f"Rep. number {rep+1}")

    # Generate data
    data = new_sci.simulate_data(n_cells=n_cells, n_nodes=n_nodes, n_bins=n_bins, n_regions=n_regions,
                            n_reads=n_reads, nu=nu, max_regions_per_node=2, min_reg_size=2*window_size)

    # Get true breakpoints
    cell_genotypes = pd.DataFrame(data['ground_truth'])
    cell_bps = cell_genotypes.diff(periods=1, axis=1)
    cell_bps = cell_bps.fillna(value=0.0)
    cell_bps[cell_bps != 0] = 1 # replace the non-zeroes by 1
    grouped_cell_bps = cell_bps.sum(axis=0)
    ground_truth = grouped_cell_bps[grouped_cell_bps > 0]
    ground_truth = ground_truth.index.tolist()

    bps_indicator = np.zeros(grouped_cell_bps.shape[0])
    bps_indicator[grouped_cell_bps>=1] = 1

    # Infer breakpoints with new and old methods
    new_bps = new_sci.detect_breakpoints(data['d_mat'], window_size=window_size, threshold=0.1)
    new_roc_curve = generate_bp_roc(new_bps, bps_indicator, threshold_coeffs, window_size, correct_close=True, add_dummy_bps=True)

    old_bps = old_sci.detect_breakpoints(data['d_mat'], window_size=window_size, threshold=0.1)
    old_roc_curve = generate_bp_roc(old_bps, bps_indicator, threshold_coeffs, window_size, correct_close=True, add_dummy_bps=True)

    # Add to list
    roc_curve_lists['new'].append(new_roc_curve)
    roc_curve_lists['old'].append(old_roc_curve)

plt.figure(figsize=(8,8))
old_fprs = []
old_tprs = []
new_fprs = []
new_tprs = []
for old_roc_curve in roc_curve_lists['old']:
    old_fprs.append(old_roc_curve['fpr'])
    old_tprs.append(old_roc_curve['tpr'])
    plt.plot(old_roc_curve['fpr'], old_roc_curve['tpr'], color="navy", alpha=0.4)
for new_roc_curve in roc_curve_lists['new']:
    new_fprs.append(new_roc_curve['fpr'])
    new_tprs.append(new_roc_curve['tpr'])
    plt.plot(new_roc_curve['fpr'], new_roc_curve['tpr'], color="darkorange", alpha=0.4)
old_auc = auc(np.append(1, np.mean(old_fprs, axis=0)), np.append(1, np.mean(old_tprs, axis=0)))
new_auc = auc(np.append(1, np.mean(new_fprs, axis=0)), np.append(1, np.mean(new_tprs, axis=0)))
plt.plot([0, 1], [0, 1], color='gray', linestyle='--')
plt.plot([0.1,0.1], [0.1, 0.1], color="darkorange", label=f'new (avg AUC={new_auc})')
plt.plot([0.1,0.1], [0.1, 0.1], color="navy", label=f'old (avg AUC={old_auc})')
plt.plot([0.1,0.1], [0.1, 0.1], color="white")
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title(f"Receiver operating characteristic (n_nodes={n_nodes})")
plt.legend(loc="lower right")
plt.show()


plt.figure(figsize=(8,8))
plt.plot(np.append(1, np.mean(old_fprs, axis=0)), np.append(1, np.mean(old_tprs, axis=0)), color="navy", alpha=0.8)
plt.plot(np.append(1, np.mean(new_fprs, axis=0)), np.append(1, np.mean(new_tprs, axis=0)), color="darkorange", alpha=0.8)
plt.plot([0, 1], [0, 1], color='gray', linestyle='--')
plt.plot([0.1,0.1], [0.1, 0.1], color="darkorange", label=f'new (avg AUC={new_auc})')
plt.plot([0.1,0.1], [0.1, 0.1], color="navy", label=f'old (avg AUC={old_auc})')
plt.plot([0.1,0.1], [0.1, 0.1], color="white")
plt.gca().set_aspect('equal', adjustable='box')
plt.grid(True)
plt.xlabel('False Positive Rate')
plt.ylabel('True Positive Rate')
plt.title(f"Mean receiver operating characteristic (n_nodes={n_nodes})")
plt.legend(loc="lower right")
plt.show()


counts = data['d_mat']
cell_genotypes = data['ground_truth']
Z = ward(pdist(counts))
hclust_index = leaves_list(Z)
counts = counts[hclust_index]
cell_genotypes = cell_genotypes[hclust_index]

fig = plt.figure(figsize=(16, 12))
ax = plt.subplot(3, 1, 1)
plt.pcolor(counts, cmap=cmap)
# plt.colorbar()
ax = plt.gca()
plt.title('True copy numbers with inferred breakpoints')
plt.xlabel('bin')
plt.ylabel('cell')
ax.vlines(np.where(new_roc_curve['bps'][5]), *ax.get_ylim(), colors='black', linestyles='dashed', linewidths=2)

ax = plt.subplot(3, 1, 2, sharex=ax)
cmap = sns.diverging_palette(220, 10, as_cmap=True)
plt.pcolor(cell_genotypes, cmap=cmap)
# plt.colorbar()
ax = plt.gca()
plt.title('True copy numbers with inferred breakpoints')
plt.xlabel('bin')
plt.ylabel('cell')
ax.vlines(np.where(new_roc_curve['bps'][5]), *ax.get_ylim(), colors='black', linestyles='dashed', linewidths=2)

plt.subplot(3, 1, 3, sharex=ax)
plt.plot(np.log(new_bps['sp_vec']))
ax = plt.gca()
# ax.vlines(np.where(new_roc_curve['bps'][5]), *ax.get_ylim(), colors='black', linestyles='dashed', linewidths=2)
plt.show()


fig = plt.figure(figsize=(16, 12))
ax = plt.subplot(3, 1, 1)
plt.pcolor(counts, cmap=cmap)
# plt.colorbar()
ax = plt.gca()
plt.title('True copy numbers with inferred breakpoints')
plt.xlabel('bin')
plt.ylabel('cell')
ax.vlines(np.where(old_roc_curve['bps'][5]), *ax.get_ylim(), colors='black', linestyles='dashed', linewidths=2)

ax = plt.subplot(3, 1, 2, sharex=ax)
cmap = sns.diverging_palette(220, 10, as_cmap=True)
plt.pcolor(cell_genotypes, cmap=cmap)
# plt.colorbar()
ax = plt.gca()
plt.title('True copy numbers with inferred breakpoints')
plt.xlabel('bin')
plt.ylabel('cell')
ax.vlines(np.where(old_roc_curve['bps'][5]), *ax.get_ylim(), colors='black', linestyles='dashed', linewidths=2)

plt.subplot(3, 1, 3, sharex=ax)
plt.plot(np.log(old_bps['sp_vec']))
ax = plt.gca()
# ax.vlines(np.where(new_roc_curve['bps'][5]), *ax.get_ylim(), colors='black', linestyles='dashed', linewidths=2)
plt.show()
