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
