"""
python {params.script} --bin_path {params.bin_path} --counts {input.d_mat} --cnvs {input.init_cnvs} --intermediate_dir {scripts.intermediate_dir} --seed {params.seed} --out {output.conet_inferred_cnvs}
"""

def tree_to_newick(g, root=None):
    if root is None:
        roots = list(filter(lambda p: p[1] == 0, g.in_degree()))
        assert 1 == len(roots)
        root = roots[0][0]
    subgs = []
    for child in g[root]:
        if len(g[child]) > 0:
            subgs.append(tree_to_newick(g, root=child))
        else:
            subgs.append(child)
    if root == "{}":
        root = "ROOT"
    return "(" + ','.join(subgs) + f"){root}"

import pandas as pd
import numpy as np
import argparse
from conet.data_converter.corrected_counts import CorrectedCounts
from conet.data_converter.data_converter import DataConverter
from conet import CONET, CONETParameters, InferenceResult
import os
from pathlib import Path
import shutil

parser = argparse.ArgumentParser()
parser.add_argument('--bin_path', type=str, required=True)
parser.add_argument('--counts', required=True, type=str)
parser.add_argument('--cnvs', required=True, type=str)
parser.add_argument('--intermediate_dir', required=True, type=str)
parser.add_argument('--delim', type=str, default=" ")
parser.add_argument('--seed', type=int, default=42)
parser.add_argument('--em_iters', type=int, default=100000)
parser.add_argument('--pt_iters', type=int, default=200000)
parser.add_argument('--out_cnvs', required=True, type=str)
parser.add_argument('--out_tree', required=True, type=str)
parser.add_argument('--out_attachments', required=True, type=str)

if __name__ == "__main__":

    args = parser.parse_args()

    bin_path = args.bin_path
    input_mat = args.counts
    input_cnvs = args.cnvs
    intermediate_dir = args.intermediate_dir
    delim = args.delim
    seed = args.seed
    em_iters = args.em_iters
    pt_iters = args.pt_iters
    out_cnvs = args.out_cnvs
    out_tree = args.out_tree
    out_attachments = args.out_attachments

    # os.mkdir(intermediate_dir)
    Path(intermediate_dir).mkdir(parents=True, exist_ok=True)

    # Create corrected counts file
    raw_counts = np.loadtxt(input_mat, delimiter=',')
    n_cells = raw_counts.shape[0]
    n_bins = raw_counts.shape[1]

    cnvs = np.loadtxt(input_cnvs, delimiter=args.delim)

    assert cnvs.shape == raw_counts.shape

    # Get candidate breakpoints
    bps = np.where(np.any(np.diff(cnvs, axis=1) != 0, axis=0))[0] #+ 1
    bps = np.sort(np.concatenate([bps, bps+1]))
    print(f"Got {len(bps)} candidate breakpoints")

    # Normalize to neutral ploidy 2
    m = np.mean(raw_counts,axis=1)
    normalized_counts = raw_counts / m[:,np.newaxis] * 2

    # Make dataframe
    add = np.full([n_bins, 5], 1.0, dtype=np.int64)
    add[:, 1] = range(0, n_bins) # bin start
    add[:, 2] = range(1, n_bins + 1) # bin end
    add[:, 4] = 0 # breakpoint markers
    add[bps, 4] = 1

    data = np.hstack([add, normalized_counts.T])
    df = pd.DataFrame(data, columns=["chr", "start", "end", "width", "candidate_brkp"] + ["cell" + str(i) for i in range(0, n_cells)])
    df["chr"] = 1
    df["start"] = range(0, n_bins)
    df["end"] = range(1, n_bins+1)
    df["width"] = 1
    df["candidate_brkp"] = 0
    df.loc[bps,"candidate_brkp"] = 1

    # Run CONET
    cc = CorrectedCounts(df)
    #cc.add_chromosome_ends(neutral_cn=args.neutral_cn, end_bin_length=150000)
    DataConverter(event_length_normalizer=3095677412).create_CoNET_input_files(out_path=intermediate_dir, corrected_counts=cc)
    conet = CONET(bin_path, intermediate_dir)
    params = CONETParameters(
        data_dir=intermediate_dir,
        param_inf_iters=em_iters,
        pt_inf_iters=pt_iters,
        counts_penalty_s1=100000,
        counts_penalty_s2=100000,
        event_length_penalty_k0=1.0,
        tree_structure_prior_k1=0.01,
        use_event_lengths_in_attachment=True,
        seed=seed,
        mixture_size=4,
        num_replicas=5,
        threads_likelihood=4,
        verbose=True,
        neutral_cn=2,
        output_dir=intermediate_dir
    )
    conet.infer_tree(params)
    result = InferenceResult(params.output_dir, cc)
    inferred_cnvs = np.transpose(result.get_inferred_copy_numbers(neutral_cn=int(params.neutral_cn)))
    np.savetxt(out_cnvs, inferred_cnvs, delimiter=",")


    result = InferenceResult(intermediate_dir, cc)

    result.dump_results_to_dir(intermediate_dir, neutral_cn=2)

    # Move results to proper output
    shutil.copy(os.path.join(intermediate_dir, "inferred_attachment"), out_attachments)
    shutil.copy(os.path.join(intermediate_dir, "inferred_tree"), out_tree)

    newick = tree_to_newick(result.inferred_tree)
    newick = newick.replace("\"", "").replace("}", ">").replace("{", "<").replace(":", "")
    with open(os.path.join(intermediate_dir, "tree_newick.nwk"), "w") as f:
        f.write(
            newick
        )
