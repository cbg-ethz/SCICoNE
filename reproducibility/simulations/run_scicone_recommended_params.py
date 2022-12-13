from scicone import SCICoNE, Tree
import numpy as np

# As in config.json
parameters = {
  "bp_detection":
  {
    "window_size":20,
    "threshold":3,
    "bp_limit":200,
  },
  "inference":
  {
    "cluster_trees":
    {
      "n_reps":10,
      "n_iters":40000,
      "move_probs":[0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.01, 0.1, 0.01, 1.0, 0.01],
      "n_tries":3
    },
    "full_trees":
    {
      "n_reps":10,
      "n_iters":200000,
      "move_probs":[0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.0, 1.0, 0.01, 0.1, 0.01, 1.0, 0.01],
    },
  }
}

binaries_path = '<your-path-to/SCICoNE/build>' # or wherever else you installed SCICoNE
output_path = './' # or wherever else you wish to put temporary files

if __name__ == '__main__':

    # Load data
    counts = np.loadtxt('<your-cells-bins-read-counts-file>', delimiter=',')

    # Create SCICoNE object
    sci = SCICoNE(binaries_path, output_path, persistence=False)

    # Detect breakpoints
    bps = sci.detect_breakpoints(counts, window_size=parameters['bp_detection']['window_size'], threshold=parameters['bp_detection']['threshold'], bp_limit=parameters['bp_detection']['bp_limit'])

    # Run cluster trees
    alpha = 1./bps['segmented_regions'].shape[0]
    gamma = 1./np.mean(np.sum(counts, axis=1))/counts.shape[0] # 1/coverage

    sci.learn_tree(counts, bps['segmented_region_sizes'], verbosity=2, n_reps=parameters['inference']['cluster_trees']['n_reps'], cluster=True, full=False, cluster_tree_n_iters=parameters['inference']['cluster_trees']['n_iters'], max_tries=parameters['inference']['cluster_trees']['n_tries'], robustness_thr=0.5, alpha=alpha, max_scoring=True, gamma=gamma)

    # Save cluster tree CNVs: cell-by-bin CNV matrix
    np.savetxt('scicone_cluster_tree_cnvs.csv', sci.best_cluster_tree.outputs['inferred_cnvs'], delimiter=',')

    # Run full trees starting from cluster tree
    sci.learn_tree(counts, bps['segmented_region_sizes'], verbosity=2, n_reps=10, cluster=False, full=True,     full_tree_n_iters=parameters['inference']['full_trees']['n_iters'], max_tries=3, robustness_thr=0.5, alpha=alpha, move_probs=parameters['inference']['full_trees']['move_probs'], max_scoring=True)

    # Save full tree CNVs: cell-by-bin CNV matrix
    np.savetxt('scicone_full_tree_cnvs.csv', sci.best_full_tree.outputs['inferred_cnvs'], delimiter=',')
