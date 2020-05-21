import glob
from tqdm import tqdm

all_folders = glob.glob("/cluster/work/bewi/members/pedrof/sc-dna/sims2020_new/results/*")

methods = ['hmm_copy', 'hclust', 'phenograph', 'cluster_tree', 'full_tree']

n_reps = 40
n_nodes = 10
n_bins = 10000
n_regions = [n_nodes, 2*n_nodes, 4*n_nodes]
n_reads = [2*n_bins, 4*n_bins, 8*n_bins]
sim_folder = [folder for folder in all_folders if "simulation" in folder]

row_list = []
for n_regions in n_regions:
    for n_reads in n_reads:
        for n_rep in range(n_reps):
            gt_path = sim_folder + f"/{n_nodes}nodes_{n_regions}regions_{n_reads}reads/" + str(rep_id) + "_ground_truth.txt"
            gt_mat = np.loadtxt(gt_path, delimiter=',')

            for method in methods:
                method_folder = [folder for folder in all_folders if method in folder][0]
                if 'tree' in method:
                    inferred_cnvs_path = method_folder + f"/{n_nodes}nodes_{n_regions}regions_{n_reads}reads/" + str(rep_id) + '_' + method + "_cnvs.csv"
                else:
                    inferred_cnvs_path = method_folder + f"/{n_nodes}nodes_{n_regions}regions_{n_reads}reads/" + str(rep_id) + '_' + method + "_inferred.txt"
                inferred_cnvs = np.loadtxt(inferred_cnvs_path, delimiter=',')

                delta = ((gt_mat - inferred_cnvs) ** 2).mean()
                tup = (rep_id, n_nodes, n_regions, n_reads, method, delta)
                row_list.append(tup)

df = pd.DataFrame(data=row_list, columns=['rep_id', 'n_nodes', 'n_regions', 'n_reads', 'method', 'delta'])
df.to_csv(f"/cluster/work/bewi/members/pedrof/sc-dna/sims2020_new/{n_nodes}nodes_deltas.csv")
