import argparse
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from sklearn.manifold import TSNE

parser = argparse.ArgumentParser()
parser.add_argument("-p", "--cluster_profile",required=True, help="mean of the cell counts inside the same cluster")
parser.add_argument("-d","--dist", required=True, help="the distance matrix used by phenograph")
parser.add_argument("-c","--communities",required=True,help="the cluster assignments")
parser.add_argument("-o","--output_path",required=False, default="./", help="path to the output")
parser.add_argument("-s", "--sample_name",required=False, default="", help="name of the sample")
parser.add_argument("-g", "--genomic_coordinates", required=True, help="chromosome stop positions")

args = parser.parse_args()

cluster_means = pd.read_csv(args.cluster_profile,sep='\t')
cluster_means.index = cluster_means["cluster_ids"].values
cluster_means = cluster_means.drop(columns="cluster_ids")

# the max val to be used in the plot
max_val = cluster_means.values.max()
max_val

dist = np.loadtxt(args.dist,delimiter=',')

communities = pd.read_csv(args.communities,sep='\t')
communities = communities['cluster'].values

tsne = TSNE(n_components=2, perplexity= 30,metric='precomputed').fit_transform(dist)
df_tsne = pd.DataFrame(tsne)
df_tsne['cluster'] = communities
ax = df_tsne.plot(kind='scatter', x=0, y=1, c='cluster', figsize=(10,8), colorbar=False,colormap='Dark2', grid=True, title='Phenograph Clusters on CNV Data')
fig = ax.get_figure()
fig.savefig(args.output_path + '/' + args.sample_name +  "_tsne_output.png")

chr_stops_df = pd.read_csv(args.genomic_coordinates, sep='\t')
chr_stops_df.columns = ["idx", "pos"]

cmap = matplotlib.cm.get_cmap('Dark2')
# use the formula below to get the distinct colors
# color = cmap(float(i)/N)
for i, cluster_idx in enumerate(cluster_means.index):
    plt.figure(figsize=(20,6))
    ax = plt.plot(cluster_means.iloc[i].values,label="cluster id: " + str(cluster_idx), color=cmap(float(i)/max(communities)))
    plt.axis([None, None, 0, max_val]) # to make the axises same
    plt.legend(loc='upper left')
    plt.xticks([], [])
    for index, row in chr_stops_df.iterrows():
        plt.text(row['idx'],-0.5,"chr " + row['pos'],rotation=90)
    plt.savefig(args.output_path + '/' + args.sample_name + "_cluster_profile_"+str(cluster_idx)+".png")


plt.figure(figsize=(20,6))
for i, cluster_idx in enumerate(cluster_means.index):
    ax = plt.plot(cluster_means.iloc[i].values, label="cluster id: " + str(cluster_idx), color=cmap(float(i)/max(communities)), alpha=0.6)
    plt.axis([None, None, 0, max_val]) # to make the axises same
    plt.legend(loc='upper left')
    plt.xticks([], [])
    for index, row in chr_stops_df.iterrows():
        plt.text(row['idx'],-0.5,"chr " + row['pos'],rotation=90)
plt.savefig(args.output_path + '/' + args.sample_name + "_cluster_profile_overlapping.png")







