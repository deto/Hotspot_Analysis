import os
import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
import json

import scanpy.api as sc
from sklearn.manifold import TSNE


latent_file = snakemake.input["latent"]
#latent_file = "pca/pca.txt"
tsne_file = snakemake.input["tsne"]
#tsne_file = "tsne/tsne.txt"

cluster_colors_file = snakemake.output['cluster_colors']
cluster_plot_file = snakemake.output['cluster_plot']
clusters_file = snakemake.output['clusters']

out_dir = os.path.dirname(clusters_file)
if not os.path.exists(out_dir):
    os.mkdir(out_dir)


try:
    n_neighbors = snakemake.params['n_neighbors']
except AttributeError:
    n_neighbors = 30

try:
    resolution = snakemake.params['resolution']
except AttributeError:
    resolution = 1

pc_data = pd.read_table(latent_file, index_col=0)

# %% Cluster Data

adata = sc.AnnData(pc_data.values,
                   obs={'smp_names': pc_data.index.tolist()},
                   var={'var_names': pc_data.columns.tolist()})

print('Running Louvain clustering with n_neighbors={} and resolution={}'.format(
    n_neighbors, resolution))

sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=0)

sc.tl.louvain(adata, resolution=resolution)

clusters = adata.obs['louvain']
clusters = ['c' + x for x in clusters]
sample_names = adata.obs_names
clusters = pd.Series(clusters, index=sample_names, name="Cluster")

if len(clusters.unique()) > 10:
    cmap = plt.get_cmap('tab20')
else:
    cmap = plt.get_cmap('tab10')

cluster_colors = {}
for i, cluster in enumerate(clusters.unique()):
    cluster_colors[cluster] = cmap.colors[i % len(cmap.colors)]


# %% Create a tSNE plot of cluster labels
plt.style.use('default')
tsne = pd.read_table(tsne_file, index_col=0)
tsne.columns = ['tsne1', 'tsne2']
tsne = tsne.join(clusters)

ms = 3
if tsne.shape[0] > 10000:
    ms = 2
if tsne.shape[0] > 30000:
    ms = 1

ff = plt.figure(figsize=(10, 5))
for clust in tsne["Cluster"].unique():
    samples = tsne.loc[tsne["Cluster"] == clust]
    color = cluster_colors[clust]
    plt.plot(samples['tsne1'], samples['tsne2'],
             'o', ms=ms, label=clust, color=color)

plt.subplots_adjust(right=.6)
plt.legend(loc='center left', bbox_to_anchor=(1, .5), markerscale=3/ms)
#plt.show()
plt.savefig(cluster_plot_file)

# Save clusters to disk
clusters.to_frame().to_csv(clusters_file, sep='\t')

# Save cluster colors to disk
json.dump(cluster_colors, open(cluster_colors_file, "w"))
