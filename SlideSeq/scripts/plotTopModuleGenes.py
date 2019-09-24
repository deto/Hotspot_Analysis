from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import loompy

plt.rcParams["svg.fonttype"] = "none"

results_hs = pd.read_table("../Puck_180819_12/hotspot/hotspot.txt", index_col=0)

positions = pd.read_table(
    "../Puck_180819_12/positions/positions.txt", index_col=0
)

# %% Load clusters

clusters = pd.read_table(
    "../Puck_180819_12/hotspot/modules.txt", index_col=0
).Cluster

cluster_map = {}
for gene, cl in clusters.iteritems():
    if cl not in cluster_map:
        cluster_map[cl] = set()

    cluster_map[cl].add(gene)

valid_clusters = [x for x in cluster_map if x != -1]


# %% Plot modules


def get_top_gene(genes):
    tg = results_hs.loc[genes].sort_values("C", ascending=False).index[0]
    return tg


top_genes = {k: get_top_gene(v) for k, v in cluster_map.items()}


ds = loompy.connect("../../data/SlideSeq/Puck_180819_12/data.loom", "r")
genes = ds.row_attrs["EnsID"][:]

from matplotlib.colors import LinearSegmentedColormap

cmap = LinearSegmentedColormap.from_list("grays", ["#E0E0E0", "#000000"])

for module, gene in top_genes.items():
    idx = np.where(genes == gene)[0][0]
    expression = np.log2(ds.layers["scaled"][idx, :] + 1)

    fig = plt.figure(figsize=(4.7, 4))
    plt.scatter(
        positions.Comp1,
        positions.Comp2,
        s=1,
        c=expression,
        vmax=np.percentile(expression, 99),
        cmap=cmap,
        edgecolors=None,
        linewidths=0,
        rasterized=True,
    )
    plt.title("Module {}: {}".format(module, gene))
    plt.xticks([])
    plt.yticks([])
    plt.colorbar()
    fig.patch.set_visible(False)
    plt.savefig("TopGene_Module_{}.svg".format(module), dpi=300)
