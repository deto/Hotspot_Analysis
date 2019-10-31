"""
This file was where these methods were initially developed
However, now they are implemented in the 'modules.py' module
of the Hotspot package in python
"""

import numpy as np
import pandas as pd


def sort_linkage(Z, node_index, node_values):
    """
    Sorts linkage by 'node_values' in place
    """

    N = Z.shape[0] + 1  # number of leaves

    if node_index < 0:
        return

    left_child = int(Z[node_index, 0] - N)
    right_child = int(Z[node_index, 1] - N)

    swap = False

    if left_child < 0 and right_child < 0:
        swap = False
    elif left_child < 0 and right_child >= 0:
        swap = True
    elif left_child >= 0 and right_child < 0:
        swap = False
    else:
        if node_values[left_child] > node_values[right_child]:
            swap = True
        else:
            swap = False

    if swap:
        Z[node_index, 0] = right_child + N
        Z[node_index, 1] = left_child + N

    sort_linkage(Z, left_child, node_values)
    sort_linkage(Z, right_child, node_values)


def calc_mean_dists(Z, node_index, out_mean_dists):
    """
    Calculates the mean density of joins
    for sub-trees underneath each node
    """

    N = Z.shape[0] + 1  # number of leaves

    left_child = int(Z[node_index, 0] - N)
    right_child = int(Z[node_index, 1] - N)

    if left_child < 0:
        left_average = 0
        left_merges = 0
    else:
        left_average, left_merges = calc_mean_dists(
            Z, left_child, out_mean_dists
        )

    if right_child < 0:
        right_average = 0
        right_merges = 0
    else:
        right_average, right_merges = calc_mean_dists(
            Z, right_child, out_mean_dists
        )

    this_height = Z[node_index, 2]
    this_merges = left_merges + right_merges + 1
    this_average = (
        left_average * left_merges + right_average * right_merges + this_height
    ) / this_merges

    out_mean_dists[node_index] = this_average

    return this_average, this_merges


def prop_label(Z, node_index, label, labels, out_clusters):
    """
    Propagates node labels downward if they are not -1
    Used to find the correct cluster label at the leaves
    """

    N = Z.shape[0] + 1  # number of leaves

    if label == -1:
        label = labels[node_index]

    left_child = int(Z[node_index, 0] - N)
    right_child = int(Z[node_index, 1] - N)

    if left_child < 0:
        out_clusters[left_child + N] = label
    else:
        prop_label(Z, left_child, label, labels, out_clusters)

    if right_child < 0:
        out_clusters[right_child + N] = label
    else:
        prop_label(Z, right_child, label, labels, out_clusters)


def prop_label2(Z, node_index, label, labels, out_clusters):

    N = Z.shape[0] + 1  # number of leaves

    parent_label = label
    this_label = labels[node_index]

    if this_label == -1:
        new_label = parent_label
    else:
        new_label = this_label

    left_child = int(Z[node_index, 0] - N)
    right_child = int(Z[node_index, 1] - N)

    if left_child < 0:
        out_clusters[left_child + N] = new_label
    else:
        prop_label2(Z, left_child, new_label, labels, out_clusters)

    if right_child < 0:
        out_clusters[right_child + N] = new_label
    else:
        prop_label2(Z, right_child, new_label, labels, out_clusters)


def assign_modules(Z, leaf_labels, offset, MIN_THRESHOLD=10):
    clust_i = 0

    labels = np.ones(Z.shape[0])*-1
    N = Z.shape[0]+1
    MIN_THRESHOLD = 10

    mean_dists = np.zeros(Z.shape[0])
    calc_mean_dists(Z, Z.shape[0]-1, mean_dists)

    for i in range(Z.shape[0]):

        ca = int(Z[i, 0])
        cb = int(Z[i, 1])

        if ca - N < 0:  # leaf node
            n_members_a = 1
            clust_a = -1
        else:
            n_members_a = Z[ca-N, 3]
            clust_a = labels[ca-N]

        if cb - N < 0:  # leaf node
            n_members_b = 1
            clust_b = -1
        else:
            n_members_b = Z[cb-N, 3]
            clust_b = labels[cb-N]

        if Z[i, 2] > offset - 3:
            new_clust_assign = -1
        elif (n_members_a >= MIN_THRESHOLD and n_members_b >= MIN_THRESHOLD):
            # don't join them
            # assign the one with the larger mean distance
            dist_a = mean_dists[ca-N]
            dist_b = mean_dists[cb-N]
            if dist_a >= dist_b:
                new_clust_assign = clust_a
            else:
                new_clust_assign = clust_b
        elif n_members_a >= MIN_THRESHOLD:
            new_clust_assign = clust_a
        elif n_members_b >= MIN_THRESHOLD:
            new_clust_assign = clust_b
        elif (n_members_b + n_members_a) >= MIN_THRESHOLD:
            # A new cluster is born!
            new_clust_assign = clust_i
            clust_i += 1
        else:
            new_clust_assign = -1  # Still too small

        labels[i] = new_clust_assign

    out_clusters = np.ones(N)*-2
    prop_label2(Z, Z.shape[0]-1, labels[-1], labels, out_clusters)

    # remap out_clusters
    clust_map = {
        x: i-1 for i, x in enumerate(np.sort(np.unique(out_clusters)))
    }
    out_clusters = [clust_map[x] for x in out_clusters]
    out_clusters = pd.Series(out_clusters, index=leaf_labels)

    return out_clusters


def assign_modules_core(Z, leaf_labels, offset, MIN_THRESHOLD=10):
    clust_i = 0

    labels = np.ones(Z.shape[0])*-1
    N = Z.shape[0]+1
    MIN_THRESHOLD = 10

    for i in range(Z.shape[0]):

        ca = int(Z[i, 0])
        cb = int(Z[i, 1])

        if ca - N < 0:  # leaf node
            n_members_a = 1
            clust_a = -1
        else:
            n_members_a = Z[ca-N, 3]
            clust_a = labels[ca-N]

        if cb - N < 0:  # leaf node
            n_members_b = 1
            clust_b = -1
        else:
            n_members_b = Z[cb-N, 3]
            clust_b = labels[cb-N]

        if (n_members_a >= MIN_THRESHOLD and n_members_b >= MIN_THRESHOLD):
            # don't join them
            new_clust_assign = -1
        elif Z[i, 2] > offset - 3:
            new_clust_assign = -1
        elif n_members_a >= MIN_THRESHOLD:
            new_clust_assign = clust_a
        elif n_members_b >= MIN_THRESHOLD:
            new_clust_assign = clust_b
        elif (n_members_b + n_members_a) >= MIN_THRESHOLD:
            # A new cluster is born!
            new_clust_assign = clust_i
            clust_i += 1
        else:
            new_clust_assign = -1  # Still too small

        labels[i] = new_clust_assign

    out_clusters = np.ones(N)*-2
    prop_label(Z, Z.shape[0]-1, labels[-1], labels, out_clusters)

    # remap out_clusters
    clust_map = {
        x: i-1 for i, x in enumerate(np.sort(np.unique(out_clusters)))
    }
    out_clusters = [clust_map[x] for x in out_clusters]
    out_clusters = pd.Series(out_clusters, index=leaf_labels)

    return out_clusters


# %% Test it out - load data

results_file = "hotspot/hotspot_pairs_z.txt.gz"
results_hs_file = "hotspot/hotspot_hvg.txt"

results = pd.read_table(results_file, index_col=0)
results_hs = pd.read_table(results_hs_file, index_col=0)

ens_map = {x: y for x, y in zip(results_hs.index, results_hs.Symbol)}

# %% Compute Linkage and Ordering

from scipy.spatial.distance import squareform
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import seaborn as sns

dd = results.copy().values
np.fill_diagonal(dd, 0)
condensed = squareform(dd)*-1
offset = condensed.min() * -1
condensed += offset
Z = linkage(condensed, method='average')

# %%

out_clusters = assign_modules_core(
    Z, offset=offset, MIN_THRESHOLD=10, leaf_labels=results.index)

out_clusters2 = assign_modules(
    Z, offset=offset, MIN_THRESHOLD=10, leaf_labels=results.index)

mean_dists = np.zeros(Z.shape[0])
calc_mean_dists(Z, Z.shape[0]-1, mean_dists)
Z_sorted = Z.copy()
sort_linkage(Z_sorted, Z.shape[0]-1, mean_dists)

# %% Plot it

colors = list(plt.get_cmap("tab10").colors)
cm = ScalarMappable(norm=Normalize(0, 0.05, clip=True), cmap="viridis")
row_colors1 = pd.Series(
    [colors[i % 10] if i != -1 else "#ffffff" for i in out_clusters],
    index=results.index,
)
row_colors2 = pd.Series(
    [colors[i % 10] if i != -1 else "#ffffff" for i in out_clusters2],
    index=results.index,
)

row_colors = pd.DataFrame({
    "Cluster-X": row_colors1,
    "Cluster-X2": row_colors2,
})

cm = sns.clustermap(
    results,
    row_linkage=Z_sorted,
    col_linkage=Z_sorted,
    vmin=-25,
    vmax=25,
    cmap="RdBu_r",
    xticklabels=False,
    #yticklabels=False,
    yticklabels=[ens_map[x] for x in results.index],
    row_colors=row_colors,
    rasterized=True,
)

fig = plt.gcf()
fig.patch.set_visible(False)
plt.sca(cm.ax_heatmap)
plt.ylabel("")
plt.xlabel("")
plt.show()

# %% Plot the modules on our UMAP

proj_file = "umap/umap_hvg.txt"
loom_file = "../data/10x_PBMC_w_proteins/cd4/data.loom"
latent_file = "scvi/hvg/latent.txt.gz"
model = "danb"
n_neighbors = 30


import loompy
import hotspot
import hotspot.modules


proj = pd.read_table(proj_file, index_col=0)

with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

latent = pd.read_table(latent_file, index_col=0)

# Have to do this because data_slideseq makes it a numpy array
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)
num_umi = pd.Series(num_umi, index=barcodes)

# Align to latent space
counts = counts.loc[:, latent.index]
num_umi = num_umi[latent.index]
proj = proj.loc[latent.index]

# need counts, latent, and num_umi

hs = hotspot.Hotspot(counts, latent, num_umi)
hs.create_knn_graph(
    weighted_graph=False, n_neighbors=n_neighbors, neighborhood_factor=3
)

# %%

module = 5
module_genes = out_clusters.index[out_clusters == module]

scores = hotspot.modules.compute_scores(
    counts.loc[module_genes].values, model, num_umi.values,
    hs.neighbors.values, hs.weights.values
)


# %% Plot scores

vmin = np.percentile(scores, 5)
vmax = np.percentile(scores, 95)

plt.figure()
plt.scatter(
    x=proj.iloc[:, 0],
    y=proj.iloc[:, 1],
    c=scores,
    vmin=vmin, vmax=vmax,
    s=2
)
plt.show()

# %% Plot scores for all modules

clusters = out_clusters2
fig, axs = plt.subplots(3, 3, figsize=(9, 9))

modules_to_plot = sorted([x for x in clusters.unique() if x != -1])

# Get the scores
module_scores = {}
for module in modules_to_plot:
    module_genes = clusters.index[clusters == module]

    scores = hotspot.modules.compute_scores(
        counts.loc[module_genes].values, model, num_umi.values,
        hs.neighbors.values, hs.weights.values
    )

    module_scores[module] = scores


for ax, module in zip(axs.ravel(), modules_to_plot):
    plt.sca(ax)

    scores = module_scores[module]

    vmin = np.percentile(scores, 5)
    vmax = np.percentile(scores, 95)

    plt.scatter(
        x=proj.iloc[:, 0],
        y=proj.iloc[:, 1],
        c=scores,
        vmin=vmin, vmax=vmax,
        s=2
    )
    plt.title(module)

plt.show()

# %%

# These modules are opposites

plt.figure()
plt.plot(module_scores[1], module_scores[2], 'o', ms=2)
plt.show()

# %% These are the same?

plt.figure()
plt.plot(module_scores[5], module_scores[8], 'o', ms=2)
plt.show()

# %% These are the same?

plt.figure()
plt.plot(module_scores[0], module_scores[7], 'o', ms=2)
plt.show()
