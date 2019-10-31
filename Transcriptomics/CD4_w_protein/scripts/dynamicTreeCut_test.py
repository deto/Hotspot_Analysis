from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
import numpy as np
import pandas as pd
import __main__ as main
if hasattr(main, '__file__'):
    import matplotlib
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import seaborn as sns
from tqdm import tqdm

plt.rcParams["svg.fonttype"] = "none"

# %%

results_file = "hotspot/hotspot_pairs_z.txt.gz"
results_hs_file = "hotspot/hotspot_hvg.txt"

MIN_CLUSTER_GENES = 50
MIN_CLUSTER_Z = 7

results = pd.read_table(results_file, index_col=0)
results_hs = pd.read_table(results_hs_file, index_col=0)

ens_map = {x: y for x, y in zip(results_hs.index, results_hs.Symbol)}

# %% Compute Linkage and Ordering

from scipy.spatial.distance import squareform
dd = results.copy().values
np.fill_diagonal(dd, 0)
condensed = squareform(dd)*-1
offset = condensed.min() * -1
condensed += offset
Z = linkage(condensed, method='average')
# Z = linkage(condensed, method='complete')

# %% Can we implement this from the linkage matrix?


# Z is a Mx4 matrix
# Each row has four elements:
# 1) Node A
# 2) Node B
# 3) Join height
# 4) New number of elements
# Each row defines a new node in the tree


# now, traverse the tree and propagate labels downward if they are not -1


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

# Try something a little different
# Instead, leaves start with -1 as an assignment
# they only get something better on merge
clust_i = 0

clusters = np.ones(Z.shape[0])*-1
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
        clust_a = clusters[ca-N]

    if cb - N < 0:  # leaf node
        n_members_b = 1
        clust_b = -1
    else:
        n_members_b = Z[cb-N, 3]
        clust_b = clusters[cb-N]

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
        new_clust_assign = -1 # Still too small

    clusters[i] = new_clust_assign

labels = clusters
out_clusters = np.ones(N)*-2
prop_label(Z, Z.shape[0]-1, labels[-1], labels, out_clusters)

# remap out_clusters
clust_map = {x: i-1 for i, x in enumerate(np.sort(np.unique(out_clusters)))}
out_clusters = [clust_map[x] for x in out_clusters]
out_clusters = pd.Series(out_clusters, index = results.index)


# %% How does this look?
# Looks much better!

colors = list(plt.get_cmap("tab10").colors)
cm = ScalarMappable(norm=Normalize(0, 0.05, clip=True), cmap="viridis")
# row_colors1 = pd.Series(
#     [colors[i % 10] if i != -1 else "#ffffff" for i in clusters],
#     index=results.index,
# )
row_colors2 = pd.Series(
    [colors[i % 10] if i != -1 else "#ffffff" for i in out_clusters],
    index=results.index,
)

row_colors = pd.DataFrame({
    #"Cluster": row_colors1,
    "Cluster-X": row_colors2,
})

cm = sns.clustermap(
    results,
    row_linkage=Z,
    col_linkage=Z,
    vmin=-25,
    vmax=25,
    cmap="RdBu_r",
    xticklabels=False,
    yticklabels=False,
    row_colors=row_colors,
    rasterized=True,
)

fig = plt.gcf()
fig.patch.set_visible(False)
plt.sca(cm.ax_heatmap)
plt.ylabel("")
plt.xlabel("")
plt.show()

# %% Next thing to do would be to try to assign clusters to their closest node
# on post-process

#  out_clusters2 = out_clusters.copy()
#  
#  
#  r_na = results.copy()
#  np.fill_diagonal(r_na.values, float('nan')) # trick to ignore self-distances
#  
#  # Mean distance to each cluster
#  c_dists = r_na.groupby(out_clusters).mean().drop(-1, axis=0).T
#  
#  to_assign = results.index[np.array(out_clusters) == -1]
#  
#  c_dists = c_dists.loc[to_assign]
#  c_dists = c_dists.loc[c_dists.max(axis=1) > 3]
#  new_assignments = c_dists.idxmax(axis=1)
#  out_clusters2[new_assignments.index] = new_assignments
#  
#  # %% How does this look?
#  
#  colors = list(plt.get_cmap("tab10").colors)
#  cm = ScalarMappable(norm=Normalize(0, 0.05, clip=True), cmap="viridis")
#  row_colors1 = pd.Series(
#      [colors[i % 10] if i != -1 else "#ffffff" for i in out_clusters],
#      index=results.index,
#  )
#  row_colors2 = pd.Series(
#      [colors[i % 10] if i != -1 else "#ffffff" for i in out_clusters2],
#      index=results.index,
#  )
#  
#  row_colors = pd.DataFrame({
#      "Cluster-X": row_colors1,
#      "Cluster-X2": row_colors2,
#  })
#  
#  cm = sns.clustermap(
#      results,
#      row_linkage=Z,
#      col_linkage=Z,
#      vmin=-25,
#      vmax=25,
#      cmap="RdBu_r",
#      xticklabels=False,
#      yticklabels=False,
#      row_colors=row_colors,
#      rasterized=True,
#  )
#  
#  fig = plt.gcf()
#  fig.patch.set_visible(False)
#  plt.sca(cm.ax_heatmap)
#  plt.ylabel("")
#  plt.xlabel("")
#  plt.show()

# This provides a solution, but the resulting visualizations
# are NOT consistent with the dendrogram

# What if we instead try to solve it by making assignments consistent with dendrogram?
# This just depends on how they draw the dendrogram...
# What if just the larger one propagates upward?

clust_i = 0

clusters = np.ones(Z.shape[0])*-1
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
        clust_a = clusters[ca-N]

    if cb - N < 0:  # leaf node
        n_members_b = 1
        clust_b = -1
    else:
        n_members_b = Z[cb-N, 3]
        clust_b = clusters[cb-N]

    if Z[i, 2] > offset - 3:
        new_clust_assign = -1
    elif (n_members_a >= MIN_THRESHOLD and n_members_b >= MIN_THRESHOLD):
        # don't join them
        if n_members_a >= n_members_b:
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
        new_clust_assign = -1 # Still too small

    clusters[i] = new_clust_assign

# Now to read off the actual assignment, just propagate labels downward always
# overwriting parent labels with child labels unless label is -1

labels = clusters
out_clusters2 = np.ones(N)*-2

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

prop_label2(Z, Z.shape[0]-1, labels[-1], labels, out_clusters2)

# remap out_clusters
clust_map = {x: i-1 for i, x in enumerate(np.sort(np.unique(out_clusters2)))}
out_clusters2 = [clust_map[x] for x in out_clusters2]
out_clusters2 = pd.Series(out_clusters2, index=results.index)

# %% Plot the difference

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
    row_linkage=Z,
    col_linkage=Z,
    vmin=-25,
    vmax=25,
    cmap="RdBu_r",
    xticklabels=False,
    yticklabels=False,
    row_colors=row_colors,
    rasterized=True,
)

fig = plt.gcf()
fig.patch.set_visible(False)
plt.sca(cm.ax_heatmap)
plt.ylabel("")
plt.xlabel("")
plt.show()

# This still looks bad
# The issue is that the dendrogram is kind of arbitrary?
# When seaborn calls it, they don't order the nodes in any particular way?
# How then, does dynamicTreeCut do it and make it look nice???

import scipy.cluster.hierarchy as hierarchy

plt.figure()
hierarchy.dendrogram(Z, labels=out_clusters, leaf_font_size=8)
plt.show()

# What if instead we just preserve the identity of the cluster that
# has the smallest mean distance - e.g., join new nodes to the cluster
# that is spread out more?

# Need to calculate mean join height under every tree node
mean_dists = np.zeros(Z.shape[0])

def calc_mean_dists(Z, node_index, out_mean_dists):

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

calc_mean_dists(Z, Z.shape[0]-1, mean_dists)

clust_i = 0

clusters = np.ones(Z.shape[0])*-1
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
        clust_a = clusters[ca-N]

    if cb - N < 0:  # leaf node
        n_members_b = 1
        clust_b = -1
    else:
        n_members_b = Z[cb-N, 3]
        clust_b = clusters[cb-N]

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
        new_clust_assign = -1 # Still too small

    clusters[i] = new_clust_assign

labels = clusters
out_clusters2 = np.ones(N)*-2

prop_label2(Z, Z.shape[0]-1, labels[-1], labels, out_clusters2)

# remap out_clusters
clust_map = {x: i-1 for i, x in enumerate(np.sort(np.unique(out_clusters2)))}
out_clusters2 = [clust_map[x] for x in out_clusters2]
out_clusters2 = pd.Series(out_clusters2, index=results.index)

# %% Plot the difference

colors = list(plt.get_cmap("tab20").colors)
cm = ScalarMappable(norm=Normalize(0, 0.05, clip=True), cmap="viridis")
row_colors1 = pd.Series(
    [colors[i % 10] if i != -1 else "#ffffff" for i in out_clusters],
    index=results.index,
)
row_colors2 = pd.Series(
    [colors[i % 20] if i != -1 else "#ffffff" for i in out_clusters2],
    index=results.index,
)

row_colors = pd.DataFrame({
    "Cluster-X": row_colors1,
    "Cluster-X2": row_colors2,
})

cm = sns.clustermap(
    results,
    row_linkage=Z,
    col_linkage=Z,
    vmin=-25,
    vmax=25,
    cmap="RdBu_r",
    xticklabels=False,
    yticklabels=False,
    row_colors=row_colors,
    rasterized=True,
)

fig = plt.gcf()
fig.patch.set_visible(False)
plt.sca(cm.ax_heatmap)
plt.ylabel("")
plt.xlabel("")
plt.show()

# This is looking pretty good!

# The last thing we are missing is a re-ordering of the linkage
# object (swap leaf order) so that more dense objects come first


def sort_linkage(Z, node_index, node_values):

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

Z_sorted = Z.copy()
sort_linkage(Z_sorted, Z.shape[0]-1, mean_dists)

# %% Plot sorted

colors = list(plt.get_cmap("tab20").colors)
cm = ScalarMappable(norm=Normalize(0, 0.05, clip=True), cmap="viridis")
row_colors1 = pd.Series(
    [colors[i % 10] if i != -1 else "#ffffff" for i in out_clusters],
    index=results.index,
)
row_colors2 = pd.Series(
    [colors[i % 20] if i != -1 else "#ffffff" for i in out_clusters2],
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

# %% Ok - now make things a bit prettier


