import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# %% Load module genes

modules_hs_df = pd.read_table(
    "../Puck_180819_12/hotspot/modules.txt", index_col=0
)

modules_hs = {}

for g, c in modules_hs_df.Cluster.iteritems():
    if c == -1: continue
    if c not in modules_hs:
        modules_hs[c] = set()

    modules_hs[c].add(g)

modules_sd_df = pd.read_table(
    "../Puck_180819_12/spatialDE/histology_windowed_10.txt", index_col=0
)


modules_sd = {}

for row in modules_sd_df.itertuples():
    if row.pattern not in modules_sd:
        modules_sd[row.pattern] = set()

    modules_sd[row.pattern].add(row.g)

# %% Compute overlap

def jaccard(setA, setB):
    setA = {x.lower() for x in setA}
    setB = {x.lower() for x in setB}

    intersection = len(setA & setB)
    union = len(setA | setB)

    return intersection/union


hs_keys = sorted(modules_hs.keys())
sd_keys = sorted(modules_sd.keys())

overlaps = pd.DataFrame(
    0, index=hs_keys, columns=sd_keys
)

for hs, hs_genes in modules_hs.items():
    for sd, sd_genes in modules_sd.items():

        overlaps.loc[hs, sd] = jaccard(hs_genes, sd_genes)

plt.figure()
sns.heatmap(overlaps)
plt.show()

# %% Plot modules for spatial DE

patterns_sd_df = pd.read_table(
    "../Puck_180819_12/spatialDE/patterns_windowed_10.txt", index_col=0
)

positions = pd.read_table(
    "../Puck_180819_12/positions/positions.txt", index_col=0
)

patterns = patterns_sd_df.columns

data = positions.join(patterns_sd_df).fillna(0)

fig, axs = plt.subplots(2, 5, figsize=(12, 5), sharex=True, sharey=True)

for pat, ax in zip(patterns, axs.ravel()):
    plt.sca(ax)
    plt.scatter(data.Comp1, data.Comp2, s=2, c=data[pat], vmin=-.5, vmax=.5)
    plt.title(pat)
    plt.xticks([])
    plt.yticks([])

plt.show()

# %% How do these modules break down between cell types?

# %% Load markers

markers_df = pd.read_table("../DropViz/edger_markers_1vAll.txt")

N = 200
markers = {}

for cluster, group in markers_df.groupby('Cluster'):
    group_markers = group.sort_values('FDR').GeneSymbol[0:N]
    markers[cluster] = set(group_markers)

# anno = pd.read_excel("../DropViz/annotation.BrainCellAtlas_Saunders_version_2018.04.01.xlsx")
# Taken from the above, but re-labeled for primary clusters (not sub-clusters)

marker_name_map = {
    1: 'Granule Cells',
    2: 'Purkinje Neurons',
    3: 'Cerebellum Basket Cells',
    4: 'Other Neurons',
    5: 'Resident Macrophages',
    6: 'Oligo/poly-dendrocytes',
    7: 'Bergmann Glia',
    8: 'Astrocytes',
    9: 'Choroid Plexus',
    10: 'Endothelial Stalk',
    11: 'Endothelial Tip / Mural',
}

markers = {marker_name_map[k]: v for k, v in markers.items()}

# Make it so we only use markers that are in only one set
marker_counts = []
for ct, genes in markers.items():
    marker_counts.extend(genes)

from collections import Counter
marker_counts = Counter(marker_counts)
marker_counts_singles = {k for k, v in marker_counts.items() if v == 1}
markers = {
    ct: genes & marker_counts_singles for
    ct, genes in markers.items()
}

# %% Plot overlaps

# What is the expected overlap?

background = set(modules_sd_df.g)
valid_clusters = modules_sd.keys()
cluster_map = modules_sd

overlaps = pd.DataFrame(
    0.0, index=valid_clusters,
    columns=markers.keys()
)
fcs = pd.DataFrame(
    0.0, index=valid_clusters,
    columns=markers.keys()
)

for celltype in markers:
    marker_genes = markers[celltype]
    marker_genes = marker_genes & background

    for cluster in valid_clusters:
        cluster_genes = cluster_map[cluster]
        cluster_genes = cluster_genes & background

        expected = len(cluster_genes) * len(marker_genes) / len(background)
        actual = len(cluster_genes & marker_genes)

        fc = (actual + 1) / (expected + 1)

        overlaps.loc[cluster, celltype] = actual / \
            len(cluster_genes | marker_genes)
        fcs.loc[cluster, celltype] = fc

cm = sns.clustermap(
    np.log2(fcs).T,
    metric='correlation',
    vmin=-4,
    vmax=4,
    cmap="RdBu_r",
    col_cluster=False,
    figsize=(7, 7),
    rasterized=True
)

fig = plt.gcf()
fig.patch.set_visible(False)
plt.sca(cm.ax_heatmap)
plt.xlabel('Gene Module')
plt.subplots_adjust(bottom=.2, right=.7)
# plt.savefig('Module_Cell_Correspondence_sd.svg', dpi=300)
plt.show()
