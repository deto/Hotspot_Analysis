import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
plt.rcParams['svg.fonttype'] = 'none'

results_hs = pd.read_table(
    "../Puck_180819_12/hotspot/hotspot.txt",
    index_col=0
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

background = set(results_hs.index)

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
plt.savefig('Module_Cell_Correspondence.svg', dpi=300)
#plt.show()

# %% Try with sub-clusters instead

anno = pd.read_excel("../DropViz/annotation.BrainCellAtlas_Saunders_version_2018.04.01.xlsx")
markers_df = pd.read_table("../DropViz/edger_markers_1vAll_subcluster.txt")

N = 200
markers = {}

for cluster, group in markers_df.groupby('Cluster'):
    group_markers = group.sort_values('FDR').GeneSymbol[0:N]
    markers[cluster] = set(group_markers)

# anno = pd.read_excel("../DropViz/annotation.BrainCellAtlas_Saunders_version_2018.04.01.xlsx")
# Taken from the above, but re-labeled for primary clusters (not sub-clusters)

marker_name_map = {
    row.subcluster: '{} {}'.format(row.subcluster, row.common_name)
    for row in anno.loc[anno.tissue == 'CB'].itertuples()
}

markers = {marker_name_map[k]: v for k, v in markers.items()}

# Remove genes in multiple cell types
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

background = set(results_hs.index)

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

sns.clustermap(
    np.log2(fcs).T,
    metric='correlation',
    vmin=-4,
    vmax=4,
    cmap="RdBu_r",
    col_cluster=False,
    figsize=(5, 8)
)

plt.subplots_adjust(bottom=.1, right=.5)
fig = plt.gcf()
fig.patch.set_visible(False)
plt.savefig('Module_SubCell_Correspondence.svg', dpi=300)
#plt.show()
