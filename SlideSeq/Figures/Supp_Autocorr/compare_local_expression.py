import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'


N = 1000 # How many genes to look at

hs_info = pd.read_table("../../Puck_180819_12/hotspot/hotspot.txt", index_col=0)
hvg_info = pd.read_table("../../Puck_180819_12/genes/hvg_info.txt", index_col=0)
danb_info = pd.read_table("../../Puck_180819_12/genes/danb_info.txt", index_col=0)
pca_info = pd.read_table("../../Puck_180819_12/genes/pca_info.txt", index_col=0)

hvg_info = hvg_info.loc[hs_info.index]
danb_info = danb_info.loc[hs_info.index]
pca_info = pca_info.loc[hs_info.index]

geary = pd.read_table("../../Puck_180819_12/autocorrelation/autocorrelation.txt", index_col=0)

selected_genes = {}

selected_genes['NBDisp'] = danb_info.sort_values('q.value').index[0:N]
selected_genes['HVG'] = hvg_info.sort_values('gene.dispersion.scaled', ascending=False).index[0:N]
selected_genes['PCA'] = pca_info.sort_values('Score', ascending=False).index[0:N]
selected_genes['Hotspot'] = hs_info.sort_values('Z', ascending=False).index[0:N]

# %%

from scipy.stats import gaussian_kde

colors = sns.color_palette('deep')
colormap = {
    'Hotspot': colors[0],
    'HVG': colors[1],
    'NBDisp': colors[2],
    'PCA': colors[3],
}

domain = np.linspace(-.01, .02, 100)

plt.figure(figsize=(6, 3))

for method in selected_genes:
    genes = selected_genes[method]

    geary_vals = geary.loc[genes]['MoranI'].values

    # bin_counts, bin_edges = np.histogram(geary_vals, bins=np.linspace(-.1, .5, 100), density=True)
    # bin_centers = bin_edges[1:]/2 + bin_edges[0:-1]/2
    # plt.bar(bin_centers, bin_counts, label=method, color=colormap[method], alpha=0.5)

    kernel = gaussian_kde(geary_vals, bw_method=.005)
    y_vals = kernel(domain)
    plt.plot(domain, y_vals, label=method, color=colormap[method])

plt.legend()
plt.grid(color='#dddddd', linestyle=(0, (5, 5)))
plt.gca().set_axisbelow(True)
plt.xlabel('Moran-I')
plt.ylabel('KDE Estimate BW=0.005')
plt.title('Distribution for Top 1000 Genes Selected by Each Method')
plt.subplots_adjust(bottom=0.2)
plt.savefig('moran_i_kde.svg', dpi=300)
# plt.show()
plt.close()

# %%

# Plot a clustermap of genes from one method or another

import loompy
import hotspot
from numba import jit
from tqdm import tqdm

loom_file = "../../../data/SlideSeq/Puck_180819_12/data.loom"

with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]


gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)
num_umi = pd.Series(num_umi, index=barcodes)


# Normalize and scale counts
norm_counts = np.log2(counts.values / num_umi.values.reshape((1, -1)) * 50 + 1)
norm_counts = pd.DataFrame(norm_counts, index=counts.index, columns=counts.columns)

# Load the positions
positions = pd.read_table("../../Puck_180819_12/positions/positions.txt", index_col=0)
positions = positions.loc[norm_counts.columns]


# %% For a given set of genes, form a clustermap

from scipy.cluster.hierarchy import leaves_list, linkage, cut_tree
from scipy.spatial.distance import pdist, squareform

def plot_comparison(method_A, method_B, n_clusters):

    genes = selected_genes[method_A].difference(selected_genes[method_B])

    norm_counts_sub = norm_counts.loc[genes]


    corr = squareform(pdist(norm_counts_sub, metric='correlation'))
    corr = 1-corr

    # cluster the rows of corr

    cluster_cmap = plt.get_cmap('tab10').colors

    Z = linkage(corr, metric='euclidean', method='ward')
    corr_ii = leaves_list(Z)

    clusters = cut_tree(Z, n_clusters=n_clusters).ravel()

    cluster_colors = [cluster_cmap[i] for i in clusters]

    # Let seaborn cluster things too so we get the dendrogram
    sns.clustermap(
        corr[:, corr_ii], vmin=-.01, vmax=.01, cmap='RdBu_r', metric='euclidean', method='ward',
        row_colors=cluster_colors, xticklabels=False, yticklabels=False, col_cluster=False, rasterized=True,
        figsize=(6, 6)
    )

    plt.xlabel('')
    plt.ylabel('')
    title = 'Genes in top 1000 by {} but not by {}\n (N={})'.format(method_A, method_B, len(genes))
    plt.suptitle(title)
    # plt.show()
    plt.savefig('moran_clustermap_{}_{}.png'.format(method_A, method_B), dpi=300)
    plt.close()

    # %% Also, plot the positions

    fig, axs = plt.subplots(
        2, 3, gridspec_kw={'height_ratios': [1, 12]}, figsize=(9, 3.5)
    )

    for i in range(3):

        plt.sca(axs[0, i])
        color = cluster_cmap[i]
        plt.imshow([[color]], aspect='auto')
        sns.despine(top=True, right=True, bottom=True, left=True)
        plt.xticks([])
        plt.yticks([])

        plt.sca(axs[1, i])

        vals = norm_counts_sub.loc[clusters == i].mean(axis=0)
        vmin, vmax = np.percentile(vals, [2, 98])
        # vmin, vmax = 0, 5

        plt.scatter(positions.Comp1, positions.Comp2, s=1, c=vals, vmin=vmin, vmax=vmax, rasterized=True, edgecolors='none')
        plt.xticks([])
        plt.yticks([])

    plt.suptitle(title)
    plt.savefig('moran_slides_{}_{}.svg'.format(method_A, method_B), dpi=300)
    plt.close()


# %%

comparison_list = [
    # (method A, method B, n_clusters, clusters_to_viz)
    ('Hotspot', 'PCA', 3, [0, 1, 2]),
    ('Hotspot', 'HVG', 3, [0, 1, 2]),
    ('Hotspot', 'NBDisp', 3, [0, 1, 2]),
    ('PCA', 'Hotspot', 3, [0, 1, 2]),
    ('HVG', 'Hotspot', 3, [0, 1, 2]),
    ('NBDisp', 'Hotspot', 3, [0, 1, 2]),
]

for method_A, method_B, n_clusters, cluster_to_viz in tqdm(comparison_list):
    plot_comparison(method_A, method_B, n_clusters)
