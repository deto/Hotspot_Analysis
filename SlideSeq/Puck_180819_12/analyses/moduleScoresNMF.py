"""
Going to test whether NMF works for module scores

It can work, but when comparing on other modules it does not
appear to perform as well as the existing PCA-based implementation
"""

import numpy as np
import pandas as pd
import loompy
import matplotlib.pyplot as plt

import hotspot
from hotspot.local_stats_pairs import create_centered_counts_row
from hotspot.utils import neighbor_smoothing_row

from sklearn.decomposition import PCA, NMF

# %% Load data
loom_file = "../../data/SlideSeq/Puck_180819_12/data.loom"
latent_file = "positions/positions.txt"
module_file = "hotspot/modules.txt"

n_neighbors = 30
gene_model = 'bernoulli'

with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

latent = pd.read_table(latent_file, index_col=0)
modules = pd.read_table(module_file, index_col=0).Cluster

# Have to do this because data_slideseq makes it a numpy array
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)
num_umi = pd.Series(num_umi, index=barcodes)

# Align to latent space
counts = counts.loc[:, latent.index]
num_umi = num_umi[latent.index]

# need counts, latent, and num_umi

hs = hotspot.Hotspot(counts, latent, num_umi)
hs.create_knn_graph(
    weighted_graph=False, n_neighbors=n_neighbors, neighborhood_factor=3
)


# %%

def neighbor_impute(counts_sub, neighbors, weights):
    cc_smooth = np.zeros_like(counts_sub, dtype=np.float64)

    for i in range(counts_sub.shape[0]):

        counts_row = counts_sub[i, :]
        smooth_row = neighbor_smoothing_row(
            counts_row, neighbors, weights, _lambda=.9)

        cc_smooth[i] = smooth_row

    return cc_smooth


def compute_scores_pca(
        counts_sub, model, num_umi, neighbors, weights):
    """
    counts_sub: row-subset of counts matrix with genes in the module
    """

    cc_smooth = np.zeros_like(counts_sub, dtype=np.float64)

    for i in range(counts_sub.shape[0]):

        counts_row = counts_sub[i, :]
        centered_row = create_centered_counts_row(counts_row, model, num_umi)
        smooth_row = neighbor_smoothing_row(
            centered_row, neighbors, weights, _lambda=.9)

        cc_smooth[i] = smooth_row

    pca_data = cc_smooth

    model = PCA(n_components=1)
    scores = model.fit_transform(pca_data.T)
    gene_scores = model.components_.ravel()

    sign = model.components_.mean()  # may need to flip
    if sign < 0:
        scores = scores * -1

    scores = scores[:, 0]

    return scores, gene_scores

# %%


module = 0

module_genes = modules.index[modules == module]

scores, gene_scores = compute_scores_pca(
    counts.loc[module_genes].values, gene_model, num_umi.values,
    hs.neighbors.values, hs.weights.values
)

gene_scores = pd.Series(
    gene_scores, index=module_genes
)


# %% Try NMF
c_sub = counts.loc[module_genes] \
    .astype('float64') \
    .divide(num_umi, axis=1).values * 50

c_sub = neighbor_impute(c_sub, hs.neighbors.values, hs.weights.values)

c_sub = np.log2(c_sub+1)

# Remove the mean component?
c_sub = c_sub / c_sub.mean(axis=1, keepdims=True)

model = NMF(n_components=1)
model.fit(c_sub.T)

scores_NMF = model.transform(c_sub.T).ravel()
gene_scores_NMF = pd.Series(
    model.components_.ravel(), index=module_genes
)

# %% Plot comparison

fig, axs = plt.subplots(1, 2, figsize=(9, 4))

all_scores = {
    'PCA': scores,
    'NMF': scores_NMF,
}

for ax, method in zip(axs.ravel(), all_scores):
    plt.sca(ax)
    sc = all_scores[method]

    vmin = np.percentile(sc, 1)
    vmax = np.percentile(sc, 99)

    plt.scatter(
        latent.iloc[:, 0],
        latent.iloc[:, 1],
        c=sc, cmap='viridis', s=2,
        vmin=vmin, vmax=vmax
    )

    plt.xticks([])
    plt.yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)

    plt.title(method)

plt.show()

# %% Compare the gene scores

from bio_utils import hover_plot
fig, ax = plt.subplots(1, 1)

hover_plot(
    gene_scores,
    gene_scores_NMF.loc[gene_scores.index],
    gene_scores.index,
    'o', ms=2, ax=ax
)
plt.show()


# %% compare cell scores
plt.figure()
plt.plot(
    scores, scores_NMF, 'o', ms=2
)
plt.show()

plt.figure()
plt.hist(scores, 30)
plt.show()

plt.figure()
plt.hist(scores_NMF, 30)
plt.show()
