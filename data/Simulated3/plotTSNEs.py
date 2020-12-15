import numpy as np
import pandas as pd
from sklearn.decomposition import PCA
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from umap import UMAP

# %% Load data

counts_file = "rep8/observed_counts.txt.gz"
cell_effects_file = "rep8/cell_effects.txt"

counts = pd.read_table(counts_file, index_col=0)
cell_effects = pd.read_table(cell_effects_file, index_col=0)


# %% Run PCA then tSNE
norm_counts = counts.divide(counts.sum(axis=0), axis=1)
norm_counts = np.log2(norm_counts+1)
norm_counts = norm_counts \
    .subtract(norm_counts.mean(axis=1), axis=0) \
    .divide(norm_counts.std(axis=1), axis=0)

norm_counts = norm_counts.loc[
    ~norm_counts.isnull().any(axis=1), :]

pca = PCA(n_components=20)

pca.fit(norm_counts.T)
pc_vals = pca.transform(norm_counts.T)


tsne = TSNE(n_components=2)
tsne_vals = tsne.fit_transform(pc_vals[:, 0:5])

umap = UMAP(n_components=2, n_neighbors=30)
umap_vals = umap.fit_transform(pc_vals[:, 0:5])

# %% Plot tSNE

fig, axs = plt.subplots(1, 5, figsize=(13, 4))

for i in range(5):
    plt.sca(axs[i])
    ce = cell_effects.iloc[:, i+1]
    vmin = np.percentile(ce, 5)
    vmax = np.percentile(ce, 95)
    if vmin == vmax:
        vmin = np.percentile(ce, .5)
        vmax = np.percentile(ce, 99.5)
    plt.scatter(tsne_vals[:, 0], tsne_vals[:, 1], s=2, c=ce, vmin=vmin, vmax=vmax)
    #plt.scatter(umap_vals[:, 0], umap_vals[:, 1], s=2, c=ce, vmin=vmin, vmax=vmax)
    plt.xticks([])
    plt.yticks([])

plt.tight_layout()
plt.show()

# %% Plot PCA components on the tSNE

fig, axs = plt.subplots(1, 5, figsize=(13, 4))

for i in range(5):
    plt.sca(axs[i])
    ce = pc_vals[:, i]
    vmin = np.percentile(ce, 5)
    vmax = np.percentile(ce, 95)
    plt.scatter(tsne_vals[:, 0], tsne_vals[:, 1], s=2, c=ce, vmin=vmin, vmax=vmax)
    #plt.scatter(umap_vals[:, 0], umap_vals[:, 1], s=2, c=ce, vmin=vmin, vmax=vmax)
    plt.xticks([])
    plt.yticks([])

plt.tight_layout()
plt.show()

# %% Plot eigenvalue spectrum
plt.figure()
plt.bar(x=np.arange(pca.explained_variance_.size)+1, height=pca.explained_variance_)
plt.xticks(np.arange(pca.explained_variance_.size)+1)
plt.show()
