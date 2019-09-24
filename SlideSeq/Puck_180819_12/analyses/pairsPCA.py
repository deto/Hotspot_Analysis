"""
Take the z-score matrix and run an SVD on it

See if this looks better than normal PCA + tSNE
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

pairs = pd.read_table("../hotspot/hotspot_pairs_z.txt.gz", index_col=0)
positions = pd.read_table("../positions/positions.txt", index_col=0)

hs_genes = pairs.index

# %% Load expression matrix and subset

import loompy

loom_file = "../../../data/SlideSeq/Puck_180819_12/data.loom"

with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    scaled = ds.layers['scaled'][:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)
scaled = pd.DataFrame(scaled, index=gene_info.index, columns=barcodes)
num_umi = pd.Series(num_umi, index=barcodes)
counts_hs = counts.loc[hs_genes]
scaled_hs = scaled.loc[hs_genes]

from hotspot.local_stats_pairs import create_centered_counts
model = 'bernoulli'
counts_hs_centered = create_centered_counts(counts_hs.values, model, num_umi.values)
counts_hs_centered = pd.DataFrame(
    counts_hs_centered, index=counts_hs.index, columns=counts_hs.columns)

# %% SVD

from sklearn.decomposition import TruncatedSVD

model = TruncatedSVD(n_components=5)

results = model.fit_transform(counts_hs_centered.values.T)
results.shape

# %% TSNE

from sklearn.manifold import TSNE

model_tsne = TSNE(n_components=2, perplexity=30, verbose=10)

res_tsne = model_tsne.fit_transform(results)

# %% Plot


plt.figure()
plt.plot(res_tsne[:, 0], res_tsne[:, 1], 'o', ms=2)
plt.show()

plt.figure()
plt.scatter(res_tsne[:, 0], res_tsne[:, 1], s=2, c=counts_hs_centered.loc['Car8'], vmin=-3, vmax=3)
plt.show()

plt.figure()
vals = counts_hs_centered.loc['Car8']
plt.scatter(
    positions.Comp1, positions.Comp2,
    s=1, c=vals,
    vmin=np.percentile(vals, 5),
    vmax=np.percentile(vals, 95))
plt.show()

plt.figure()
plt.scatter(positions.Comp1, positions.Comp2, s=2, c=counts_hs_centered.loc['Fth1'], vmin=0, vmax=3)
plt.show()



# What if we just do the log-counts PCA but on the reduced features?

# %% PCA

from sklearn.decomposition import PCA

model = PCA(n_components=15)

results = model.fit_transform(np.log2(scaled_hs.T+1))
results.shape

# %% TSNE

from sklearn.manifold import TSNE

model_tsne = TSNE(n_components=2, perplexity=30, verbose=10)

res_tsne = model_tsne.fit_transform(results)
