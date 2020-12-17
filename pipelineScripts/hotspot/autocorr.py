import loompy
import numpy as np
import pandas as pd
import hotspot
from numba import jit
from tqdm import tqdm

loom_file = snakemake.input['loom']
latent_file = snakemake.input['latent']
out_file = snakemake.output['results']


@jit(nopython=True)
def geary(X, N, W):
    """
    X: values (length N)
    N: neighbor indices (length N x Kneighbors)
    W: neighbor weights (length N x Kneighbors)
    """

    wtot = 0.0
    num_tot = 0.0
    denom_tot = 0.0

    mu = X.mean()

    for i in range(N.shape[0]):

        xi = X[i]

        for k in range(N.shape[1]):

            j = N[i, k]
            xj = X[j]

            wij = W[i, k]

            num_tot += wij * (xi-xj)**2
            wtot += wij

        denom_tot += (xi-mu)**2

    M = N.shape[0]

    return (M-1) * num_tot / (2*wtot*denom_tot)


@jit(nopython=True)
def moran(X, N, W):
    """
    X: values (length N)
    N: neighbor indices (length N x Kneighbors)
    W: neighbor weights (length N x Kneighbors)
    """

    wtot = 0.0
    num_tot = 0.0
    denom_tot = 0.0

    mu = X.mean()

    for i in range(N.shape[0]):

        xi = X[i]

        for k in range(N.shape[1]):

            j = N[i, k]
            xj = X[j]

            wij = W[i, k]

            num_tot += wij * (xi - mu) * (xj - mu)
            wtot += wij

        denom_tot += (xi-mu)**2

    M = N.shape[0]

    return M/wtot * num_tot / denom_tot


try:
    n_neighbors = int(snakemake.params['n_neighbors'])
except AttributeError:
    n_neighbors = 30

try:
    n_cells_min = snakemake.params['n_cells_min']
except AttributeError:
    n_cells_min = 50

try:
    weighted_graph = snakemake.params['weighted_graph']
except AttributeError:
    weighted_graph = False

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

# need counts, latent, and num_umi

valid_genes = (counts > 0).sum(axis=1) >= n_cells_min
counts = counts.loc[valid_genes]

hs = hotspot.Hotspot(counts, model='danb', latent=latent, umi_counts=num_umi)

hs.create_knn_graph(
    weighted_graph=weighted_graph,
    n_neighbors=n_neighbors,
    neighborhood_factor=3
)

neighbors = hs.neighbors.values
weights = hs.weights.values

norm_counts = counts.values / num_umi.values.reshape((1, -1)) * 10000
norm_counts = np.log(norm_counts + 1)

from hotspot.local_stats_pairs import create_centered_counts

norm_counts = create_centered_counts(counts.values, hs.model, num_umi.values)


geary_results = []
moran_results = []
for i in tqdm(range(norm_counts.shape[0])):

    exp_vals = norm_counts[i, :]

    g = geary(exp_vals, neighbors, weights)
    geary_results.append(g)

    m = moran(exp_vals, neighbors, weights)
    moran_results.append(m)

results = pd.DataFrame(
    {
        'GearyC': geary_results,
        'MoranI': moran_results,
    }, index=counts.index
)

results = gene_info.join(results, how='right')

results.to_csv(out_file, sep="\t")
