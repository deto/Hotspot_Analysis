"""
This file is no longer needed as the methods here have been
fully integrated into the Hotspot python code-base
"""

# %% Define some methods

from numba import jit
import numpy as np


@jit(nopython=True)
def conditional_eg2(x, neighbors, weights):
    """
    Computes EG2 for the conditional correlation
    """
    out_j = np.zeros(len(x))

    for i in range(len(x)):

        xi = x[i]

        for K in range(neighbors.shape[1]):

            j = neighbors[i, K]
            w_ij = weights[i, K]

            out_j[j] += (w_ij * xi)

    out_eg2 = (out_j**2).sum()

    return out_eg2


@jit(nopython=True)
def local_cov_pair_cond(x, y, neighbors, weights):
    """
    Test statistic for local pair-wise autocorrelation
    Computes the statistic assuming x values are fixed and
    y values are allowed to vary.
    """
    out_xy = 0
    out_yx = 0

    for i in range(len(x)):
        for k in range(neighbors.shape[1]):

            j = neighbors[i, k]
            w_ij = weights[i, k]

            xi = x[i]
            xj = x[j]

            yi = y[i]
            yj = y[j]

            out_xy += w_ij*xi*yj
            out_yx += w_ij*xj*yi

    return out_xy, out_yx


@jit(nopython=True)
def _compute_hs_pairs_inner_centered_cond_sym(
    rowpair, counts, neighbors, weights, eg2s
):
    """
    This version assumes that the counts have already been modeled
    and centered
    """
    row_i, row_j = rowpair

    vals_x = counts[row_i]
    vals_y = counts[row_j]

    lc_xy, lc_yx = local_cov_pair_cond(vals_x, vals_y, neighbors, weights)

    # Compute xy
    EG, EG2 = 0, eg2s[row_i]

    stdG = (EG2 - EG ** 2) ** 0.5

    Zxy = (lc_xy - EG) / stdG

    # Compute yx
    EG, EG2 = 0, eg2s[row_j]

    stdG = (EG2 - EG ** 2) ** 0.5

    Zyx = (lc_yx - EG) / stdG

    if abs(Zxy) < abs(Zyx):
        Z = Zxy
    else:
        Z = Zyx

    lc = (lc_xy + lc_yx) / 2

    return (lc, Z)






# %% Actual pipeline code
import loompy
import pandas as pd
import hotspot
from hotspot import local_stats_pairs
from tqdm import tqdm

loom_file = snakemake.input['loom']
latent_file = snakemake.input['latent']
hs_results_file = snakemake.input['hs_results']

out_file_lc = snakemake.output['results_lc']
out_file_lcz = snakemake.output['results_z']

model = snakemake.params['model']

fdrThresh = snakemake.params['fdrThresh']
n_neighbors = snakemake.params['n_neighbors']

try:
    topN = int(snakemake.params['topN'])
except AttributeError:
    topN = None

try:
    highXMeanCutoff = float(snakemake.params['highXMeanCutoff'])
except AttributeError:
    highXMeanCutoff = None

with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

latent = pd.read_table(latent_file, index_col=0)
hs_results = pd.read_table(hs_results_file, index_col=0)

try:
    genes_file = snakemake.input['genes']
except AttributeError:
    genes_file = None

# Have to do this because data_slideseq makes it a numpy array
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)
num_umi = pd.Series(num_umi, index=barcodes)

# Align to latent space
counts = counts.loc[:, latent.index]
num_umi = num_umi[latent.index]

# need counts, latent, and num_umi

hs = hotspot.Hotspot(counts, latent=latent, umi_counts=num_umi)

hs.create_knn_graph(
    weighted_graph=False, n_neighbors=n_neighbors, neighborhood_factor=3
)

if highXMeanCutoff is not None:

    scaled = counts.divide(counts.sum(axis=0), axis=1)*10000
    gene_means = scaled.mean(axis=1)
    valid_genes = gene_means.index[gene_means < highXMeanCutoff]
    hs_results = hs_results.loc[valid_genes & hs_results.index]


if topN is None:
    hs_genes = hs_results.index[hs_results.FDR < fdrThresh]
else:
    hs_genes = hs_results.sort_values('Z').tail(topN).index

# Overrides the hs results if desired
if genes_file is not None:
    hs_genes = pd.Index(pd.read_table(genes_file, header=None).iloc[:, 0].tolist())

hs_genes = hs_genes & counts.index

counts = counts.loc[hs_genes]

c_counts = local_stats_pairs.create_centered_counts(
    counts.values, model, num_umi.values)

eg2s = np.array(
    [
        conditional_eg2(c_counts[i], hs.neighbors.values, hs.weights.values)
        for i in range(c_counts.shape[0])
    ]
)

def _map_fun_parallel_pairs_centered(rowpair):
    global g_neighbors
    global g_weights
    global g_counts
    global g_eg2s
    return _compute_hs_pairs_inner_centered_cond_sym(
        rowpair, g_counts, g_neighbors, g_weights, eg2s
    )


def initializer():
    global g_neighbors
    global g_weights
    global g_counts
    global g_eg2s
    g_counts = c_counts
    g_neighbors = hs.neighbors.values
    g_weights = hs.weights.values
    g_eg2s = eg2s

import itertools
import multiprocessing

jobs = 20

pairs = list(itertools.combinations_with_replacement(range(c_counts.shape[0]), 2))

with multiprocessing.Pool(
        processes=jobs, initializer=initializer) as pool:

    results = list(
        tqdm(
            pool.imap(_map_fun_parallel_pairs_centered, pairs),
            total=len(pairs)
        )
    )

N = counts.shape[0]
pairs = np.array(pairs)
vals_lc = np.array([x[0] for x in results])
vals_z = np.array([x[1] for x in results])
lcs = local_stats_pairs.expand_pairs(pairs, vals_lc, N)
lc_zs = local_stats_pairs.expand_pairs(pairs, vals_z, N)

res_Z = pd.DataFrame(
    lc_zs, index=counts.index, columns=counts.index
)

res_lc = pd.DataFrame(
    lcs, index=counts.index, columns=counts.index
)


res_lc.to_csv(out_file_lc, sep="\t", compression="gzip")
res_Z.to_csv(out_file_lcz, sep="\t", compression="gzip")
