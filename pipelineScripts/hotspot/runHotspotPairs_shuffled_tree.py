# %% Define some methods

from numba import jit
import numpy as np

from hotspot.local_stats_pairs import local_cov_pair, conditional_eg2

@jit(nopython=True)
def _compute_hs_pairs_inner_centered_cond_sym_shuff(
    rowpair, counts, neighbors, weights, eg2s
):
    """
    This version assumes that the counts have already been modeled
    and centered
    """
    row_i, row_j = rowpair

    vals_x = counts[row_i]
    vals_y = counts[row_j]

    vals_y = np.random.permutation(vals_y)

    lc = local_cov_pair(vals_x, vals_y, neighbors, weights)*2

    # Compute xy
    EG, EG2 = 0, eg2s[row_i]

    stdG = (EG2 - EG ** 2) ** 0.5

    Zxy = (lc - EG) / stdG

    # Compute yx
    eg2s_y_shuff = conditional_eg2(vals_y, neighbors, weights)
    EG, EG2 = 0, eg2s_y_shuff

    stdG = (EG2 - EG ** 2) ** 0.5

    Zyx = (lc - EG) / stdG

    if abs(Zxy) < abs(Zyx):
        Z = Zxy
    else:
        Z = Zyx

    return (lc, Z)


@jit(nopython=True)
def expand_pairs(pairs, vals, N):

    out = np.zeros((N, N))

    for i in range(len(pairs)):

        x = pairs[i, 0]
        y = pairs[i, 1]
        v = vals[i]

        out[x, y] = v

    return out


# %% Actual pipeline code
import loompy
import pandas as pd
import hotspot
from hotspot import local_stats_pairs
import hotspot.knn
from ete3 import Tree
from tqdm import tqdm

loom_file = snakemake.input['loom']
tree_file = snakemake.input['tree_file']
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

t = Tree(tree_file, format=1)
hs_results = pd.read_table(hs_results_file, index_col=0)

# Have to do this because data_slideseq makes it a numpy array
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)
num_umi = pd.Series(num_umi, index=barcodes)

valid_cells = set()
for x in t:
    if x.is_leaf():
        valid_cells.add(x.name)
valid_cells = pd.Index(valid_cells)

# Align to latent space
counts = counts.loc[:, valid_cells]
num_umi = num_umi[valid_cells]

# need counts, latent, and num_umi

latent = pd.DataFrame(0, index=counts.columns, columns=range(10))
hs = hotspot.Hotspot(counts, latent=latent, umi_counts=num_umi)

neighbors, weights = hotspot.knn.tree_neighbors_and_weights(
    t, n_neighbors, counts)

weights = hotspot.knn.make_weights_non_redundant(
    neighbors.values, weights.values
)
weights = pd.DataFrame(
    weights, index=neighbors.index, columns=neighbors.columns)

hs.weights = weights
hs.neighbors = neighbors

if highXMeanCutoff is not None:

    scaled = counts.divide(counts.sum(axis=0), axis=1)*10000
    gene_means = scaled.mean(axis=1)
    valid_genes = gene_means.index[gene_means < highXMeanCutoff]
    hs_results = hs_results.loc[valid_genes & hs_results.index]


if topN is None:
    hs_genes = hs_results.index[hs_results.FDR < fdrThresh]
else:
    hs_genes = hs_results.sort_values('Z').tail(topN).index

hs_genes = hs_genes & counts.index

counts = counts.loc[hs_genes]

c_counts = local_stats_pairs.create_centered_counts(
    counts.values, model, num_umi.values)

# Run the process

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
    return _compute_hs_pairs_inner_centered_cond_sym_shuff(
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

pairs = list(itertools.product(
    range(c_counts.shape[0]),
    range(c_counts.shape[0]),
))

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
lcs = expand_pairs(pairs, vals_lc, N)
lc_zs = expand_pairs(pairs, vals_z, N)

new_index = counts.index.tolist()

res_Z = pd.DataFrame(
    lc_zs, index=new_index, columns=new_index,
)

res_lc = pd.DataFrame(
    lcs, index=new_index, columns=new_index,
)


res_lc.to_csv(out_file_lc, sep="\t", compression="gzip")
res_Z.to_csv(out_file_lcz, sep="\t", compression="gzip")
