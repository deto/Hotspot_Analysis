import ete3
from ete3 import Tree
import loompy
import numpy as np
import pandas as pd
import hotspot
import hotspot.knn
from numba import njit
from tqdm import tqdm

loom_file = snakemake.input['loom']
tree_file = snakemake.input['tree_file']
out_file = snakemake.output['results']

model = snakemake.params['model']


try:
    n_neighbors = snakemake.params['n_neighbors']
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

try:
    use_umi = bool(snakemake.params['use_umi'])
except AttributeError:
    use_umi = True

with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

t = Tree(tree_file, format=1)

# Have to do this because data_slideseq makes it a numpy array
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)

if use_umi:
    num_umi = pd.Series(num_umi, index=barcodes)
else:
    num_umi = pd.Series(1.0, index=barcodes)

valid_cells = set()
for x in t:
    if x.is_leaf():
        valid_cells.add(x.name)
valid_cells = pd.Index(valid_cells)

# Align to latent space
counts = counts.loc[:, valid_cells]
num_umi = num_umi[valid_cells]

# need counts, latent, and num_umi

valid_genes = (counts > 0).sum(axis=1) >= n_cells_min
counts = counts.loc[valid_genes]

# Compute distance matrix

latent = pd.DataFrame(0, index=counts.columns, columns=range(10))
hs = hotspot.Hotspot(counts, model=model, latent=latent, umi_counts=num_umi)

neighbors, weights = hotspot.knn.tree_neighbors_and_weights(
    t, n_neighbors, counts)


if weighted_graph:
    raise ValueError("Can't do weighted graph on tree")

weights = hotspot.knn.make_weights_non_redundant(
    neighbors.values, weights.values
)
weights = pd.DataFrame(weights, index=neighbors.index, columns=neighbors.columns)

hs.weights = weights
hs.neighbors = neighbors

results = hs.compute_hotspot(jobs=5)

results = gene_info.join(results, how='right')

results.to_csv(out_file, sep="\t")
