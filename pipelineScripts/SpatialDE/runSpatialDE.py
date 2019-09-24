import os
import numpy as np
import pandas as pd
import loompy

import NaiveDE
import SpatialDE


# Load some data

print('Loading Data...')
latent_file = snakemake.input['latent']
latent = pd.read_csv(latent_file, index_col=0, sep="\t")

loom_file = snakemake.input['loom']
out_file = snakemake.output['results']

with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')

print('Filtering Data...')
# Do some filtering
valid_genes = (counts > 0).sum(axis=1) > 50
counts = counts[valid_genes, :]
gene_info = gene_info.loc[valid_genes]

counts = counts.T

counts = pd.DataFrame(counts, columns=gene_info.index)

# prepare 'sample_info' like they want it
sample_info = latent.copy()
sample_info['total_counts'] = num_umi

# Some normalization


print('Normalize Data...')
norm_expr = NaiveDE.stabilize(counts.T).T

resid_expr = NaiveDE.regress_out(
    sample_info, norm_expr.T, 'np.log(total_counts)').T

X = sample_info[['Comp1', 'Comp2']].values

# Override the ranges that spatialDE uses
# These mess up on some of the pucks and result in ridiculous values
from SpatialDE.base import get_l_limits
l_min, l_max = get_l_limits(X)

try:
    l_min = snakemake.params['l_min']
except AttributeError:
    pass

print(l_min, l_max)
kernel_space = {
    'SE': np.logspace(np.log10(l_min), np.log10(l_max), 10),
    'const': 0
}


print('Running SpatialDE...')
results = SpatialDE.run(X, resid_expr, kernel_space=kernel_space)

print('Saving Results...')
results = results.set_index('g')
results.to_csv(out_file, sep="\t")
print('Done!')
