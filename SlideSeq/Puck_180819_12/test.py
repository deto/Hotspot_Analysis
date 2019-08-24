import os
import numpy as np
import pandas as pd
import loompy

import NaiveDE
import SpatialDE

import sys


def trace(frame, event, arg):
    if event != "line":
        print("%s, %s:%d" % (event, frame.f_code.co_filename, frame.f_lineno))
    return trace


sys.settrace(trace)

# Load some data

print('Loading Data...')
# latent_file = snakemake.input['latent']
latent_file = "positions/positions.txt"
latent = pd.read_csv(latent_file, index_col=0, sep="\t")

loom_file = "../data/SlideSeq/data.loom"

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


print('Running SpatialDE...')
results = SpatialDE.run(X, resid_expr)

print('Saving Results...')
results = results.set_index('g')
#results.to_csv(out_file, sep="\t")
print('Done!')



# -------------------------------------------
# This part to reproduce the error I was getting (segfault).  Appears to come from linalg.eigh
# -------------------------------------------
from SpatialDE.base import get_l_limits
#l_min, l_max = get_l_limits(X)
l_min, l_max = (.25495, 10423.78)
kernel_space = {
    'SE': np.logspace(np.log10(l_min), np.log10(l_max), 10),
    'const': 0
}
US_mats = []

# for lengthscale in kernel_space['SE']:
#     K = SE_kernel(X, lengthscale)
#     U, S = factor(K)
#     gower = gower_scaling_factor(K)
#     UT1 = get_UT1(U)
#     US_mats.append({
#         'model': 'SE',
#         'M': 4,
#         'l': lengthscale,
#         'U': U,
#         'S': S,
#         'UT1': UT1,
#         'Gower': gower
#     })

lengthscale = kernel_space['SE'][0]
from SpatialDE.base import SE_kernel

K = SE_kernel(X, lengthscale)
S, U = np.linalg.eigh(K)

# More minimal example?
import numpy as np

X = np.random.randn(30000, 30000)

out = np.linalg.eigh(X)
