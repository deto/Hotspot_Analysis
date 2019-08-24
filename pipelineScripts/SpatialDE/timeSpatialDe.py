import time
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

os.makedirs(os.path.dirname(out_file), exist_ok=True)

N_CELLS = int(snakemake.params['N_CELLS'])
N_GENES = int(snakemake.params['N_GENES'])
N_JOBS = snakemake.threads

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

# Subsample cells
cells_sub = np.random.choice(counts.shape[0], N_CELLS, replace=False)
counts = counts.iloc[cells_sub, :]
sample_info = sample_info.iloc[cells_sub, :]

# Subsample Genes
# Need to samples genes that are expressed in the subset
valid_genes = np.nonzero((counts > 0).sum(axis=0) > 5)[0]
genes_sub = np.random.choice(valid_genes, N_GENES, replace=False)
counts = counts.iloc[:, genes_sub]

# Some normalization

t1 = time.time()

print('Normalize Data...')
norm_expr = NaiveDE.stabilize(counts.T).T

resid_expr = NaiveDE.regress_out(
    sample_info, norm_expr.T, 'np.log(total_counts)').T


X = sample_info[['Comp1', 'Comp2']].values

t2 = time.time()

print('Running SpatialDE...')
results = SpatialDE.run(X, resid_expr)

t3 = time.time()

with open(out_file, 'w') as fout:
    fout.write("Genes\tCells\tNormSeconds\tSDESeconds\n")
    fout.write("{}\t{}\t{:.3f}\t{:.3f}\n".format(
        counts.shape[1],
        counts.shape[0],
        t2-t1,
        t3-t2,
    ))
