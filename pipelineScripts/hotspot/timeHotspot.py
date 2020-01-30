import time
import os
import loompy
import numpy as np
import pandas as pd
import hotspot

loom_file = snakemake.input['loom']
latent_file = snakemake.input['latent']
out_file = snakemake.output['results']

os.makedirs(os.path.dirname(out_file), exist_ok=True)

N_CELLS = int(snakemake.params['N_CELLS'])
N_GENES = int(snakemake.params['N_GENES'])
N_NEIGHBORS = int(snakemake.params['n_neighbors'])
N_JOBS = snakemake.threads

model = snakemake.params['model']

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

valid_genes = (counts > 0).sum(axis=1) > 50
counts = counts.loc[valid_genes]

# Subsample cells

cells_sub = np.random.choice(counts.shape[1], N_CELLS, replace=False)
counts = counts.iloc[:, cells_sub]
latent = latent.iloc[cells_sub, :]
num_umi = num_umi.iloc[cells_sub]

# Subsample Genes
# Need to samples genes that are expressed in the subset
valid_genes = np.nonzero((counts > 0).sum(axis=1) > 5)[0]
genes_sub = np.random.choice(valid_genes, N_GENES, replace=False)
counts = counts.iloc[genes_sub, :]

start = time.time()
hs = hotspot.Hotspot(counts, latent=latent, umi_counts=num_umi)

hs.create_knn_graph(
    weighted_graph=False, n_neighbors=N_NEIGHBORS, neighborhood_factor=3
)

results = hs.compute_hotspot(model=model, jobs=N_JOBS, centered=True)

stop = time.time()

with open(out_file, 'w') as fout:
    fout.write("Genes\tCells\tElapsedSeconds\n")
    fout.write("{}\t{}\t{:.3f}\n".format(
        counts.shape[0],
        counts.shape[1],
        stop-start,
    ))
