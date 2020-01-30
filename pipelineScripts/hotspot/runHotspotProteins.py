import loompy
import numpy as np
import pandas as pd
import hotspot

proteins_file = snakemake.input['proteins']
latent_file = snakemake.input['latent']
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

latent = pd.read_table(latent_file, index_col=0)

counts = pd.read_table(proteins_file, index_col=0).T

num_umi = counts.sum(axis=0)

counts = np.log2(counts + 1)
num_umi = np.log2(num_umi)

# num_umi = pd.Series(1.0, index=counts.columns)

# Align to latent space
counts = counts.loc[:, latent.index]
num_umi = num_umi[latent.index]

# need counts, latent, and num_umi
hs = hotspot.Hotspot(counts, latent=latent, umi_counts=num_umi)

hs.create_knn_graph(
    weighted_graph=False, n_neighbors=n_neighbors, neighborhood_factor=3
)

results = hs.compute_hotspot(model=model, jobs=5, centered=True)

results.to_csv(out_file, sep="\t")
