import numpy as np
import pandas as pd
import hotspot.modules
import loompy
import hotspot
import hotspot.modules


loom_file = snakemake.input["loom"]
latent_file = snakemake.input["latent"]
module_file = snakemake.input["modules"]

n_neighbors = snakemake.params['n_neighbors']
model = snakemake.params['model']

out_file = snakemake.output['scores']

with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

latent = pd.read_table(latent_file, index_col=0)
modules = pd.read_table(module_file, index_col=0).Cluster

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

# %% Plot scores for all modules
modules_to_compute = sorted([x for x in modules.unique() if x != -1])

# Get the scores
module_scores = {}
for module in modules_to_compute:
    module_genes = modules.index[modules == module]

    scores = hotspot.modules.compute_scores(
        counts.loc[module_genes].values, model, num_umi.values,
        hs.neighbors.values, hs.weights.values
    )

    module_scores[module] = scores

module_scores = pd.DataFrame(module_scores)
module_scores.index = counts.columns

module_scores.to_csv(out_file, sep="\t")
