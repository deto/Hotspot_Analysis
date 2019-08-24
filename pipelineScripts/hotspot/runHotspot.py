import loompy
import pandas as pd
import hotspot

loom_file = snakemake.input['loom']
latent_file = snakemake.input['latent']
out_file = snakemake.output['results']

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

hs = hotspot.Hotspot(counts, latent, num_umi)

hs.create_knn_graph(
    weighted_graph=False, n_neighbors=30, neighborhood_factor=3
)

results = hs.compute_hotspot(model=model, jobs=2, centered=True)

results = gene_info.join(results, how='right')

results.to_csv(out_file, sep="\t")
