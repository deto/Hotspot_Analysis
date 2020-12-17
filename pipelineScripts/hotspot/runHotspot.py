import loompy
import pandas as pd
import hotspot

loom_file = snakemake.input['loom']
latent_file = snakemake.input['latent']
out_file = snakemake.output['results']

model = snakemake.params['model']


try:
    n_neighbors = int(snakemake.params['n_neighbors'])
except AttributeError:
    n_neighbors = 30

try:
    n_cells_min = snakemake.params['n_cells_min']
except AttributeError:
    n_cells_min = 50

try:
    use_umi = bool(snakemake.params['use_umi'])
except AttributeError:
    use_umi = True

try:
    weighted_graph = snakemake.params['weighted_graph']
except AttributeError:
    weighted_graph = False

try:
    n_bins_bernoulli = int(snakemake.params['n_bins_bernoulli'])
    from hotspot import bernoulli_model
    bernoulli_model.N_BIN_TARGET = n_bins_bernoulli
except AttributeError:
    pass

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

if use_umi:
    num_umi = pd.Series(num_umi, index=barcodes)
else:
    num_umi = pd.Series(1.0, index=barcodes)

# Align to latent space
counts = counts.loc[:, latent.index]
num_umi = num_umi[latent.index]

# need counts, latent, and num_umi

valid_genes = (counts > 0).sum(axis=1) >= n_cells_min
counts = counts.loc[valid_genes]

hs = hotspot.Hotspot(counts, model=model, latent=latent, umi_counts=num_umi)

hs.create_knn_graph(
    weighted_graph=weighted_graph, n_neighbors=n_neighbors, neighborhood_factor=3
)

results = hs.compute_hotspot(jobs=5)

results = gene_info.join(results, how='right')

results.to_csv(out_file, sep="\t")
