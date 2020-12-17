import loompy
import pandas as pd
import hotspot
import hotspot.knn
from ete3 import Tree

loom_file = snakemake.input['loom']
tree_file = snakemake.input['tree_file']
hs_results_file = snakemake.input['hs_results']

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
hs_results = pd.read_table(hs_results_file, index_col=0)

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

# need counts, distances, and num_umi

latent = pd.DataFrame(0, index=counts.columns, columns=range(10))
hs = hotspot.Hotspot(counts, model=model, latent=latent, umi_counts=num_umi)

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

lcz = hs.compute_local_correlations(hs_genes, jobs=20)

lcz.to_csv(out_file_lcz, sep="\t", compression="gzip")
