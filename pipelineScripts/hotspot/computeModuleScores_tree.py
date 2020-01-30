import pandas as pd
import hotspot.modules
import loompy
import hotspot
import hotspot.modules

from ete3 import Tree


loom_file = snakemake.input["loom"]
tree_file = snakemake.input['tree_file']
module_file = snakemake.input["modules"]

n_neighbors = snakemake.params['n_neighbors']
model = snakemake.params['model']

try:
    use_umi = bool(snakemake.params['use_umi'])
except AttributeError:
    use_umi = True

out_file = snakemake.output['scores']

with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

t = Tree(tree_file, format=1)
modules = pd.read_table(module_file, index_col=0).Cluster

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
hs = hotspot.Hotspot(counts, latent=latent, umi_counts=num_umi)

neighbors, weights = hotspot.knn.tree_neighbors_and_weights(
    t, n_neighbors, counts)

weights = hotspot.knn.make_weights_non_redundant(
    neighbors.values, weights.values
)
weights = pd.DataFrame(
    weights, index=neighbors.index, columns=neighbors.columns)

hs.weights = weights
hs.neighbors = neighbors

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
