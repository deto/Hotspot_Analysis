import numpy as np
import pandas as pd
import hotspot.modules
import loompy
import hotspot
import hotspot.modules
from numba import njit
from tqdm import tqdm


loom_file = snakemake.input["loom"]
cm_file = snakemake.input['cm_file']  # Character-matrix file
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

cm = pd.read_table(cm_file, index_col=0)
modules = pd.read_table(module_file, index_col=0).Cluster

# Have to do this because data_slideseq makes it a numpy array
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)

if use_umi:
    num_umi = pd.Series(num_umi, index=barcodes)
else:
    num_umi = pd.Series(1.0, index=barcodes)

# Align to cm index
counts = counts.loc[:, cm.index]
num_umi = num_umi[cm.index]

# Compute distance matrix

cm[cm == '-'] = float('nan')
for x in cm.columns:
    cm[x] = cm[x].astype('float64')


@njit
def dist_fun(r1, r2):
    N = len(r1)
    tot = 0.0
    recovered = 0.0
    for i in range(N):
        x = r1[i]
        y = r2[i]

        if np.isnan(x) or np.isnan(y):
            continue

        recovered += 1

        if x > 0 and y == 0:
            tot += 1
        if x == 0 and y > 0:
            tot += 1
        if x > 0 and y > 0 and x != y:
            tot += 2

    if recovered == 0:
        return 2

    return tot/recovered


dist_mat = np.zeros((cm.shape[0], cm.shape[0]))

cm_np = cm.values

for i in tqdm(range(cm.shape[0])):
    for j in range(cm.shape[0]):

        dd = dist_fun(cm_np[i, :], cm_np[j, :])
        dist_mat[i, j] = dd

dist_mat = pd.DataFrame(dist_mat, index=cm.index, columns=cm.index)

# need counts, distances, and num_umi

hs = hotspot.Hotspot(counts, distances=dist_mat, umi_counts=num_umi)
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
