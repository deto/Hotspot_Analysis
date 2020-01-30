import loompy
import numpy as np
import pandas as pd
import hotspot
from numba import njit
from tqdm import tqdm

loom_file = snakemake.input['loom']
cm_file = snakemake.input['cm_file']  # Character-matrix file
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

try:
    weighted_graph = snakemake.params['weighted_graph']
except AttributeError:
    weighted_graph = False

try:
    use_umi = bool(snakemake.params['use_umi'])
except AttributeError:
    use_umi = True

with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

cm = pd.read_table(cm_file, index_col=0)

# Have to do this because data_slideseq makes it a numpy array
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)

if use_umi:
    num_umi = pd.Series(num_umi, index=barcodes)
else:
    num_umi = pd.Series(1.0, index=barcodes)


# Align to latent space
counts = counts.loc[:, cm.index]
num_umi = num_umi[cm.index]

# need counts, latent, and num_umi

valid_genes = (counts > 0).sum(axis=1) >= n_cells_min
counts = counts.loc[valid_genes]

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

hs = hotspot.Hotspot(counts, distances=dist_mat, umi_counts=num_umi)

hs.create_knn_graph(
    weighted_graph=weighted_graph, n_neighbors=n_neighbors,
    neighborhood_factor=3
)

results = hs.compute_hotspot(model=model, jobs=5, centered=True)

results = gene_info.join(results, how='right')

results.to_csv(out_file, sep="\t")
