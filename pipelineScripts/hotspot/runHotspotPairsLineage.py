import loompy
import numpy as np
import pandas as pd
import hotspot
from numba import njit
from tqdm import tqdm

loom_file = snakemake.input['loom']
cm_file = snakemake.input['cm_file']  # Character-matrix file
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

cm = pd.read_table(cm_file, index_col=0)
hs_results = pd.read_table(hs_results_file, index_col=0)

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

hs = hotspot.Hotspot(counts, model=model, distances=dist_mat, umi_counts=num_umi)

hs.create_knn_graph(
    weighted_graph=False, n_neighbors=n_neighbors, neighborhood_factor=3
)

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
