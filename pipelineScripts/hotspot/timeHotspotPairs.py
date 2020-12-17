import time
import os
import loompy
import numpy as np
import pandas as pd
import hotspot
import hotspot.modules

loom_file = snakemake.input['loom']
latent_file = snakemake.input['latent']
hs_results_file = snakemake.input['hs_results']
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

hs_results = pd.read_table(hs_results_file, index_col=0)

# Have to do this because data_slideseq makes it a numpy array
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)
num_umi = pd.Series(num_umi, index=barcodes)

# Align to latent space
counts = counts.loc[:, latent.index]
num_umi = num_umi[latent.index]

# Subsample cells

# Find a window around the middle of the space with N_CELLS in it
def find_window(latent, N_CELLS):
    X = latent.iloc[:, 0].values
    Y = latent.iloc[:, 1].values
    midX = np.median(X)
    midY = np.median(Y)

    scaleX = X.max() - X.min()
    scaleY = Y.max() - Y.min()

    widthXMin = 0
    widthXMax = scaleX
    widthYMin = 0
    widthYMax = scaleY

    for _ in range(100):

        testX = (widthXMax + widthXMin)/2
        testY = (widthYMax + widthYMin)/2

        inWin = (
            (X < midX + testX/2) & (X > midX - testX/2) &
            (Y < midY + testY/2) & (Y > midY - testY/2)
        )

        cellsInWin = inWin.sum()

        if cellsInWin > N_CELLS:
            widthXMax = testX
            widthYMax = testY
        elif cellsInWin < N_CELLS:
            widthXMin = testX
            widthYMin = testY
        else:
            break

        if abs(cellsInWin - N_CELLS) <= 1:
            break

    else:
        raise Exception("Bad window convergence")

    return inWin


cellsInWin = find_window(latent, N_CELLS)

counts = counts.iloc[:, cellsInWin]
latent = latent.iloc[cellsInWin, :]
num_umi = num_umi.iloc[cellsInWin]

# Subsample Genes
# Need to remove genes that are note expressed in the subset
valid_genes = (counts > 0).sum(axis=1) > 0
counts = counts.loc[valid_genes, :]
hs_results = hs_results.loc[counts.index & hs_results.index]

genes_sub = hs_results.sort_values('Z', ascending=False).index[0:N_GENES]
counts = counts.loc[genes_sub, :]

hs = hotspot.Hotspot(counts, model=model, latent=latent, umi_counts=num_umi)

hs.create_knn_graph(
    weighted_graph=False, n_neighbors=N_NEIGHBORS, neighborhood_factor=3
)


start = time.time()

lcz = hs.compute_local_correlations(genes_sub, jobs=N_JOBS)

modules, Z = hotspot.modules.compute_modules(
    lcz, min_gene_threshold=20,
    z_threshold=1.65, core_only=False
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

stop = time.time()

with open(out_file, 'w') as fout:
    fout.write("Genes\tCells\tElapsedSeconds\n")
    fout.write("{}\t{}\t{:.3f}\n".format(
        counts.shape[0],
        N_CELLS,
        stop-start,
    ))
