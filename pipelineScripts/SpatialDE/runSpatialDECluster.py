import pandas as pd
import loompy

import NaiveDE
import SpatialDE


# Load some data

print('Loading Data...')
latent_file = snakemake.input['latent']
latent = pd.read_csv(latent_file, index_col=0, sep="\t")

loom_file = snakemake.input['loom']
sde_results_file = snakemake.input['sde_results']

out_file_hist = snakemake.output['hist']
out_file_patterns = snakemake.output['patterns']

with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')

sde_results = pd.read_table(sde_results_file)

print('Filtering Data...')
# Do some filtering
valid_genes = (counts > 0).sum(axis=1) > 50
counts = counts[valid_genes, :]
gene_info = gene_info.loc[valid_genes]

counts = counts.T

counts = pd.DataFrame(counts, columns=gene_info.index)

# prepare 'sample_info' like they want it
sample_info = latent.copy()
sample_info['total_counts'] = num_umi

# Some normalization


print('Normalize Data...')
norm_expr = NaiveDE.stabilize(counts.T).T

resid_expr = NaiveDE.regress_out(
    sample_info, norm_expr.T, 'np.log(total_counts)').T


X = sample_info[['Comp1', 'Comp2']].values

# Pick results subset
sde_results_sub = sde_results.sort_values('LLR').tail(500)

# Pick l=350 as average is around 300
# Pick C=5 as there seems to be at least 5 components

print('Running AES...')
histology_results, patterns = SpatialDE.aeh.spatial_patterns(
    X, resid_expr, sde_results_sub, C=5, l=350, verbosity=999)

print('Saving Results...')
histology_results.to_csv(out_file_hist, sep="\t")
patterns.to_csv(out_file_patterns, sep="\t")
print('Done!')
