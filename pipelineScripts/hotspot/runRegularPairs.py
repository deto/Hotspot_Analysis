import loompy
import pandas as pd
import numpy as np
from hotspot.local_stats_pairs import create_centered_counts

loom_file = snakemake.input['loom']
hs_results_file = snakemake.input['hs_results']

out_file_lc = snakemake.output['results_lc']

model = snakemake.params['model']

fdrThresh = snakemake.params['fdrThresh']

try:
    genes_file = snakemake.input['genes']
except AttributeError:
    genes_file = None

try:
    topN = int(snakemake.params['topN'])
except AttributeError:
    topN = None

try:
    highXMeanCutoff = float(snakemake.params['highXMeanCutoff'])
except AttributeError:
    highXMeanCutoff = None

with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

hs_results = pd.read_table(hs_results_file, index_col=0)

# Have to do this because data_slideseq makes it a numpy array
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)
num_umi = pd.Series(num_umi, index=barcodes)

if highXMeanCutoff is not None:

    scaled = counts.divide(counts.sum(axis=0), axis=1)*10000
    gene_means = scaled.mean(axis=1)
    valid_genes = gene_means.index[gene_means < highXMeanCutoff]
    hs_results = hs_results.loc[valid_genes & hs_results.index]


if topN is None:
    hs_genes = hs_results.index[hs_results.FDR < fdrThresh]
else:
    hs_genes = hs_results.sort_values('Z').tail(topN).index

# Overrides the hs results if desired
if genes_file is not None:
    hs_genes = pd.Index(pd.read_table(genes_file, header=None).iloc[:, 0].tolist())

hs_genes = hs_genes & counts.index

counts = counts.loc[hs_genes]

c_counts = create_centered_counts(
    counts.values, model, num_umi.values)

lc = np.corrcoef(c_counts)
lc = pd.DataFrame(lc, index=counts.index, columns=counts.index)

lc.to_csv(out_file_lc, sep="\t", compression="gzip")
