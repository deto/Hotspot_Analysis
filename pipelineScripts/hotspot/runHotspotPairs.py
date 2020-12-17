import loompy
import pandas as pd
import hotspot

loom_file = snakemake.input['loom']
latent_file = snakemake.input['latent']
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

try:
    genes_file = snakemake.input['genes']
except AttributeError:
    genes_file = None

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

if use_umi:
    num_umi = pd.Series(num_umi, index=barcodes)
else:
    num_umi = pd.Series(1.0, index=barcodes)

# Align to latent space
counts = counts.loc[:, latent.index]
num_umi = num_umi[latent.index]

# need counts, latent, and num_umi

hs = hotspot.Hotspot(counts, model=model, latent=latent, umi_counts=num_umi)

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

# Overrides the hs results if desired
if genes_file is not None:
    hs_genes = pd.Index(pd.read_table(genes_file, header=None).iloc[:, 0].tolist())

hs_genes = hs_genes & counts.index

lcz = hs.compute_local_correlations(hs_genes, jobs=20)
lcz = lcz.fillna(0)

lcz.to_csv(out_file_lcz, sep="\t", compression="gzip")
