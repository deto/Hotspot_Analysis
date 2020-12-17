import numpy as np
import pandas as pd
from sklearn.decomposition import FastICA
import loompy


loom_file = snakemake.input['loom']
genes = snakemake.input['genes']
n_components = snakemake.params['n_components']
ica_out = snakemake.output['components']

try:
    ica_cell_loadings = snakemake.output['cell_loadings']
except AttributeError:
    ica_cell_loadings = None

with loompy.connect(loom_file, 'r') as ds:
    expression = ds.layers['scaled'][:]
    ens_ids = ds.ra['EnsID'][:]
    barcodes = ds.ca['Barcode'][:]

expression = pd.DataFrame(expression, index=ens_ids, columns=barcodes)
expression = np.log2(expression + 1)

expression_std = expression \
    .subtract(expression.mean(axis=1), axis=0) \
    .divide(expression.std(axis=1), axis=0)

# Gene filtering
genes = pd.read_table(genes, header=None).iloc[:, 0].tolist()
expression_std = expression_std.loc[genes]

model = FastICA(n_components=n_components, random_state=0, max_iter=2000)

# Old way

# model.fit(expression_std.values.T)
# gene_components = model.components_.T

# New way

source = model.fit_transform(expression_std.values)
gene_components = source

# Save outputs
gene_components = pd.DataFrame(
    gene_components, index=expression_std.index,
    columns=['ICA{}'.format(i+1) for i in range(gene_components.shape[1])]
)


gene_components.to_csv(ica_out, sep="\t")

if ica_cell_loadings is not None:
    cell_loadings = model.components_.T
    cell_loadings = pd.DataFrame(
        cell_loadings, index=expression_std.columns, columns=gene_components.columns
    )
    cell_loadings.to_csv(ica_cell_loadings, sep="\t")
