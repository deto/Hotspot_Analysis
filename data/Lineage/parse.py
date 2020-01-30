import loompy
import h5py
import os
from scipy.sparse import csc_matrix
import numpy as np
import pandas as pd

h5_file = snakemake.input['dge']
cm_file = snakemake.input['characterMatrix']

out_file = snakemake.output['loom']
out_cm_file = snakemake.output['characterMatrix']


os.makedirs(os.path.dirname(out_file), exist_ok=True)

# %% Load in the counts and character matrix

cm = pd.read_table(cm_file, index_col=0)

fin = h5py.File(h5_file)
cells = fin['mm10']['barcodes'][:].astype('str')
data = fin['mm10']['data'][:]
indices = fin['mm10']['indices'][:]
indptr = fin['mm10']['indptr'][:]
shape = fin['mm10']['shape'][:]
genes = fin['mm10']['genes'][:].astype('str')
gene_names = fin['mm10']['gene_names'][:].astype('str')
fin.close()

counts_sp = csc_matrix((data, indices, indptr), shape=shape)
genes = pd.DataFrame({
    'GeneSymbol': gene_names,
    'EnsID': genes
}).set_index('EnsID')

valid_cells = {x for x in cells if x in cm.index}

is_valid = [x in valid_cells for x in cells]
is_valid_i = [i for i, x in enumerate(is_valid) if x]

cells = cells[is_valid]
counts_sp = counts_sp[:, is_valid_i]

counts = pd.DataFrame(
    counts_sp.todense(),
    index=genes.index, columns=cells
)

# Re-order counts to match cm
counts = counts.loc[:, cm.index]

num_umi = counts.sum(axis=0)

# Gene filtering
n_detects = (counts > 0).sum(axis=1)
gene_passes = n_detects >= 10
counts = counts.loc[gene_passes]
genes = genes.loc[gene_passes]

scaled = counts \
    .divide(num_umi, axis=1) \
    * 10000

# Save results

row_attrs = {
    "Symbol": genes.GeneSymbol.values.astype('str'),
    "EnsID": genes.index.values.astype('str')
}

col_attrs = {
    "Barcode": counts.columns.values.astype('str'),
    "NumUmi": num_umi.values,
}

layers = {
    '': counts.values,
    'scaled': scaled.values,
}

loompy.create(out_file, layers, row_attrs, col_attrs)
cm.to_csv(out_cm_file, sep="\t")
