import loompy
import scipy.sparse as sparse
import numpy as np
import pandas as pd

dge_file = snakemake.input['dge']
pos_file = snakemake.input['pos']
out_file = snakemake.output['loom']

pos = pd.read_csv(pos_file, index_col=0)
dge = pd.read_csv(dge_file, index_col=0)

dge = dge.loc[:, pos.index]
barcodes = pos.index.values

matrix = sparse.csc_matrix(dge.values)


num_umi = matrix.sum(axis=0)
num_umi = np.array(num_umi).ravel()

# Gene filtering
num_umi_mat = sparse.lil_matrix((len(num_umi), len(num_umi)))
num_umi_mat.setdiag(1/num_umi)

scaled = matrix.dot(num_umi_mat) * np.median(num_umi)

gene_counts = np.array((matrix > 0).sum(axis=1)).ravel()
valid_genes = gene_counts >= 10

names = dge.index[valid_genes].values
ens_ids = names  # since we don't have ensids for these

matrix = matrix[valid_genes, :]
scaled = scaled[valid_genes, :]

# Fix positioning
# We swap x'=y and y'=-x to match the slides in the paper
pos2 = pd.DataFrame(
    {
        'X': pos.ycoord,
        'Y': pos.xcoord*-1,
    }, index=pos.index
)


# Save results

row_attrs = {
    "Symbol": names,
    "EnsID": ens_ids
}

col_attrs = {
    "Barcode": barcodes,
    "NumUmi": num_umi,
    "Position": pos2.values
}

layers = {
    '': matrix,
    'scaled': scaled
}

loompy.create(out_file, layers, row_attrs, col_attrs)
