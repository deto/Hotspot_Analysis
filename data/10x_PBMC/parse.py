import h5py
import loompy
import scipy.sparse as sparse
import numpy as np


in_file = snakemake.input['h5']
out_file = snakemake.output['loom']

# in_file = "raw/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
# out_file = "data.loom"

f = h5py.File(in_file, 'r')

barcodes = f['matrix/barcodes'][:]

names = f['matrix/features/name'][:]
ens_ids = f['matrix/features/id'][:]

data = f['matrix/data'][:]
indices = f['matrix/indices'][:]
indptr = f['matrix/indptr'][:]
shape = f['matrix/shape'][:]


matrix = sparse.csc_matrix(
    (data, indices, indptr),
    shape=shape
)


num_umi = matrix.sum(axis=0)
num_umi = np.array(num_umi).ravel()

# Gene filtering
num_umi_mat = sparse.lil_matrix((len(num_umi), len(num_umi)))
num_umi_mat.setdiag(1/num_umi)

scaled = matrix.dot(num_umi_mat) * np.median(num_umi)

gene_counts = np.array((matrix > 0).sum(axis=1)).ravel()
valid_genes = gene_counts > 10

names = names[valid_genes]
ens_ids = ens_ids[valid_genes]

matrix = matrix[valid_genes, :]
scaled = scaled[valid_genes, :]

# Save results

row_attrs = {
    "Symbol": names,
    "EnsID": ens_ids
}

col_attrs = {
    "Barcode": barcodes,
    "NumUmi": num_umi
}

layers = {
    '': matrix,
    'scaled': scaled
}

loompy.create(out_file, layers, row_attrs, col_attrs)
