import h5py
import loompy
import scipy.sparse as sparse
import numpy as np
import pandas as pd


in_file = snakemake.input['h5']
out_file = snakemake.output['loom']
out_file_ab = snakemake.output['ab']

# in_file = "raw/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
# out_file = "data.loom"

f = h5py.File(in_file, 'r')

barcodes = f['matrix/barcodes'][:]

names = f['matrix/features/name'][:]
ens_ids = f['matrix/features/id'][:]
f_type = f['matrix/features/feature_type'][:]
genome = f['matrix/features/genome'][:]
pattern = f['matrix/features/id'][:]
read = f['matrix/features/read'][:]
sequence = f['matrix/features/sequence'][:]

features = pd.DataFrame({
    'Name': names,
    'Id': ens_ids,
    'Type': f_type,
    'Genome': genome,
    'Pattern': pattern,
    'Read': read,
    'Sequence': sequence,
})

data = f['matrix/data'][:]
indices = f['matrix/indices'][:]
indptr = f['matrix/indptr'][:]
shape = f['matrix/shape'][:]


matrix = sparse.csc_matrix(
    (data, indices, indptr),
    shape=shape
)

# Split the matrix into the mRNA and the antibody counts
mRNA_rows = features.Type == b"Gene Expression"
i_mRNA_rows = mRNA_rows.values.nonzero()[0]
protein_rows = features.Type == b"Antibody Capture"
i_protein_rows = protein_rows.values.nonzero()[0]

proteins = matrix[i_protein_rows, :]
mRNA = matrix[i_mRNA_rows, :]

del matrix

features_genes = features.loc[mRNA_rows]
features_proteins = features.loc[protein_rows]

del features

# Gene filtering
num_umi = mRNA.sum(axis=0)
num_umi = np.array(num_umi).ravel()

num_umi_mat = sparse.lil_matrix((len(num_umi), len(num_umi)))
num_umi_mat.setdiag(1/num_umi)

scaled = mRNA.dot(num_umi_mat) * np.median(num_umi)

gene_counts = np.array((mRNA > 0).sum(axis=1)).ravel()
valid_genes = gene_counts > 10

features_genes = features_genes.loc[valid_genes]

mRNA = mRNA[valid_genes, :]
scaled = scaled[valid_genes, :]

# Antibody stats
num_ab = proteins.sum(axis=0)
num_ab = np.array(num_ab).ravel()

proteins_df = pd.DataFrame(
    proteins.toarray(), index=features_proteins.Id.values.astype('str'),
    columns=barcodes.astype('str')
).T

# Save results

row_attrs = {
    "Symbol": features_genes.Name.values.astype('str'),
    "EnsID": features_genes.Id.values.astype('str'),
}

col_attrs = {
    "Barcode": barcodes,
    "NumUmi": num_umi,
    "NumAb": num_ab,
}

layers = {
    '': mRNA,
    'scaled': scaled
}

loompy.create(out_file, layers, row_attrs, col_attrs)

# Not going to store this in the loom file because I'd have to make them part
# of all the layers.
proteins_df.to_csv(out_file_ab, sep="\t", compression='gzip')
