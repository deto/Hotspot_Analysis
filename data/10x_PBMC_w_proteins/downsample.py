import numpy as np
import pandas as pd
import h5py
import loompy
from tqdm import tqdm
import scipy.sparse as sparse

mol_info_file = snakemake.input['mol_info']
loom_file_in = snakemake.input['loom']

out_file = snakemake.output['loom']
out_file_ab = snakemake.output['ab']

downsample_rate = int(snakemake.params['rate'])/100

def load_mol_info(mol_info):
    hf = h5py.File(mol_info, "r")

    # barcode groups
    umis = pd.DataFrame({
        "barcode_idx": hf["barcode_idx"][:],
        "count": hf["count"][:],
        "feature_idx": hf["feature_idx"][:],
        "gem_group": hf["gem_group"][:],
        "library_idx": hf["library_idx"][:],
        "umi": hf["umi"][:],
    })

    # Feature group
    features = pd.DataFrame({
        "feature_type": hf["features/feature_type"][:].astype('str'),
        "genome": hf["features/genome"][:].astype('str'),
        "id": hf["features/id"][:].astype('str'),
        "name": hf["features/name"][:].astype('str'),
        "pattern": hf["features/pattern"][:].astype('str'),
        "read": hf["features/read"][:].astype('str'),
        "sequence": hf["features/sequence"][:].astype('str'),
    })

    barcodes = hf["barcodes"][:].astype('str')

    hf.close()

    return umis, features, barcodes


def umis_to_feature_mat(umis, features, barcodes, threshold=1):
    umis = umis.loc[umis['count'] >= threshold]
    feature_mat = umis \
        .groupby(['barcode', 'feature_idx']).size() \
        .reset_index() \
        .pivot(columns='barcode', index='feature_idx', values=0) \
        .fillna(0)

    missing_features = pd.Index(
        np.arange(features.shape[0])
    ).difference(feature_mat.index)

    missing_feature_mat = pd.DataFrame(
        0.0, index=missing_features, columns=feature_mat.columns
    )
    missing_feature_mat.index.name = 'feature_idx'

    feature_mat = pd.concat((feature_mat, missing_feature_mat), axis=0)

    feature_mat.index = features['id'][feature_mat.index]
    feature_mat = feature_mat.loc[features.id]
    return feature_mat


umis, features, barcodes = load_mol_info(mol_info_file)
umis['barcode'] = [
    "{}-{}".format(barcodes[bb], gg)
    for bb, gg in tqdm(zip(umis.barcode_idx, umis.gem_group))
]


# Load the barcode list for CD4 cells from the loom file
with loompy.connect(loom_file_in, 'r') as ds:
    cd4_barcodes = set(ds.col_attrs['Barcode'])

# Subset to only the cd4 barcodes
is_cd4_umi = umis['barcode'].isin(cd4_barcodes)

umis = umis.loc[is_cd4_umi].copy()

# Now, downsample the reads
from scipy.stats import binom
umis['count'] = binom.rvs(n=umis['count'].values, p=downsample_rate)

umis = umis.loc[umis['count'] > 0]

feature_mat = umis_to_feature_mat(umis, features, barcodes, threshold=1)

# %% Now, common code to create the new loom file

matrix = sparse.csc_matrix(feature_mat.values)

# Split the matrix into the mRNA and the antibody counts
mRNA_rows = features.feature_type == "Gene Expression"
i_mRNA_rows = mRNA_rows.values.nonzero()[0]
protein_rows = features.feature_type == "Antibody Capture"
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
    proteins.toarray(), index=features_proteins.id.values.astype('str'),
    columns=feature_mat.columns,
).T

# Other Stats

mito_genes = [x.lower().startswith('mt-') for x in features_genes.name.values.astype('str')]
mito_genes = np.array(mito_genes)

mito_count = np.array(mRNA[mito_genes].sum(axis=0)).ravel()
mito_percent = mito_count / num_umi * 100


# Save results

row_attrs = {
    "Symbol": features_genes.name.values.astype('str'),
    "EnsID": features_genes.id.values.astype('str'),
}

col_attrs = {
    "Barcode": feature_mat.columns.values.astype('str'),
    "NumUmi": num_umi,
    "NumAb": num_ab,
    "MitoPercent": mito_percent,
}

layers = {
    '': mRNA,
    'scaled': scaled
}

loompy.create(out_file, layers, row_attrs, col_attrs)

# Not going to store this in the loom file because I'd have to make them part
# of all the layers.
proteins_df.to_csv(out_file_ab, sep="\t", compression='gzip')
