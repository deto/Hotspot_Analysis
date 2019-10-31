import numpy as np
import pandas as pd
import loompy

counts_file = snakemake.input['counts']
adt_raw_file = snakemake.input['adt']

out_file = snakemake.output['loom']
out_file_ab = snakemake.output['ab']

counts = pd.read_csv(counts_file, index_col=0)

# Mouse vs. Human counts?
# Remove mouse cells and remove mouse genes

counts['species'] = [
    'mouse' if 'mouse' in x.lower() else 'human' for x in counts.index]

zz = counts.groupby('species').sum()
ratio = np.log2(zz.loc['human'] / zz.loc['mouse'])

human_cell = (ratio > 0)  # Peaks are easily separable

counts = counts.loc[counts.species == 'human'] \
    .drop('species', axis=1) \
    .loc[:, human_cell]

num_umi = counts.sum(axis=0)

# Filter out genes barely expressed

cell_counts = (counts > 0).sum(axis=1)
counts = counts.loc[cell_counts >= 5, :]

# Rename genes - remove _HUMAN

counts.index = [x.replace("HUMAN_", "") for x in counts.index]

# Get MitoPercent
mito_genes = [x for x in counts.index if 'mt-' in x.lower()]
mito_count = counts.loc[mito_genes].sum(axis=0)
mito_percent = mito_count / num_umi * 100

# This dataset is already filtered by umi

# Load ADTs

adts = pd.read_csv(adt_raw_file, index_col=0)
adts = adts.T.loc[counts.columns]

num_adt = adts.sum(axis=1)

# Scale expression data

scaled = counts.divide(num_umi, axis=1) * 10000


# Save results

row_attrs = {
    "Symbol": counts.index.values.astype('str'),
    "EnsID": counts.index.values.astype('str'),
}

col_attrs = {
    "Barcode": counts.columns.values.astype('str'),
    "NumUmi": num_umi.values,
    "NumAb": num_adt.values,
    "MitoPercent": mito_percent.values,
}

layers = {
    '': counts.values,
    'scaled': scaled.values,
}

loompy.create(out_file, layers, row_attrs, col_attrs)

# Not going to store this in the loom file because I'd have to make them part
# of all the layers.
adts.to_csv(out_file_ab, sep="\t", compression='gzip')
