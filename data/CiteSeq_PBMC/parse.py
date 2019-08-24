import numpy as np
import pandas as pd
import feather

counts_file = snakemake.input['counts']
adt_raw_file = snakemake.input['adt']
adt_clr_file = snakemake.input['adt_clr']

counts_out = snakemake.output['counts']
meta_out = snakemake.output['meta']
scaled_out = snakemake.output['scaled']
adt_out = snakemake.output['adt']
adt_clr_out = snakemake.output['adt_clr']

counts = pd.read_csv(counts_file, index_col=0)

# Mouse vs. Human counts?
counts['species'] = [
    'mouse' if 'mouse' in x.lower() else 'human' for x in counts.index]

zz = counts.groupby('species').sum()
ratio = np.log2(zz.loc['human'] / zz.loc['mouse'])

human_cell = (ratio > 0)

counts = counts.loc[counts.species == 'human'] \
    .drop('species', axis=1) \
    .loc[:, human_cell]

num_umi = counts.sum(axis=0)

# This dataset is already filtered by umi

adts = pd.read_csv(adt_raw_file, index_col=0)
adts_clr = pd.read_csv(adt_clr_file, index_col=0)

adts = adts.T.loc[counts.columns]
adts_clr = adts_clr.T.loc[counts.columns]

num_adt = adts.sum(axis=1)

meta = pd.DataFrame({
    'umi_rna': num_umi,
    'umi_adt': num_adt,
})

meta.to_csv(meta_out, sep="\t", compression="gzip")


cell_counts = (counts > 0).sum(axis=1)
counts = counts.loc[cell_counts >= 5, :]

scaled = counts.divide(num_umi, axis=1)

counts['index'] = counts.index
feather.write_dataframe(counts, counts_out)

scaled['index'] = scaled.index
feather.write_dataframe(scaled, scaled_out)


adts.to_csv(adt_out, sep="\t", compression="gzip")
adts_clr.to_csv(adt_clr_out, sep="\t", compression="gzip")
