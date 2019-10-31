import loompy
import numpy as np
import pandas as pd
import os


in_loom = snakemake.input['loom']
in_ab = snakemake.input['ab']

out_loom = snakemake.output['loom']
out_ab = snakemake.output['ab']

out_dir = os.path.dirname(out_loom)
os.makedirs(out_dir, exist_ok=True)

# in_file = "raw/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
# out_file = "data.loom"

ab = pd.read_csv(in_ab, index_col=0, sep="\t")

# Further filter by MitoPercent and NumUmi
with loompy.connect(in_loom, 'r') as ds:
    num_umi = ds.ca['NumUmi'][:]
    mito_percent = ds.ca['MitoPercent'][:]

is_healthy = (
    (num_umi > 3000) &
    (mito_percent < 16)
)

ab_healthy = ab.loc[is_healthy]

with loompy.connect(in_loom, mode='r') as ds:
    with loompy.new(out_loom) as ds_out:
        view = ds.view[:, is_healthy]
        ds_out.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)

ab_healthy.to_csv(out_ab, sep="\t", compression='gzip')
