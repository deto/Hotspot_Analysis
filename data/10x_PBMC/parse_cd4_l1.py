import loompy
import numpy as np
import pandas as pd
import os


in_loom = snakemake.input['loom']
in_clusters = snakemake.input['clusters']

out_loom = snakemake.output['loom']

selected_clusters = snakemake.params['selected_clusters']

out_dir = os.path.dirname(out_loom)
os.makedirs(out_dir, exist_ok=True)

# in_file = "raw/pbmc_10k_v3_filtered_feature_bc_matrix.h5"
# out_file = "data.loom"

clusters = pd.read_table(in_clusters, index_col=0)
valid_bcs = clusters.index[[c in selected_clusters for c in clusters.Cluster]]
valid_bcs = set(valid_bcs)

with loompy.connect(in_loom, mode='r') as ds:
    barcodes = ds.ca['Barcode'][:]

    cell_passes = np.array([b in valid_bcs for b in barcodes])

    with loompy.new(out_loom) as ds_out:
        view = ds.view[:, cell_passes]
        ds_out.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)
