import numpy as np
import pandas as pd
import h5py
import loompy
from tqdm import tqdm
import scipy.sparse as sparse

loom_file_in = snakemake.input['loom']

loom_test_out = snakemake.output['loom_test']
loom_train_out = snakemake.output['loom_train']


# Load the barcode list for cells from the loom file
ds = loompy.connect(loom_file_in, 'r')

test_ii = np.random.choice(ds.shape[1], size=ds.shape[1]//2, replace=False)

train_ii = np.setdiff1d(np.arange(ds.shape[1]), test_ii)

test_ii = np.sort(test_ii)
train_ii = np.sort(train_ii)

ds_test_out = loompy.new(loom_test_out)

for (ix, selection, view) in ds.scan(items=test_ii, axis=1):
    ds_test_out.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)

ds_train_out = loompy.new(loom_train_out)

for (ix, selection, view) in ds.scan(items=train_ii, axis=1):
    ds_train_out.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)

ds.close()
ds_test_out.close()
ds_train_out.close()
