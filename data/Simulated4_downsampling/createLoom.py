import loompy
import pandas as pd

out_file = snakemake.output['loom']
in_file = snakemake.input['obs']

oc = pd.read_table(in_file, index_col=0)

scaled = oc.values
scaled = scaled / scaled.sum(axis=0, keepdims=True) * 10000

num_umi = oc.values.sum(axis=0)

row_attrs = {
    "Symbol": oc.index.values.astype('str'),
    "EnsID": oc.index.values.astype('str'),
}

col_attrs = {
    "Barcode": oc.columns.values.astype('str'),
    "NumUmi": num_umi,
}

layers = {
    '': oc.values,
    'scaled': scaled
}

loompy.create(out_file, layers, row_attrs, col_attrs)
