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


def clr(x):

    x = np.log(x+1)
    x = x.subtract(x.mean(axis=1), axis=0)

    return x


ab_clr = clr(ab)

# import matplotlib.pyplot as plt
# plt.figure()
# plt.plot(ab_clr['CD3'], ab_clr['CD4'], 'o', ms=2)
# #plt.plot(ab_clr['CD14'], ab_clr['CD4'], 'o', ms=2)
# plt.show()
#
# plt.figure()
# plt.hist(np.log10(ab['CD4']+1), 30)
# plt.show()


is_cd4 = (
    (ab_clr['CD14'] < 2) &
    (ab_clr['CD4'] > 3) &
    (ab_clr['CD3'] > 2.5)
)


ab_cd4 = ab.loc[is_cd4]

with loompy.connect(in_loom, mode='r') as ds:
    with loompy.new(out_loom) as ds_out:
        view = ds.view[:, is_cd4.values]
        ds_out.add_columns(view.layers, col_attrs=view.ca, row_attrs=view.ra)

ab_cd4.to_csv(out_ab, sep="\t", compression='gzip')
