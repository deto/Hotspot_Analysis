import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
import json

import scanpy.api as sc

latent_file = snakemake.input['latent']
out_file = snakemake.output['out']

out_dir = os.path.dirname(out_file)
if not os.path.exists(out_dir):
    os.mkdir(out_dir)

umap_plot_file = os.path.join(out_dir, 'umap.png')

data = pd.read_table(latent_file, index_col=0)

adata = sc.AnnData(data.values,
                   obs={'smp_names': data.index.tolist()},
                   var={'var_names': data.columns.tolist()})


try:
    n_neighbors = snakemake.params['n_neighbors']
except AttributeError:
    n_neighbors = 10

print('Running UMAP dimensionality reduction with n_neighbors={}'.format(n_neighbors))

sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=0)
sc.tl.umap(adata)
umap = adata.obsm['X_umap']

umap = pd.DataFrame(umap, index=adata.obs_names, columns=['umap1', 'umap2'])

# Plot results
ms = 3
if umap.shape[0] > 10000:
    ms = 2
if umap.shape[0] > 30000:
    ms = 1

ff = plt.figure(figsize=(5, 5))
plt.style.use('default')
plt.plot(umap.umap1, umap.umap2, 'o', ms=ms*.5)
plt.savefig(umap_plot_file)

# Save umap to disk
umap.to_csv(out_file, sep='\t')

# 3 components?
try:
    out_3_file = snakemake.output['out3']

    sc.tl.umap(adata, n_components=3)
    umap3 = adata.obsm['X_umap']

    umap3 = pd.DataFrame(umap3, index=adata.obs_names,
                         columns=['umap1', 'umap2', 'umap3'])

    umap3.to_csv(out_3_file, sep='\t')
except AttributeError:
    pass
