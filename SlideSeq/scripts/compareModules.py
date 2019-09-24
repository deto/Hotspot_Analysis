from scipy.cluster.hierarchy import linkage, leaves_list
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns


# %% Load HS Results for each puck

puck_analysis_dirs = {
    '9': '../Puck_180819_9/',
    '10': '../Puck_180819_10/',
    '11': '../Puck_180819_11/',
    '12': '../Puck_180819_12/',
}
puck_data_dirs = {
    '9': '../../data/SlideSeq/Puck_180819_9',
    '10': '../../data/SlideSeq/Puck_180819_10',
    '11': '../../data/SlideSeq/Puck_180819_11',
    '12': '../../data/SlideSeq/Puck_180819_12',
}


def load_hs_genes(analysis_dir):
    results = pd.read_table(
        os.path.join(analysis_dir, "hotspot/hotspot.txt"),
        index_col=0
    )
    return results


def load_hs_pairs(analysis_dir):
    if "_12" in analysis_dir:
        file_name = "hotspot/hotspot_pairs_z.txt.gz"
    else:
        file_name = "hotspot/hotspot_pairs_z_12.txt.gz"

    results = pd.read_table(
        os.path.join(analysis_dir, file_name),
        index_col=0
    )
    return results


hs_pairs = {k: load_hs_pairs(v) for k, v in puck_analysis_dirs.items()}

# %%  compute ordering


z0 = hs_pairs['12']

ii = leaves_list(
    linkage(z0, method='average', metric='cosine')
)

ig = z0.index[ii]

# %% plot clustering

fig, axs = plt.subplots(2, 2, figsize=(9, 9))
cbar_ax = fig.add_axes([.93, .11, .01, .1])

for ax, sample in zip(axs.ravel(), hs_pairs):
    z = hs_pairs[sample]
    igx = ig & z.index

    plt.sca(ax)

    if sample == "12":
        cbar = True
    else:
        cbar = False

    sns.heatmap(z.loc[igx, igx], vmin=-5, vmax=5, cmap='RdBu_r',
                xticklabels=False, yticklabels=False, cbar_ax=cbar_ax,
                cbar=cbar, cbar_kws={'ticks': [-5, 0, 5], 'label': 'Z-score'},
                rasterized=True)
    plt.xlabel('')
    plt.ylabel('')
    plt.title('Puck_180819_{}'.format(sample))

plt.show()
#plt.savefig('moduleReproducibility.svg', dpi=300)
