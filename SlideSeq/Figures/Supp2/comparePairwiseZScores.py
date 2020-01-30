from scipy.cluster.hierarchy import linkage, leaves_list
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'


# %% Load HS Results for each puck

puck_analysis_dirs = {
    '9': '../../Puck_180819_9/',
    '10': '../../Puck_180819_10/',
    '11': '../../Puck_180819_11/',
    '12': '../../Puck_180819_12/',
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
    file_name = "hotspot/hotspot_pairs_z.txt.gz"

    results = pd.read_table(
        os.path.join(analysis_dir, file_name),
        index_col=0
    )
    return results


hs_pairs = {k: load_hs_pairs(v) for k, v in puck_analysis_dirs.items()}

# %%  compute ordering

x2 = '12'

x1s = ['9', '10', '11']

fig, axs = plt.subplots(3, 1, figsize=(3.5, 8),
                        gridspec_kw=dict(
                            top=.95,
                            left=.2,
                            hspace=.4,
                        ))


for ax, x1 in zip(axs.ravel(), x1s):

    plt.sca(ax)

    pairs1 = hs_pairs[x1]
    pairs2 = hs_pairs[x2]

    common = pairs1.index & pairs2.index

    pairs1 = pairs1.loc[common, common]
    pairs2 = pairs2.loc[common, common]

    vals1 = pairs1.values.ravel()
    vals2 = pairs2.values.ravel()

    plt.plot(
        vals2, vals1,
        'o', ms=1, rasterized=True,
    )
    ax.set_axisbelow(True)

    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    
    plt.xlim(-40, 65)
    plt.ylim(-40, 65)

    plt.xticks([-25, 0, 25, 50])
    plt.yticks([-25, 0, 25, 50])

    plt.ylabel("Puck_180819_{}".format(x1))
    plt.xlabel("Puck_180819_{}".format(x2))

    plt.axhline(xmin=-40, xmax=65, color='black', lw=.5)
    plt.axvline(ymin=-40, ymax=65, color='black', lw=.5)


# plt.show()
plt.savefig('Pairwise_ZScores.svg', dpi=300)

