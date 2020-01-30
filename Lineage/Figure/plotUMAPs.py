import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


plt.rcParams['svg.fonttype'] = 'none'

# %% Load data

tsne = pd.read_table("../Embryo3/umap/umap.txt", index_col=0)

scores = pd.read_table(
    "../Embryo3/hotspot/module_scores_lineage_tree.txt.gz",
    index_col=0
)
scores = scores.loc[tsne.index]

scores = scores.divide(scores.std(axis=0), axis=1)

# %% Plot module scores

fig, axs = plt.subplots(1, 5, figsize=(12, 3))
cax = fig.add_axes(
    [.92, .1, .007, .2]
)

for ax, mod in zip(axs.ravel(), scores.columns):
    sc = scores[mod]
    # vmin = np.percentile(sc, 5)
    # vmax = np.percentile(sc, 95)
    vmin = -1
    vmax = 1
    plt.sca(ax)
    sc = plt.scatter(
        tsne.iloc[:, 0], tsne.iloc[:, 1],
        s=3, c=sc, vmin=vmin, vmax=vmax,
        rasterized=True)
    plt.xticks([])
    plt.yticks([])
    plt.title(mod)
    for sp in ax.spines.values():
        sp.set_visible(False)

plt.colorbar(sc, cax=cax, ticks=[vmin, vmax])
plt.subplots_adjust(left=0.02, right=0.9, top=0.8)
# plt.show()
plt.savefig('moduleUMAPS.svg', dpi=300)
