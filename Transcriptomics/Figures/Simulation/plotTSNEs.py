import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

plt.rcParams['svg.fonttype'] = 'none'

# %%

tsne = pd.read_table("../../Simulation4/rep2/tsne/tsne_pca.txt", index_col=0)

ce = pd.read_table(
    "../../../data/Simulated4/rep2/cell_effects.txt", index_col=0
)

data = tsne.join(ce)

# %%


fig, axs = plt.subplots(1, ce.shape[1]-1, figsize=(9, 2.5))
cbar_ax = fig.add_axes([.9, .1, .01, .2])

import matplotlib.colors
cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
    'grays', ['#D2D2D2', '#000000']
)

for ax, col in zip(
    axs.ravel(), ce.columns[1:]
):

    plt.sca(ax)
    vals = data[col]
    if 'evf1' in col or 'evf5' in col:
        vals = vals.subtract(vals.mean()).divide(vals.std())
        vmin = -1
        vmax = 1
    else:
        vmin = 0
        vmax = 1.5

    sc = plt.scatter(
        x=data.tsne1, y=data.tsne2, c=vals,
        s=1, vmin=vmin, vmax=vmax, cmap=cmap,
        rasterized=True,
    )
    plt.xticks([])
    plt.yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    plt.title(col[-1])

cb = plt.colorbar(sc, cax=cbar_ax, ticks=[-1, 1])
cb.set_label('Cell-Effect', labelpad=10, rotation=0, size=9, verticalalignment='center')
cbar_ax.set_yticklabels(['Low', 'High'], size=7)

plt.subplots_adjust(right=0.95, left=0.05, top=0.75)
plt.suptitle('Simulated Components')
# plt.show()
plt.savefig('Simulation_tSNES.svg', dpi=300)
