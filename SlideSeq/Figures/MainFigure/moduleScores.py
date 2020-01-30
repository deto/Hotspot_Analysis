import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json
import loompy
plt.rcParams['svg.fonttype'] = 'none'

# %% Load things

results_hs = pd.read_table(
    "../../Puck_180819_12/hotspot/hotspot.txt",
    index_col=0
)

positions = pd.read_table(
    "../../Puck_180819_12/positions/positions.txt", index_col=0
)

scores = pd.read_table(
    "../../Puck_180819_12/hotspot/module_scores.txt.gz", index_col=0
)

modules = pd.read_table(
    "../../Puck_180819_12/hotspot/modules.txt", index_col=0
).Cluster

positions = positions.loc[scores.index]

name_map_r = json.load(open("module_name_map.json"))
name_map = {v: k for k, v in name_map_r.items()}

loom_file = "../../../data/SlideSeq/Puck_180819_12/data.loom"
with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    scaled = ds.layers['scaled'][:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

# Have to do this because data_slideseq makes it a numpy array
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
scaled = pd.DataFrame(scaled, index=gene_info.index, columns=barcodes)

# %% plot things

all_modules = [name_map[int(x)] for x in scores.columns]
all_modules = sorted(all_modules)

fig, axs = plt.subplots(2, len(all_modules), figsize=(15, 5))

from matplotlib.colors import LinearSegmentedColormap
gene_cmap = LinearSegmentedColormap.from_list(
    "grays", ['#dddddd', '#000000']
)

for i, m in enumerate(all_modules):

    ax = axs[0, i]
    plt.sca(ax)
    mod = name_map_r[m]

    vals = scores[str(mod)]
    vmin = np.percentile(vals, 5)
    vmax = np.percentile(vals, 95)
    plt.scatter(
        x=positions.Comp1,
        y=positions.Comp2,
        c=vals, vmin=vmin, vmax=vmax,
        s=.1, alpha=0.5, rasterized=True,
    )
    for sp in ax.spines.values():
        sp.set_visible(False)

    plt.xticks([])
    plt.yticks([])
    plt.title(m)

    # Plot top gene
    ax = axs[1, i]
    plt.sca(ax)

    mod_genes = modules.index[modules == mod]
    top_gene = results_hs.loc[mod_genes].sort_values('Z').index[-1]
    vals = np.log2(scaled.loc[top_gene]+1)

    #vmin = np.percentile(vals, 5)
    #vmax = np.percentile(vals, 95)
    #if vmin == vmax:
    vmin = 0
    vmax = np.percentile(vals, 99)

    plt.scatter(
        x=positions.Comp1,
        y=positions.Comp2,
        c=vals, vmin=vmin, vmax=vmax, cmap=gene_cmap,
        s=.2, alpha=0.3, rasterized=True,
    )
    for sp in ax.spines.values():
        sp.set_visible(False)

    plt.xticks([])
    plt.yticks([])
    plt.title(top_gene)

# colorbars

from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

cbax = fig.add_axes([.93, .6, .01, .2])
plt.colorbar(
    ScalarMappable(norm=Normalize(), cmap="viridis"),
    cax=cbax,
    ticks=[-1, 0, 1],
    label="Module Score"
)

cbax = fig.add_axes([.93, .2, .01, .2])
plt.colorbar(
    ScalarMappable(norm=Normalize(), cmap=gene_cmap),
    cax=cbax,
    ticks=[-1, 0, 1],
    label="Normalized Log\nExpression",
)

plt.show()
# plt.savefig("Module_Scores.svg", dpi=300)
