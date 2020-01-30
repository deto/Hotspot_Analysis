import numpy as np
import pandas as pd
import hotspot.modules
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable

plt.rcParams['svg.fonttype'] = 'none'

# %% Load data

z_scores = pd.read_table(
    "../Embryo3/hotspot/hotspot_pairs_lineage_tree_z.txt.gz",
    index_col=0
)

hs_results = pd.read_table(
    "../Embryo3/hotspot/hotspot_lineage_tree.txt",
    index_col=0
)

ens_map = {x: y for x, y in zip(hs_results.index, hs_results.Symbol)}


# %% Load Modules

modules = pd.read_table(
    "../Embryo3/hotspot/modules_lineage_tree.txt",
    index_col=0
).Cluster

Z = pd.read_table(
    "../Embryo3/hotspot/linkage_lineage_tree.txt",
    header=None
).values

# %% Plot Modules

colors = list(plt.get_cmap("tab10").colors)
module_colors = {i: colors[(i-1) % len(colors)] for i in modules.unique()}
module_colors[-1] = '#ffffff'

cm = ScalarMappable(norm=Normalize(0, 0.05, clip=True), cmap="viridis")
row_colors1 = pd.Series(
    [module_colors[i] for i in modules],
    index=z_scores.index,
)

row_colors = pd.DataFrame({
    "Modules": row_colors1,
})

zvals = z_scores.values.ravel()
vmax = 8
vmin = -8

cm = sns.clustermap(
    z_scores,
    row_linkage=Z,
    col_linkage=Z,
    vmin=vmin,
    vmax=vmax,
    cmap="RdBu_r",
    xticklabels=False,
    #yticklabels=False,
    yticklabels=[ens_map[x] for x in z_scores.index],
    row_colors=row_colors,
    rasterized=True,
)

fig = plt.gcf()
fig.patch.set_visible(False)
plt.sca(cm.ax_heatmap)
plt.ylabel("")
plt.xlabel("")
plt.yticks([])
# plt.show()
plt.savefig('Lineage_ModuleCorelation.svg', dpi=300)

# %% Plot a color legend
import matplotlib.patches
from natsort import natsorted

handles = []
labels = []
for i in natsorted(modules.unique()):

    if i == -1: continue

    patch = matplotlib.patches.Patch(
        facecolor=module_colors[i]
    )
    handles.append(patch)
    labels.append(str(i))

fig = plt.figure()
plt.legend(handles, labels)
fig.patch.set_visible(False)
# plt.show()
plt.savefig('Lineage_ModuleCorelation_Legend.svg')
