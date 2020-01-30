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
    "../../CD4_w_protein/hotspot/hotspot_pairs_z.txt.gz",
    index_col=0
)

hs_results = pd.read_table(
    "../../CD4_w_protein/hotspot/hotspot_hvg.txt",
    index_col=0
)

ens_map = {x: y for x, y in zip(hs_results.index, hs_results.Symbol)}


# %% Load Modules

modules = pd.read_table(
    "../../CD4_w_protein/hotspot/modules.txt",
    index_col=0
).Cluster

Z = pd.read_table(
    "../../CD4_w_protein/hotspot/linkage.txt",
    header=None
).values

scores = pd.read_table(
    "../../CD4_w_protein/hotspot/module_scores.txt.gz",
    index_col=0
)

proj = pd.read_table(
    "../../CD4_w_protein/umap/umap_hvg.txt",
    index_col=0
)

# %% Plot Modules

colors = list(plt.get_cmap("tab20").colors)
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
plt.savefig('CD4_ModuleCorelation.svg', dpi=300)

# %% Plot a color legend
import matplotlib.patches
from natsort import natsorted

handles = []
labels = []
for i in natsorted(modules.unique()):
    if i == -1: continue;

    patch = matplotlib.patches.Patch(
        facecolor=module_colors[i]
    )
    handles.append(patch)
    labels.append(str(i))

plt.figure()
plt.legend(handles, labels)
#plt.show()
plt.savefig('CD4_ModuleCorelation_legend.svg')
    

# %% Plot Module Scores
# %% Load data

fig, axs = plt.subplots(
    3, 4, figsize=(9, 7),
    gridspec_kw=dict(
        left=.05, right=.95, top=.95, bottom=.05
    ),
)

cbar_ax = fig.add_axes(
    [.95, .05, .01, .1]
)

for ax, module in zip(axs.ravel(), scores.columns):
    plt.sca(ax)

    vals = scores[module]

    vmin = np.percentile(vals, 5)
    vmax = np.percentile(vals, 95)
    vmin = -2
    vmax = 2

    sc = plt.scatter(
        x=proj.iloc[:, 0],
        y=proj.iloc[:, 1],
        c=vals,
        vmin=vmin, vmax=vmax,
        s=2, rasterized=True
    )
    plt.xticks([])
    plt.yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    plt.title(module)

plt.colorbar(sc, cax=cbar_ax)

# plt.show()
plt.savefig('CD4_ModuleScore_UMAPs.svg', dpi=300)
