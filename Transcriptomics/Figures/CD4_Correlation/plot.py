import numpy as np
import pandas as pd
import hotspot.modules
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


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

# %% Compute Modules

modules, Z = hotspot.modules.compute_modules(
    z_scores, min_gene_threshold=10, core_only=True
)

# %% Plot Modules

colors = list(plt.get_cmap("tab10").colors)
cm = ScalarMappable(norm=Normalize(0, 0.05, clip=True), cmap="viridis")
row_colors1 = pd.Series(
    [colors[i % 10] if i != -1 else "#ffffff" for i in modules],
    index=z_scores.index,
)

row_colors = pd.DataFrame({
    "Modules": row_colors1,
})

zvals = z_scores.values.ravel()
vmax = np.percentile(zvals, 90)
vmin = -1 * vmax

cm = sns.clustermap(
    z_scores,
    row_linkage=Z,
    col_linkage=Z,
    vmin=vmin,
    vmax=vmax,
    cmap="RdBu_r",
    xticklabels=False,
    yticklabels=False,
    #yticklabels=[ens_map[x] for x in z_scores.index],
    row_colors=row_colors,
    rasterized=True,
)

fig = plt.gcf()
fig.patch.set_visible(False)
plt.sca(cm.ax_heatmap)
plt.ylabel("")
plt.xlabel("")
#plt.show()
plt.savefig('CD4_ModuleCorelation.svg', dpi=300)

# %% Plot a color legend
import matplotlib.patches

handles = []
labels = []
for i in modules.unique():
    if i == -1: continue;

    patch = matplotlib.patches.Patch(
        facecolor=colors[i % 10]
    )
    handles.append(patch)
    labels.append(str(i))

plt.figure()
plt.legend(handles, labels)
plt.show()
    

# %% Plot Module Scores
# %% Load data

proj_file = "../../CD4_w_protein/umap/umap_hvg.txt"
loom_file = "../../../data/10x_PBMC_w_proteins/cd4/data.loom"
latent_file = "../../CD4_w_protein/scvi/hvg/latent.txt.gz"
model = "danb"
n_neighbors = 30


import loompy
import hotspot
import hotspot.modules


proj = pd.read_table(proj_file, index_col=0)

with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

latent = pd.read_table(latent_file, index_col=0)

# Have to do this because data_slideseq makes it a numpy array
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)
num_umi = pd.Series(num_umi, index=barcodes)

# Align to latent space
counts = counts.loc[:, latent.index]
num_umi = num_umi[latent.index]
proj = proj.loc[latent.index]

# need counts, latent, and num_umi

hs = hotspot.Hotspot(counts, latent, num_umi)
hs.create_knn_graph(
    weighted_graph=False, n_neighbors=n_neighbors, neighborhood_factor=3
)

# %% Plot scores for all modules

fig, axs = plt.subplots(3, 3, figsize=(9, 9))

modules_to_plot = sorted([x for x in modules.unique() if x != -1])

# Get the scores
module_scores = {}
for module in modules_to_plot:
    module_genes = modules.index[modules == module]

    scores = hotspot.modules.compute_scores(
        counts.loc[module_genes].values, model, num_umi.values,
        hs.neighbors.values, hs.weights.values
    )

    module_scores[module] = scores


for ax, module in zip(axs.ravel(), modules_to_plot):
    plt.sca(ax)

    scores = module_scores[module]

    vmin = np.percentile(scores, 5)
    vmax = np.percentile(scores, 95)

    plt.scatter(
        x=proj.iloc[:, 0],
        y=proj.iloc[:, 1],
        c=scores,
        vmin=vmin, vmax=vmax,
        s=2, rasterized=True
    )
    plt.xticks([])
    plt.yticks([])
    for sp in ax.spines.values():
        sp.set_visible(False)
    plt.title(module)

plt.show()
plt.savefig('CD4_ModuleScore_UMAPs.svg', dpi=300)
