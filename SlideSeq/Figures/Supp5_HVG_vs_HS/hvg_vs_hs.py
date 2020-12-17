import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import loompy
from bio_utils.plots import hover_plot
from matplotlib.colors import LinearSegmentedColormap
plt.rcParams['svg.fonttype'] = 'none'

# %% Load results

hs = pd.read_table("../../Puck_180819_12/hotspot/hotspot.txt", index_col=0)
hvg = pd.read_table("../../Puck_180819_12/genes/hvg_info.txt", index_col=0)
hvg = hvg.loc[hs.index]


# %% Plot one vs other

plot_data = hs[['Z']].join(hvg[['gene.dispersion.scaled']])

plt.figure(figsize=(5, 5))

plt.plot(
    plot_data['gene.dispersion.scaled'],
    plot_data['Z'],
    'o', ms=2, rasterized=True
)

plt.gca().set_axisbelow(True)
plt.xlabel('Scaled Dispersion')
plt.ylabel('Autocorrelation Z')
plt.grid(color='#BBBBBB', ls=(0, (5, 5)), lw=.5)
plt.subplots_adjust(left=0.15, right=1-0.15, bottom=0.15, top=1-0.15)

plt.ylim(-10, 260)  # Zoom in to see the points better
plt.xlim(-4, 8)
plt.savefig('autocorrelation_vs_hvg.svg', dpi=300)

# %%

# hover_plot(
#     plot_data['gene.dispersion.scaled'],
#     plot_data['Z'],
#     plot_data.index,
#     'o', ms=2
# )
# plt.show()

# %% Plot rank vs rank


# %%

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

positions = pd.read_table(
    "../../Puck_180819_12/positions/positions.txt", index_col=0
)

positions = positions.loc[scaled.columns]

# %%
gene_cmap = LinearSegmentedColormap.from_list(
    "grays", ['#dddddd', '#000000']
)


def plot_gene(gene, positions, scaled, ax=None):

    if ax is None:
        ax = plt.gca()
    else:
        plt.sca(ax)

    vals = np.log2(scaled.loc[gene]+1)

    vmin = 0
    vmax = np.percentile(vals[vals > 0], 80)

    plot_data = pd.DataFrame({
        'X': positions.Comp1,
        'Y': positions.Comp2,
        'Expression': vals,
    }).sort_values('Expression')

    plot_data['S'] = 0.5
    plot_data.loc[plot_data.Expression > 0, 'S'] = 1

    plt.scatter(
        x=plot_data.X,
        y=plot_data.Y,
        c=plot_data.Expression,
        s=plot_data.S,
        vmin=vmin, vmax=vmax, cmap=gene_cmap,
        alpha=0.7, rasterized=True, edgecolors='none',
    )
    for sp in ax.spines.values():
        sp.set_visible(False)

    plt.xticks([])
    plt.yticks([])
    plt.title(gene)


# %%

# top_genes_hs = hs.Z.sort_values(ascending=False).index[0:6]

# top_genes_hvg = hvg['gene.dispersion.scaled'].sort_values(ascending=False).index[0:6]
# top_genes_hvg = hvg.loc[hvg['gene.mean'] > 1] \
#     .sort_values('gene.dispersion.scaled', ascending=False).index[0:6]


top_genes_hs = ['Plp1', 'Enpp2', 'Cartpt', 'Pcp2', 'Car8', 'Calb1']
top_genes_hvg = ['Ube2q1', 'Rexo4', 'Letm1', 'Cnksr2', 'Xrn2', 'Prmt5']


fig, axs = plt.subplots(2, 3, figsize=(6.5, 4))


for gene, ax in zip(top_genes_hs, axs.ravel()):

    plot_gene(gene, positions, scaled, ax=ax)

# plt.show()
plt.savefig('top_genes_hs.svg', dpi=600)

# %%

fig, axs = plt.subplots(2, 3, figsize=(6.5, 4))


for gene, ax in zip(top_genes_hvg, axs.ravel()):

    plot_gene(gene, positions, scaled, ax=ax)

plt.savefig('top_genes_hvg.svg', dpi=600)
