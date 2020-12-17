import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import loompy
from bio_utils.plots import hover_plot
from matplotlib.colors import LinearSegmentedColormap
import seaborn as sns
plt.rcParams['svg.fonttype'] = 'none'

# %% Load results

hs = pd.read_table("../../CD4_w_protein/hotspot/hotspot_hvg.txt", index_col=0)
hvg = pd.read_table("../../CD4_w_protein/genes/hvg_info.txt", index_col=0)
hvg = hvg.loc[hs.index]


# %% Plot one vs other

# plot_data = hs[['Symbol', 'Z']].join(
#     hvg[['gene.dispersion.scaled', 'gene.mean']]
# )
# 
# plt.figure(figsize=(5, 5))
# 
# plt.plot(
#     plot_data['gene.dispersion.scaled'],
#     plot_data['Z'],
#     'o', ms=2, rasterized=True
# )
# 
# plt.gca().set_axisbelow(True)
# plt.xlabel('Scaled Dispersion')
# plt.ylabel('Autocorrelation Z')
# plt.grid(color='#BBBBBB', ls=(0, (5, 5)), lw=.5)
# plt.subplots_adjust(left=0.15, right=1-0.15, bottom=0.15, top=1-0.15)
# 
# plt.show()
# plt.savefig('autocorrelation_vs_hvg.svg', dpi=300)


# %%

loom_file = "../../../data/10x_PBMC_w_proteins/cd4/data.loom"
with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    scaled = ds.layers['scaled'][:, :]
    counts = ds.layers[''][:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

# Have to do this because data_slideseq makes it a numpy array
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
scaled = pd.DataFrame(scaled, index=gene_info.index, columns=barcodes)
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)

umap = pd.read_table(
    "../../CD4_w_protein/umap/umap_hvg.txt", index_col=0
)

umap = umap.loc[scaled.columns]

embedding = pd.read_table(
    "../../CD4_w_protein/scvi/hvg/latent.txt.gz", index_col=0)

embedding = embedding.loc[scaled.columns]

# %%
gene_cmap = LinearSegmentedColormap.from_list(
    "grays", ['#cccccc', '#000000']
)


def plot_gene(gene, gene_info, umap, scaled, ax=None):

    if ax is None:
        ax = plt.gca()
    else:
        plt.sca(ax)

    gene_ens = gene_info.loc[gene_info.Symbol == gene].index[0]
    vals = np.log2(scaled.loc[gene_ens]+1)

    vmin = 0
    vmax = np.percentile(vals, 80)

    plot_data = pd.DataFrame({
        'X': umap.umap1,
        'Y': umap.umap2,
        'Expression': vals,
    }).sample(frac=1)

    plt.scatter(
        x=plot_data.X,
        y=plot_data.Y,
        c=plot_data.Expression,
        s=4,
        vmin=vmin, vmax=vmax, cmap=gene_cmap,
        alpha=0.7, rasterized=True, edgecolors='none',
    )
    for sp in ax.spines.values():
        sp.set_visible(False)

    plt.xticks([])
    plt.yticks([])
    plt.title(gene)


# %%

plot_data = hs[['Symbol', 'Z']].join(hvg['gene.dispersion.scaled'])
plot_data = plot_data.rename({
    'Z': 'Hotspot',
    'gene.dispersion.scaled': 'HVG',
}, axis=1)

to_showcase = [
    'TCF7', 'TGFB1', 'FOSB', 'SELL', 'CD81', 'ANXA1'
]

plt.figure()
sns.distplot(plot_data['Hotspot'], bins=50,
             hist_kws={'range': (-5, 40)}, kde_kws={'bw': .3, 'gridsize': 200})
plt.xlim(-5, 40)

plt.autoscale(False)

ymax = plt.gca().get_ylim()[1]

for gene in to_showcase:
    x = plot_data.loc[plot_data['Symbol'] == gene, 'Hotspot'][0]
    plt.vlines(x, -1e9, ymax*.5)
    plt.annotate(
        gene,
        (x, ymax*.5), (0, 30),
        va='center', ha='center', textcoords='offset pixels',
        arrowprops={'arrowstyle': '-'},
    )

sns.despine(top=True, right=True, left=True)
plt.yticks([])
plt.xlabel('Hotspot Z-score')

plt.savefig('hs_vs_hvg_hs_score.svg', dpi=300)

# %%

plt.figure()
sns.distplot(plot_data['HVG'], bins=50,
             hist_kws={'range': (-3, 15)}, kde_kws={'bw': .3, 'gridsize': 200})
plt.xlim(-3, 15)

plt.autoscale(False)

ymax = plt.gca().get_ylim()[1]

for gene in to_showcase:
    x = plot_data.loc[plot_data['Symbol'] == gene, 'HVG'][0]
    plt.vlines(x, -1e9, ymax*.5)
    plt.annotate(
        gene,
        (x, ymax*.5), (0, 30),
        va='center', ha='center', textcoords='offset pixels',
        arrowprops={'arrowstyle': '-'},
    )

sns.despine(top=True, right=True, left=True)
plt.yticks([])
plt.xlabel('HVG Scaled Dispersion')

plt.savefig('hs_vs_hvg_hvg_score.svg', dpi=300)

# %%

fig, axs = plt.subplots(3, 2, figsize=(6.2, 9.2))

for gene, ax in zip(to_showcase, axs.ravel()):

    plot_gene(gene, gene_info, umap, scaled, ax=ax)

plt.subplots_adjust(top=.9, right=.9, left=.1, bottom=.1)
plt.savefig('hs_vs_hvg_umaps.svg', dpi=300)
# plt.show()
