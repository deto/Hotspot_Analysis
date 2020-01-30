import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import loompy
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'

hs = pd.read_table("../../CD4_w_protein/hotspot/hotspot_hvg.txt", index_col=0)
hvg_ens = pd.read_table("../../CD4_w_protein/genes/hvg.txt", header=None)[0].tolist()
hvg_info = pd.read_table("../../CD4_w_protein/genes/hvg_info.txt", index_col=0)
danb_info = pd.read_table("../../CD4_w_protein/genes/danb_info.txt", index_col=0)
pca_info = pd.read_table("../../CD4_w_protein/genes/pca_info.txt", index_col=0)
gr_pre = pd.read_table("../../CD4_w_protein/evaluation/geneRelevanceScores.txt", index_col=0)

hs_th = pd.read_table("../../CD4_w_protein/hotspot/hotspot_threshold.txt", index_col=0)
hs_pc = pd.read_table("../../CD4_w_protein/hotspot/hotspot_pca_threshold.txt", index_col=0)

# load gene map for ens -> symbol
ds = loompy.connect("../../../data/10x_PBMC_w_proteins/cd4/data.loom", "r")
gene_info = pd.DataFrame({
    "EnsID": ds.ra["EnsID"][:],
    "Symbol": ds.ra["Symbol"][:]
}).set_index("EnsID")
ens_map = {k: v for k, v in gene_info.itertuples()}
ds.close()

# %% Load background of gene_mean?
ds = loompy.connect("../../../data/10x_PBMC_w_proteins/cd4/data.loom", "r")
counts = pd.DataFrame(ds.layers['scaled'][:, :], index=gene_info.index)
ds.close()

gene_mean = counts.mean(axis=1)
gene_var = counts.var(axis=1)
gene_fano = gene_var / gene_mean
high_genes_ens = gene_mean.sort_values().index[-500:]
gene_mean = gene_mean / 5500 * 10000  # Make it counts/10k

# %% Load gene relevance

background = set(hs.Symbol.str.upper())

gr = {
    k.upper(): v for k, v in gr_pre.itertuples() if k.upper() in background
}

for x in background:
    if x not in gr:
        gr[x] = 0

gr_ens = {x: gr[ens_map[x].upper()] for x in hs.index}

# %% What about the mean ECDF as a function of threshold?

# np.exp(hvg_info['gene.mean']-1) = Mean in counts/10k

hvg_info = hvg_info.loc[hs.index]
danb_info = danb_info.loc[hs.index]
pca_info = pca_info.loc[hs.index]
data = pd.DataFrame({
    'Hotspot': hs.Z,
    'HSth': hs_th.Z,
    'HSpc': hs_pc.Z,
    'HVG': hvg_info['gene.dispersion.scaled'],
    'NBDisp': -1*np.log10(danb_info['q.value']+1e-300),
    'PCA': pca_info['Score'],
})

data['GR'] = [gr_ens[x] for x in data.index]
data = data.join(gene_mean.rename("Mean"))

# Remove highly-expressed genes
# Default behavior for Seurat-hvg
data = data.loc[data.Mean <= 20]

# %% First plot - All genes in at least 10 cells

plt.figure()

colors = sns.color_palette('deep')
colormap = {
    'Hotspot': colors[0],
    'HVG': colors[1],
    'NBDisp': colors[2],
    'PCA': colors[3],
}

def plot_individual(data, variable):
    d_sub = data[[variable, 'GR']].sort_values(variable, ascending=False)
    d_sub['Rank'] = np.arange(d_sub.shape[0])+1
    d_sub['GR_mean'] = d_sub['GR'].cumsum() / d_sub['Rank']

    plt.plot(d_sub['Rank'], d_sub['GR_mean'], '-',
             label=variable, color=colormap[variable])

plot_individual(data, 'Hotspot')
plot_individual(data, 'HVG')
plot_individual(data, 'NBDisp')
plot_individual(data, 'PCA')

# plt.xlim(5, data.shape[0])
plt.xlim(5, 3000)
plt.ylim(10, 50)
plt.xlabel("# of Genes")
plt.ylabel("Mean Gene-Relevance Score")

plt.legend()
plt.title('GR Scores\nAll Genes in at least 10 cells')
plt.grid(color='#dddddd', linestyle=(0, (5, 5)))
plt.gca().set_axisbelow(True)
#plt.show()
plt.savefig('GeneRelevance_noThresh.svg')

# %% Second plot - only genes above a threshold of np.exp(0.1) - 1

plt.figure()

thresh = np.exp(0.1) - 1  # To match default Seurat thresh

dataT = data.loc[data.Mean > thresh]

def plot_individual(data, variable):
    d_sub = data[[variable, 'GR']].sort_values(variable, ascending=False)
    d_sub['Rank'] = np.arange(d_sub.shape[0])+1
    d_sub['GR_mean'] = d_sub['GR'].cumsum() / d_sub['Rank']

    plt.plot(d_sub['Rank'], d_sub['GR_mean'], '-', label=variable,
             color=colormap[variable])

plot_individual(dataT, 'Hotspot')
plot_individual(dataT, 'HVG')
plot_individual(dataT, 'NBDisp')
plot_individual(dataT, 'PCA')


plt.xlim(5, dataT.shape[0])
plt.xlabel("# of Genes")
plt.ylabel("Mean Gene-Relevance Score")

plt.legend()
plt.title('GR Scores\nAll Genes With Mean Expression > 0.1 CP10K')
#plt.show()
plt.savefig('GeneRelevance_Thresh.svg')

# %% Compare different HS results here too:
# They are essentially all the same

plt.figure()

thresh = np.exp(0.1) - 1  # To match default Seurat thresh

dataT = data.loc[data.Mean > thresh]

def plot_individual(data, variable):
    d_sub = data[[variable, 'GR']].sort_values(variable, ascending=False)
    d_sub['Rank'] = np.arange(d_sub.shape[0])+1
    d_sub['GR_mean'] = d_sub['GR'].cumsum() / d_sub['Rank']

    plt.plot(d_sub['Rank'], d_sub['GR_mean'], '-', label=variable)

plot_individual(dataT, 'Hotspot')
plot_individual(dataT, 'HSth')
plot_individual(dataT, 'HSpc')

plt.xlim(5, dataT.shape[0])
plt.xlabel("# of Genes")
plt.ylabel("Mean Gene-Relevance Score")

plt.legend()
plt.savefig('GeneRelevance_HSComparison.svg')
