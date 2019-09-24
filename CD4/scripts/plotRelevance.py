import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import loompy

hs = pd.read_table("../hotspot/hotspot.txt", index_col=0)
hvg_ens = pd.read_table("../genes/hvg.txt", header=None)[0].tolist()
gr_pre = pd.read_table("../evaluation/geneRelevanceScores.txt", index_col=0)

# load gene map for ens -> symbol
ds = loompy.connect("../../data/10x_PBMC_w_proteins/cd4/data.loom", "r")
gene_info = pd.DataFrame({
    "EnsID": ds.ra["EnsID"][:],
    "Symbol": ds.ra["Symbol"][:]
}).set_index("EnsID")
ens_map = {k: v for k, v in gene_info.itertuples()}
ds.close()

# %% Load background of gene_mean?
ds = loompy.connect("../../data/10x_PBMC_w_proteins/cd4/data.loom", "r")
counts = pd.DataFrame(ds.layers['scaled'][:, :], index=gene_info.index)
ds.close()

gene_mean = counts.mean(axis=1)
high_genes_ens = gene_mean.sort_values().index[-500:]

# %% Create sets

background = set(hs.Symbol.str.upper())
hs_genes = hs.loc[hs.FDR < .05].Symbol.str.upper()
#hs_genes = hs.loc[hs.FDR < 1e-25].Symbol.str.upper()
hvg_genes = {ens_map[x].upper() for x in hvg_ens}
high_genes = {ens_map[x].upper() for x in high_genes_ens}

gr = {
    k.upper(): v for k, v in gr_pre.itertuples() if k.upper() in background
}

for x in background:
    if x not in gr:
        gr[x] = 0

# %% Plot ECDFs


def get_ecdf(genes, gr):

    vals = pd.Series({x: gr[x] for x in genes})
    vals = vals.sort_values()

    x = vals.values

    N = len(x)
    y = np.linspace(1/N, 1, N)

    return x, y


plt.figure()
x, y = get_ecdf(hs_genes, gr)
plt.plot(1-y, x, label='Hotspot')

x, y = get_ecdf(hvg_genes, gr)
plt.plot(1-y, x, label='Variable Genes')

x, y = get_ecdf(high_genes, gr)
plt.plot(1-y, x, label='High')

x, y = get_ecdf(background, gr)
plt.plot(1-y, x, label='Background')

plt.legend()

plt.ylim(0, max(gr.values()))
plt.xlim(0, 1.001)
plt.ylabel("Number of Overlapping Subsets")
plt.xlabel("ECDF")

plt.show()


# %% How is gene relevance related to mean expression?

x = gene_mean.sort_values()
y = []
for g in gene_mean.index:
    sym = ens_map[g]
    try:
        gr_score = gr[sym]
    except KeyError:
        gr_score = 0
    y.append(gr_score)

data = pd.DataFrame({
    'Mean': np.log10(x+1e-4),
    'Relevance': y,
})

data['Bin'] = pd.qcut(data['Mean'], 30)

bin_means = data.groupby('Bin')['Relevance'].mean()
mids = [x.mid for x in bin_means.index]

plt.figure()
plt.plot(data['Mean'], data['Relevance'], 'o')
plt.show()

plt.figure()
plt.plot(mids, bin_means, 'o')
plt.show()
