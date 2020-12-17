import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import loompy
import matplotlib.ticker as ticker
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'

k_values = [
    5, 10, 30, 50, 100, 300, 500, 1000
]

hs_results = {
    k: pd.read_table("../../CD4_w_protein/k_sensitivity/k_{k}/hotspot.txt".format(k=k), index_col=0)
    for k in k_values
}

gr_pre = pd.read_table("../../CD4_w_protein/evaluation/geneRelevanceScores.txt", index_col=0)

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
gene_mean = gene_mean / 5500 * 10000  # Make it counts/10k
gene_mean = gene_mean.rename('GeneMean')

# %% Load gene relevance

background = set(hs_results[k_values[0]].Symbol.str.upper())

gr = {
    k.upper(): v for k, v in gr_pre.itertuples() if k.upper() in background
}

for x in background:
    if x not in gr:
        gr[x] = 0

gr_ens = {x: gr[ens_map[x].upper()] for x in hs_results[k_values[0]].index}


def gr_score(genes):
    scores = [gr_ens[x] for x in genes]
    return np.mean(scores)


def gr_score_k(k, N, hs_results, gene_mean):
    """
    For the hotspot with k=k results, take the top
    N genes and return the average gene relevance score
    (after excluding highly-expressed genes)
    """

    hs = hs_results[k]
    hs = hs.join(gene_mean)
    hs = hs.loc[hs['GeneMean'] <= 20]
    hs = hs.sort_values('Z', ascending=False)

    genes = hs.index[0:N]
    return gr_score(genes)


results = pd.Series(
    {
        k: gr_score_k(k, N=1000, hs_results=hs_results, gene_mean=gene_mean)
        for k in k_values
    }
)


plt.figure()
ax = plt.gca()
plt.plot(results.index, results.values, 'o-')
plt.xscale('log')
plt.xticks(results.index)
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax.xaxis.set_ticks([], minor=True)
plt.xlabel('# of Neighbors')
plt.ylabel('GR Score of top\n100 Genes (Mean)')
plt.gca().set_axisbelow(True)
plt.grid(color='#CCCCCC', ls=(0, (5, 5)), lw=.5)
plt.savefig('GeneRelevance_K_sensitivity.svg', dpi=300)

# %% Also, just the percent overlap with K=30 top 1000 genes result

def extract_top_genes(k, N, hs_results, gene_mean):
    """
    Take the top N genes for the hotspot results with
    run with num_neighbors=k
    """

    hs = hs_results[k]
    hs = hs.join(gene_mean)
    hs = hs.loc[hs['GeneMean'] <= 20]
    hs = hs.sort_values('Z', ascending=False)

    genes = hs.index[0:N]
    return genes

# target_k = 100
# overlap = {}
# 
# target_genes = extract_top_genes(target_k, N=1000, hs_results=hs_results, gene_mean=gene_mean)
# 
# for k in k_values:
# 
#     k_genes = extract_top_genes(k, N=1000, hs_results=hs_results, gene_mean=gene_mean)
# 
#     intersection = len(set(k_genes) & set(target_genes))
#     union = len(set(k_genes) | set(target_genes))
# 
#     overlap[k] = intersection / union
# 
# overlap = pd.Series(overlap)
# 
# 
# plt.figure()
# ax = plt.gca()
# plt.plot(overlap.index, overlap.values, 'o')
# plt.xscale('log')
# plt.xticks(overlap.index)
# ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
# ax.xaxis.set_ticks([], minor=True)
# plt.xlabel('# of Neighbors')
# plt.ylabel('Jaccard Overlap with k={}'.format(target_k))
# plt.show()


# %% Or, what if we look at the average rank of the top 1000 genes for k=100 vs k=k

def average_rank(k, N, target_genes, hs_results, gene_mean):
    """
    Take the top N genes for the hotspot results with
    run with num_neighbors=k
    """

    hs = hs_results[k]
    hs = hs.join(gene_mean)
    hs = hs.loc[hs['GeneMean'] <= 20]
    hs = hs.sort_values('Z', ascending=False)
    hs['Rank'] = hs.Z.rank(ascending=False)

    return hs.loc[target_genes]['Rank'].mean()

target_k = 100

avg_rank = {}

target_genes = extract_top_genes(target_k, N=1000, hs_results=hs_results, gene_mean=gene_mean)

for k in k_values:

    k_genes = extract_top_genes(k, N=1000, hs_results=hs_results, gene_mean=gene_mean)

    avg_rank[k] = average_rank(
        k=k, N=1000, target_genes=target_genes, hs_results=hs_results, gene_mean=gene_mean
    )

avg_rank = pd.Series(avg_rank)

plt.figure()
ax = plt.gca()
plt.plot(avg_rank.index, avg_rank.values, 'o-')
plt.xscale('log')
plt.xticks(avg_rank.index)
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax.xaxis.set_ticks([], minor=True)
plt.xlabel('# of Neighbors')
plt.ylabel('Average Rank of\nHS(k={}) genes'.format(target_k))
plt.gca().set_axisbelow(True)
plt.grid(color='#CCCCCC', ls=(0, (5, 5)), lw=.5)
plt.savefig('AverageRank_K_sensitivity.svg', dpi=300)
