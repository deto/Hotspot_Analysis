"""
This file plots the tree next to the gene modules' expression
It doesn't look that good, though.  These associations
don't visualize that well as they are subtle
"""
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
from scipy.cluster.hierarchy import linkage, leaves_list
from hotspot.local_stats_pairs import create_centered_counts
from Bio import Phylo
import loompy


# %% Load data
bptree = Phylo.read(
    "../../data/Lineage/Embryo3/tree.txt",
    format='newick'
)

loom_file = "../../data/Lineage/Embryo3/data.loom"
module_file = "../Embryo3/hotspot/modules_lineage_tree.txt"

with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

modules = pd.read_table(module_file, index_col=0).Cluster

# Have to do this because data_slideseq makes it a numpy array
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)
num_umi = pd.Series(num_umi, index=barcodes)

# %%

# Align counts to tree and center data
leaf_order = [x.name for x in bptree.clade.get_terminals()]
counts = counts.loc[:, leaf_order]
num_umi = num_umi[leaf_order]

cplot = create_centered_counts(
    counts.values, 'danb',
    num_umi.values
)
cdata = pd.DataFrame(
    cplot.T, index=counts.columns, columns=counts.index
)

# %% Sort the tree a bit better when order is ambiguous

def order_clade(clade):

    if clade.is_terminal():
        return (1, cdata.loc[clade.name])

    results = [order_clade(x) for x in clade.clades]
    
    if len(results) == 1:
        return results[0]

    cl_size = pd.Series(
        [r[0] for r in results],
        index=[r[1].name for r in results]
    )
    cl_exp = pd.concat([r[1] for r in results], axis=1).T

    ii = leaves_list(linkage(cl_exp.values, metric='cosine', method='average'))

    new_clades = [clade.clades[i] for i in ii]
    clade.clades = new_clades

    new_size = cl_size.sum()
    new_exp = cl_exp.multiply(cl_size, axis=0).mean(axis=0) / new_size
    new_exp.name = results[0][1].name


    return new_size, new_exp


final_size, final_exp = order_clade(bptree.clade)

# %%

zdata = np.log2(counts.divide(counts.sum(axis=0), axis=1)*10000 + 1)
zdata = zdata.subtract(zdata.mean(axis=1), axis=0)
zdata = zdata.divide(zdata.std(axis=1), axis=0)

# %%

x3 = [x.name for x in bptree.clade.get_terminals()]
gene_order = []

for mod in np.sort(modules.unique()):
    if mod == -1: continue

    mod_genes = modules.index[modules == mod]
    sampled_genes = np.random.choice(
        mod_genes, min(50, len(mod_genes)), replace=False
    )

    # Order within the group:
    ii = leaves_list(linkage(zdata.loc[sampled_genes, :], metric='cosine'))

    gene_order.extend(sampled_genes[ii])

fig, axs = plt.subplots(
    1, 2, sharey=True, gridspec_kw=dict(width_ratios=[1, 3])
)

matplotlib.rcParams['lines.linewidth'] = .5
plt.sca(axs[0])
Phylo.draw(bptree, do_show=False, show_confidence=False,
           label_func=lambda x: None,
           axes=axs[0])
plt.xlim(-1, 13)
plt.xticks([])
plt.yticks([])
for sp in plt.gca().spines.values():
    sp.set_visible(False)
plt.xlabel('')
plt.ylabel('')

# need to flip y axes .... I think
plt.sca(axs[1])

pdata = zdata.loc[gene_order, x3].T.values

plt.pcolormesh(
    np.arange(pdata.shape[1] + 1) + 0.5,
    np.arange(pdata.shape[0] + 1) + 0.5,
    pdata,
    vmin=-1,
    vmax=1,
    cmap="RdBu_r",
)
plt.xticks([])
plt.yticks([])
for sp in plt.gca().spines.values():
    sp.set_visible(False)
# plt.yticks(np.arange(cplot.shape[1])+1, leaf_order)
plt.show()

# %% 
import seaborn as sns
sns.clustermap(
    zdata.loc[gene_order, x3].T,
    vmin=-1, vmax=1, cmap='RdBu_r',
    col_cluster=False, row_cluster=False
)
plt.show()
