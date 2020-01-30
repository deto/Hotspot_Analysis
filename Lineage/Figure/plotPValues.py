"""
Plotting p-value distributions for hotspot - on transcriptional data
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
from scipy.stats import norm

plt.rcParams['svg.fonttype'] = 'none'

tab10 = plt.get_cmap('tab10').colors

hs_data = pd.read_table(
    "../Embryo3/hotspot/hotspot_lineage_tree.txt", index_col=0
)

hs_data_shuffled = pd.read_table(
    "../Embryo3/hotspot/hotspot_tree_shuffled.txt", index_col=0
)


# This looks better now with the shuffled data
# Need to also limit on the high end?
#   No that's not it.  What's wrong?

# %% Q-Q plot instead?


def qq_plot(vals):
    sorted_pvals = np.sort(vals)

    N = len(sorted_pvals)
    expected = np.linspace(1/(2*N), 1-1/(2*N), N)

    plt.plot(expected, sorted_pvals, 'o-', ms=1)
    plt.plot(expected, expected, '--', linewidth=1)


def log_qq_plot(vals):
    sorted_pvals = np.sort(vals)

    min_nnz = sorted_pvals[sorted_pvals > 0].min()
    sorted_pvals[sorted_pvals == 0] = min_nnz

    N = len(sorted_pvals)
    expected = np.linspace(1/(2*N), 1-1/(2*N), N)

    plt.plot(np.log10(expected), np.log10(sorted_pvals), 'o-', ms=3, rasterized=True)
    plt.plot(np.log10(expected), np.log10(expected), '--', linewidth=1, rasterized=True)


def z_qq_plot(vals):
    sorted_zvals = np.sort(vals)

    N = len(sorted_zvals)
    expected = np.linspace(1/(2*N), 1-1/(2*N), N)
    expected = norm.ppf(expected)

    plt.plot(expected, sorted_zvals, 'o-', ms=3, rasterized=True)
    plt.plot(expected, expected, '--', linewidth=1, rasterized=True)


from statsmodels.stats.multitest import multipletests
false_pos = multipletests(hs_data_shuffled.Pval.values, alpha=0.05, method='fdr_bh')[0].sum()

plt.figure(figsize=(4.5, 4.5))
z_qq_plot(hs_data_shuffled.Z.values)
plt.xlabel('Theoretical Quantiles')
plt.ylabel('Actual Quantiles')
plt.title('Hotspot Feature Selection\nQ-Q Plot')
plt.text(-3, 3, 'sum(FDR < .05) = {}'.format(false_pos))
plt.subplots_adjust(left=.2, bottom=.2)
# plt.show()
plt.savefig('hotspot_qq_fs.svg', dpi=300)

# %% Look at pair-wise values then


hs_pairs = pd.read_table(
    "../Embryo3/hotspot/hotspot_pairs_shufled_z.txt.gz", index_col=0
)

hs_pairs_z = hs_pairs.values.ravel()
hs_pairs_p = norm.sf(hs_pairs_z)

plt.figure(figsize=(4.5, 4.5))
# log_qq_plot(hs_pairs_p)
z_qq_plot(hs_pairs_z)
plt.xlabel('Theoretical Quantiles')
plt.ylabel('Actual Quantiles')
plt.title('Hotspot Pairwise Statistic\nQ-Q Plot')
plt.subplots_adjust(left=.2, bottom=.2)
# plt.show()
plt.savefig('hotspot_qq_pairs.svg', dpi=300)

# # %% What's going on, need to compare with data
# import loompy
# loom_file = "../../../data/10x_PBMC_w_proteins/cd4/data.loom"
# 
# with loompy.connect(loom_file, 'r') as ds:
#     barcodes = ds.ca['Barcode'][:]
#     counts = ds[:, :]
#     gene_info = ds.ra['EnsID', 'Symbol']
#     num_umi = ds.ca['NumUmi'][:]
# 
# # Have to do this because data_slideseq makes it a numpy array
# gene_info = pd.DataFrame(
#     gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
# counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)
# num_umi = pd.Series(num_umi, index=barcodes)
# 
# # Look at the bad shuffled values
# ii = hs_data_shuffled.sort_values('Z').index[-1]
# cc = counts.loc[ii]
# 
# # %% how do these agree with the number of umis?
# 
# pos = cc > 0
# 
# kwargs = dict(
#     bins=50,
#     range=(11, 15),
#     alpha=.5,
#     density=True
# )
# 
# plt.figure()
# plt.hist(np.log2(num_umi[pos]), **kwargs)
# plt.hist(np.log2(num_umi[~pos]), **kwargs)
# plt.show()
# 
# # %% What if we try many shufflings?  What is it about these values?
# 
# from hotspot import sim_data
# import hotspot
# 
# c_perm = sim_data.generate_permutation_null(cc, 1000)
# c_perm = pd.DataFrame(
#     c_perm,
#     columns=counts.columns,
#     index=['Gene{}'.format(i) for i in range(c_perm.shape[0])]
# )
# 
# latent = pd.read_table(
#     "../../CD4_w_protein/scvi/threshold/latent_shuffled.txt",
#     index_col=0
# )
# 
# 
# hs = hotspot.Hotspot(c_perm, latent, num_umi)
# 
# hs.create_knn_graph(
#     weighted_graph=True, n_neighbors=30, neighborhood_factor=3
# )
# 
# results = hs.compute_hotspot(model='danb', jobs=5, centered=True)
# 
# # This is fine.  But, this breaks any association between num UMI and the counts
# # %% What if we shuffle them together
# 
# 
# res = []
# 
# for _ in range(50):
#     shuf_i = np.random.permutation(cc.size)
# 
#     c_perm = cc.iloc[shuf_i].to_frame().T
# 
#     num_umi_perm = num_umi.iloc[shuf_i]
# 
#     latent = pd.read_table(
#         "../../CD4_w_protein/scvi/threshold/latent_shuffled.txt",
#         index_col=0
#     )
# 
#     hs = hotspot.Hotspot(c_perm, latent, num_umi_perm)
# 
#     hs.create_knn_graph(
#         weighted_graph=True, n_neighbors=30, neighborhood_factor=3
#     )
# 
#     results = hs.compute_hotspot(model='danb', jobs=5, centered=True)
#     res.append(results)
# 
# zs = [x.Z[0] for x in res]
# 
# plt.figure()
# plt.hist(zs, 30)
# plt.show()
# 
# shuf_old = hs_data_shuffled
# shuf_old2 = hs_data_shuffled
# 
# # %%
# 
# dat = pd.concat(
#     (
#         shuf_old.Z.rename('z1'),
#         shuf_old2.Z.rename('z2').loc[shuf_old.index]
#     ), axis=1
# )
# 
# plt.figure()
# plt.plot(dat.z1, dat.z2, 'o')
# plt.xlabel('One shuffling')
# plt.ylabel('Another shuffling')
# plt.show()
# 
# # %% bin the genes
# 
# zz = (counts.loc[shuf_old2.index] > 0).sum(axis=1)
# valid = zz.index[zz > 22]
# #gene_means = counts.loc[shuf_old2.index].mean(axis=1)
# gene_means = counts.loc[valid].mean(axis=1)
# gene_mean_bins = pd.qcut(gene_means, 9)
# 
# fig, axs = plt.subplots(3, 3)
# 
# 
# for ax, cat in zip(
#         axs.ravel(), gene_mean_bins.cat.categories):
# 
#     plt.sca(ax)
#     print(cat)
#     genes_ii = gene_mean_bins.index[gene_mean_bins == cat]
#     #plt.plot(dat.z1[genes_ii], dat.z2[genes_ii], 'o')
#     log_qq_plot(shuf_old2.loc[genes_ii].Pval.values)
#     plt.title(cat)
# 
# plt.show()
# 
# # %% What if we look at genes that are more highly expressed but still inflated?
# 
# shuf_old.loc[
#     gene_means.index[gene_means > .08]
# ].sort_values('Z')
# 
# cat = gene_mean_bins.cat.categories[6]
# genes_ii = gene_mean_bins.index[gene_mean_bins == cat]
