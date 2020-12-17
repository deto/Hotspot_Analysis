import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from tqdm import tqdm
from scipy.stats import spearmanr

plt.rcParams['svg.fonttype'] = 'none'

# SIMULATION_DIR = "../../../data/Simulated4_downsampling/"
# SIMULATION_RESULTS_DIR = "../../Simulation4_downsampling/"
# SIMULATION_DIR = "../../../data/Simulated5/"
# SIMULATION_RESULTS_DIR = "../../Simulation5/"
SIMULATION_DIR = "../../../data/Simulated6/"
SIMULATION_RESULTS_DIR = "../../Simulation6_downsampling/"
# SIMULATION_DIR = "../../../data/Simulated_Tree_Continuous/"
# SIMULATION_RESULTS_DIR = "../../Simulation_Tree_Continuous_downsampling/"

# %% First, need to compute the 'true' correlation matrix


true_counts = pd.read_table(
    SIMULATION_DIR + "true_counts/true_counts.txt.gz",
    index_col=0
)

gene_effects = pd.read_table(
    SIMULATION_DIR + "true_counts/gene_effects_s.txt",
    index_col=0
)

ge_size = gene_effects.iloc[:, 1:].abs().max(axis=1)
ge_size = ge_size.loc[ge_size != 0]
genes = ge_size.index.tolist()

true_counts_sub = true_counts.loc[genes]

true_corr = true_counts_sub.T.corr()
true_corr.index.name = 'GeneX'
true_corr.columns.name = 'GeneY'

# %% Plot true correlation - does it look reasonable?


# sns.clustermap(true_corr, vmin=-.2, vmax=.2, cmap="RdBu_r")
# plt.show()
# 
# sns.clustermap(true_corr.abs(), vmin=0, vmax=.2)
# plt.show()
# 
# # %% Plot this, but with the gene effects next to it
# 
# from matplotlib.colors import Normalize
# from matplotlib.cm import ScalarMappable
# 
# cmap = plt.get_cmap('RdBu_r')
# norm = Normalize(vmin=-2, vmax=2)
# sm = ScalarMappable(norm, cmap)
# 
# 
# ge_sub = gene_effects.loc[true_corr.index]
# 
# row_colors = pd.DataFrame(index=ge_sub.index)
# for col in ge_sub:
#     row_colors[col] = [sm.to_rgba(x) for x in ge_sub[col]]
# 
# sns.clustermap(
#     true_corr, vmin=-.2, vmax=.2, cmap="RdBu_r",
#     row_colors=row_colors
# )
# plt.show()

# %% Ok, load some downsampled correlations


downsample_corr = pd.read_table(
    SIMULATION_RESULTS_DIR + "/5e3/hotspot/regular_pairs_lc.txt.gz",
    index_col=0
)
downsample_corr.index.name = 'GeneX'
downsample_corr.columns.name = 'GeneY'

hotspot_corr = pd.read_table(
    SIMULATION_RESULTS_DIR + "5e3/hotspot/hotspot_pairs_z.txt.gz",
    index_col=0
)
hotspot_corr.index.name = 'GeneX'
hotspot_corr.columns.name = 'GeneY'


# %% Compare true vs downsampled
def compare_spearmanr(seq_depth):

    downsample_corr = pd.read_table(
        (SIMULATION_RESULTS_DIR + "{}/hotspot/regular_pairs_lc.txt.gz")
        .format(seq_depth), index_col=0
    )
    downsample_corr.index.name = 'GeneX'
    downsample_corr.columns.name = 'GeneY'

    hotspot_corr = pd.read_table(
        (SIMULATION_RESULTS_DIR + "{}/hotspot/hotspot_pairs_z.txt.gz")
        .format(seq_depth), index_col=0
    )
    hotspot_corr.index.name = 'GeneX'
    hotspot_corr.columns.name = 'GeneY'

    downsample_corr_long = downsample_corr.reset_index() \
        .melt('GeneX', value_name='DownsampleCorr')
    hotspot_corr_long = hotspot_corr.reset_index() \
        .melt('GeneX', value_name='HotspotCorr')
    true_corr_long = true_corr.reset_index() \
        .melt('GeneX', value_name='TrueCorr')

    plot_data = (
        true_corr_long
        .merge(downsample_corr_long)
        .merge(hotspot_corr_long)
    )

    corr_result = spearmanr(plot_data['TrueCorr'], plot_data['DownsampleCorr'])
    hs_result = spearmanr(plot_data['TrueCorr'], plot_data['HotspotCorr'])

    return pd.Series({
        'RegCorrR': corr_result.correlation,
        'RegCorrP': corr_result.pvalue,
        'HsCorrR': hs_result.correlation,
        'HsCorrP': hs_result.pvalue,
        'SeqDepth': seq_depth
        })


seq_depths = ["1e3", "2e3", "5e3", "10e3", "20e3", "40e3"]

results = pd.concat(
    [compare_spearmanr(x) for x in tqdm(seq_depths)],
    axis=1).T

reg_results = results[['RegCorrR', 'RegCorrP', 'SeqDepth']].copy()
reg_results['Method'] = 'Pearsons'

hs_results = results[['HsCorrR', 'HsCorrP', 'SeqDepth']].copy()
hs_results['Method'] = 'Hotspot'

reg_results.columns = ['CorrR', 'CorrP', 'SeqDepth', 'Method']
hs_results.columns = ['CorrR', 'CorrP', 'SeqDepth', 'Method']

plot_data = pd.concat((reg_results, hs_results), axis=0)

# %%

plt.figure(figsize=(6.5, 5))
sns.barplot(
    data=plot_data, x='SeqDepth', y='CorrR', hue='Method',
    order=seq_depths[::-1], hue_order=['Hotspot', 'Pearsons'],
    alpha=0.9)
plt.xlabel('Simulated Sequencing Depth\n(Thousands of Reads/Cell)')
plt.ylabel("Comparison to True Correlation\nSpearman's $\\rho$")
plt.xticks(np.arange(len(seq_depths)), [x.split('e')[0] for x in seq_depths[::-1]])
plt.subplots_adjust(bottom=0.15)
plt.gca().set_axisbelow(True)
plt.grid(axis='y', color='#AAAAAA', ls=(0, (5, 5)), zorder=1)
# plt.show()
plt.savefig('downsampling_correlation.svg', dpi=300)

# %% How to the matrices look - head-to-head

from scipy.cluster.hierarchy import leaves_list, linkage

downsample_corr = pd.read_table(
    SIMULATION_RESULTS_DIR + "5e3/hotspot/regular_pairs_lc.txt.gz",
    index_col=0
)
downsample_corr.index.name = 'GeneX'
downsample_corr.columns.name = 'GeneY'

hotspot_corr = pd.read_table(
    SIMULATION_RESULTS_DIR + "5e3/hotspot/hotspot_pairs_z.txt.gz",
    index_col=0
)
hotspot_corr.index.name = 'GeneX'
hotspot_corr.columns.name = 'GeneY'


Z_true = linkage(true_corr.values, method='average', metric='euclidean')
gene_ii = leaves_list(Z_true)
genes = true_corr.index[gene_ii].tolist()

# %% Plot all three side-by-side

fig, axs = plt.subplots(
    1, 3, figsize=(13, 4.5),
    # sharex=True, sharey=True
)

cbaxs = [
    fig.add_axes([.29+.95/3*i, .1, .01, .2]) for i in range(3)
]

cbar_kws = {
    'fraction': .1, 'shrink': .5, 'panchor': (1.0, 1.0),
}


def add_title(cbar_kws, title):
    cbar_kws = cbar_kws.copy()
    cbar_kws['label'] = title
    return cbar_kws


plt.sca(axs[0])
vmax = np.percentile(true_corr.abs().values.ravel(), 95)
sns.heatmap(
    true_corr.iloc[gene_ii, gene_ii], vmin=vmax*-1, vmax=vmax, cmap="RdBu_r",
    cbar_kws=add_title(cbar_kws, r"Pearson's $\rho$"), cbar_ax=cbaxs[0],
    rasterized=True,
)
plt.xticks([])
plt.yticks([])
plt.xlabel('')
plt.ylabel('')
plt.title('Correlations: True Counts')

plt.sca(axs[1])
vmax = np.percentile(downsample_corr.abs().values.ravel(), 95)
sns.heatmap(
    downsample_corr.loc[genes, genes], vmin=-1*vmax, vmax=vmax, cmap="RdBu_r",
    cbar_kws=add_title(cbar_kws, r"Pearson's $\rho$"), cbar_ax=cbaxs[1],
    rasterized=True,
)
plt.xticks([])
plt.yticks([])
plt.xlabel('')
plt.ylabel('')
plt.title("Correlations: Downsampled Counts")

plt.sca(axs[2])
vmax = np.percentile(hotspot_corr.abs().values.ravel(), 95)
sns.heatmap(
    hotspot_corr.loc[genes, genes], vmin=-1*vmax, vmax=vmax, cmap="RdBu_r",
    cbar_kws=add_title(cbar_kws, 'Z-Score'), cbar_ax=cbaxs[2],
    rasterized=True,
)
plt.xticks([])
plt.yticks([])
plt.xlabel('')
plt.ylabel('')
plt.title("Hotspot Z: Downsampled Counts")
plt.subplots_adjust(left=0.05, right=0.9, wspace=0.45)
plt.savefig('downsampled_correlation_matrices.svg', dpi=300)
# plt.show()
