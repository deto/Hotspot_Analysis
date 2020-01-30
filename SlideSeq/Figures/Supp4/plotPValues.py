"""
Plotting p-value distributions for hotspot and spatialDE

Both the general distributions and the distribution of 'positive'
genes taken from the DropViz data
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
from scipy.stats import norm

plt.rcParams['svg.fonttype'] = 'none'

tab10 = plt.get_cmap('tab10').colors

hs_data = pd.read_table(
    "../../Puck_180819_12/hotspot/hotspot.txt", index_col=0
)

hs_data_shuffled = pd.read_table(
    "../../Puck_180819_12/hotspot/hotspot_shuffled.txt", index_col=0
)

# %% Load positives:
pos_data = pd.read_table(
    "../../DropViz/edger_markers_1vAll.txt"
)
pos_data.head()

N = 100 # Number of genes per cluster

markers = {}
for c_id, group in pos_data.groupby('Cluster'):
    group = group.loc[group.FDR < .1]
    m = group.sort_values('LR', ascending=False).GeneSymbol[0:N]
    m = set(m)
    markers[c_id] = m

# Collapse markers
from functools import reduce

markers_all = reduce(set.union, markers.values(), set())
markers_all = {x.lower() for x in markers_all}


# %% Plot all distributions

is_marker = [
    x.lower() in markers_all for x in hs_data.Symbol
]

plt.figure()
plt.hist(hs_data.Pval, 100, range=(0, 1), density=False, alpha=0.6, label='All')
plt.hist(hs_data.Pval[is_marker], 100, range=(0, 1), density=False, alpha=0.6, label='Markers')
plt.hist(hs_data_shuffled.Pval, 100, range=(0, 1), density=False, alpha=0.6, label='Shuffled')
plt.legend()
plt.xlabel('P-Value')
plt.ylabel('PDF')
# plt.show()
plt.savefig('p_value_histograms.svg')

# %% What about for SDE?
sd_data = pd.read_table(
    "../../Puck_180819_12/spatialDE/spatialDE.txt", index_col=0
)

sd_data_shuffled = pd.read_table(
    "../../Puck_180819_12/spatialDE/spatialDE_shuffled.txt", index_col=0
)

is_marker = [
    x.lower() in markers_all for x in sd_data.index
]

plt.figure()
plt.hist(sd_data.pval, 100, range=(0, 1), density=False, alpha=0.6, label='All')
plt.hist(sd_data.pval[is_marker], 100, range=(0, 1), density=False, alpha=0.6, label='Markers')
plt.hist(sd_data_shuffled.pval, 100, range=(0, 1), density=False, alpha=0.6, label='Shuffled')
plt.legend()
plt.xlabel('P-Value')
plt.ylabel('PDF')
# plt.show()
plt.savefig('p_value_histograms_sde.svg')

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

# %%

from statsmodels.stats.multitest import multipletests

false_pos = multipletests(hs_data_shuffled.Pval.values, alpha=0.05, method='fdr_bh')[0].sum()

plt.figure(figsize=(4.5, 4.5))
z_qq_plot(hs_data_shuffled.Z.values)
plt.xlabel('Theoretical Quantiles')
plt.ylabel('Sample Quantiles')
plt.title('Hotspot Feature Selection\nQ-Q Plot')
plt.text(-3, 3, 'sum(FDR < .05) = {}'.format(false_pos))
plt.subplots_adjust(left=.2, bottom=.2)
# plt.show()
plt.savefig('hotspot_qq_fs.svg', dpi=300)

# %%

plt.figure(figsize=(4.5, 4.5))
log_qq_plot(sd_data_shuffled.pval.values)
plt.xlabel('Theoretical Quantiles\n($-log_{10}$ p-value)')
plt.ylabel('Actual Quantiles\n($-log_{10}$ p-value)')
plt.title('SpatialDE Feature Selection\nQ-Q Plot')
plt.subplots_adjust(left=.2, bottom=.2)
# plt.show()
plt.savefig('sde_qq_fs.svg', dpi=300)


# %% Now, what about the spectrum and heatmaps for hotspot pairwise values?

hs_pairs = pd.read_table(
    "../../Puck_180819_12/hotspot/hotspot_pairs_shufled_z.txt.gz", index_col=0
)

hs_pairs_z = hs_pairs.values.ravel()
hs_pairs_p = norm.sf(hs_pairs_z)

plt.figure(figsize=(4.5, 4.5))
z_qq_plot(hs_pairs_z)
plt.xlabel('Theoretical Quantiles')
plt.ylabel('Sample Quantiles')
plt.title('Hotspot Pairwise Statistic\nQ-Q Plot')
plt.subplots_adjust(left=.2, bottom=.2)
# plt.show()
plt.savefig('hotspot_qq_pairs.svg', dpi=300)
