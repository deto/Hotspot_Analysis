import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import loompy

plt.rcParams['svg.fonttype'] = 'none'

# For the Monocyte results
ref_loom_file = "../../../data/10x_PBMC_w_proteins/mono/data.loom"
results_dir = "../../Mono_w_protein_downsampling/"


# %% load all the data

ref_loom = loompy.connect(
    ref_loom_file, "r"
)

ref_loom.row_attrs.keys()
gene_info = pd.DataFrame(
    {x: ref_loom.row_attrs[x] for x in ref_loom.row_attrs}
).set_index('EnsID')
ref_loom.close()

rates = ["100", "60", "40", "20", "10"]

hs_z = {}
pearson_r = {}

for rate in rates:

    hs = pd.read_table(
        os.path.join(
            results_dir, "{}/hotspot/hotspot_pairs_z.txt.gz".format(rate)
        ),
        index_col=0
    )

    np.fill_diagonal(hs.values, 0)

    pearson_results = pd.read_table(
        os.path.join(
            results_dir, "{}/hotspot/regular_pairs_lc.txt.gz".format(rate)
        ),
        index_col=0
    )

    np.fill_diagonal(pearson_results.values, 0)

    hs_z[rate] = hs
    pearson_r[rate] = pearson_results


# Also grab the held-out pearson
pearson_r_held_out = pd.read_table(
    os.path.join(results_dir, "test_data/hotspot/regular_pairs_lc.txt.gz"),
    index_col=0
)

np.fill_diagonal(pearson_r_held_out.values, 0)

# %%

from scipy.stats import spearmanr

def compare_spearman(m1, m2):

    # Sometimes indices differ if we lose genes entirely
    common = m1.index.intersection(m2.index)
    m1 = m1.loc[common, common]
    m2 = m2.loc[common, common]

    x = m1.values.ravel()
    y = m2.values.ravel()

    return spearmanr(x, y).correlation

results = []
for rate in rates:
    corr = compare_spearman(pearson_r[rate], pearson_r_held_out)
    results.append(
        ['Pearson', rate, corr]
    )

for rate in rates:
    corr = compare_spearman(hs_z[rate], pearson_r_held_out)
    results.append(
        ['Hotspot', rate, corr]
    )

results = pd.DataFrame(results, columns=['Method', 'Rate', 'Corr'])


# %% Plot

plt.figure(figsize=(6.5, 5))
sns.barplot(
    data=results, x='Rate', y='Corr', hue='Method',
    order=rates, hue_order=['Hotspot', 'Pearson'], alpha=0.9
)
plt.gca().set_axisbelow(True)
plt.grid(axis='y', color='#AAAAAA', ls=(0, (5, 5)), zorder=1)
plt.xlabel('Downsampling Rate (%)')
plt.ylabel('Spearman Correlation vs. Held-Out Data')
plt.savefig('Downsampled_Correlation_Monocytes.svg', dpi=300)


# %% Repeat with the CD4 Results

ref_loom_file = "../../../data/10x_PBMC_w_proteins/cd4/data.loom"
results_dir = "../../CD4_w_protein_downsampling"

# %% load all the data

ref_loom = loompy.connect(
    ref_loom_file, "r"
)

ref_loom.row_attrs.keys()
gene_info = pd.DataFrame(
    {x: ref_loom.row_attrs[x] for x in ref_loom.row_attrs}
).set_index('EnsID')
ref_loom.close()

rates = ["100", "60", "40", "20", "10"]

hs_z = {}
pearson_r = {}

for rate in rates:

    hs = pd.read_table(
        os.path.join(
            results_dir, "{}/hotspot/hotspot_pairs_z.txt.gz".format(rate)
        ),
        index_col=0
    )

    np.fill_diagonal(hs.values, 0)

    pearson_results = pd.read_table(
        os.path.join(
            results_dir, "{}/hotspot/regular_pairs_lc.txt.gz".format(rate)
        ),
        index_col=0
    )

    np.fill_diagonal(pearson_results.values, 0)

    hs_z[rate] = hs
    pearson_r[rate] = pearson_results


# Also grab the held-out pearson
pearson_r_held_out = pd.read_table(
    os.path.join(results_dir, "test_data/hotspot/regular_pairs_lc.txt.gz"),
    index_col=0
)

np.fill_diagonal(pearson_r_held_out.values, 0)

# %%

from scipy.stats import spearmanr

def compare_spearman(m1, m2):

    # Sometimes indices differ if we lose genes entirely
    common = m1.index.intersection(m2.index)
    m1 = m1.loc[common, common]
    m2 = m2.loc[common, common]

    x = m1.values.ravel()
    y = m2.values.ravel()

    return spearmanr(x, y).correlation

results = []
for rate in rates:
    corr = compare_spearman(pearson_r[rate], pearson_r_held_out)
    results.append(
        ['Pearson', rate, corr]
    )

for rate in rates:
    corr = compare_spearman(hs_z[rate], pearson_r_held_out)
    results.append(
        ['Hotspot', rate, corr]
    )

results = pd.DataFrame(results, columns=['Method', 'Rate', 'Corr'])


# %% Plot

plt.figure(figsize=(6.5, 5))
sns.barplot(
    data=results, x='Rate', y='Corr', hue='Method',
    order=rates, hue_order=['Hotspot', 'Pearson'], alpha=0.9
)
plt.gca().set_axisbelow(True)
plt.grid(axis='y', color='#AAAAAA', ls=(0, (5, 5)), zorder=1)
plt.xlabel('Downsampling Rate (%)')
plt.ylabel('Spearman Correlation vs. Held-Out Data')
plt.savefig('Downsampled_Correlation_CD4.svg', dpi=300)
