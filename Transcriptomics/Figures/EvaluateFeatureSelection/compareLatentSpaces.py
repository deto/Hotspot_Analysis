"""
Compare latent spaces using the protein autocorrelation
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

plt.rcParams['svg.fonttype'] = 'none'

# %% Load data

data = {
    "Threshold": pd.read_table(
        "../../CD4_w_protein/hotspot/hotspot_threshold_proteins.txt",
        index_col=0,
    ),
    "HVG": pd.read_table(
        "../../CD4_w_protein/hotspot/hotspot_hvg_proteins.txt", index_col=0
    ),
    "NBDisp": pd.read_table(
        "../../CD4_w_protein/hotspot/hotspot_danb_proteins.txt", index_col=0
    ),
    "PCA": pd.read_table(
        "../../CD4_w_protein/hotspot/hotspot_pca_proteins.txt", index_col=0
    ),
    "Hotspot": pd.read_table(
        "../../CD4_w_protein/hotspot/hotspot_hs_proteins.txt", index_col=0
    ),
}

# %% Load Protein Data
proteins = pd.read_table("../../../data/10x_PBMC_w_proteins/cd4/ab.txt.gz", index_col=0)

proteins_all = pd.read_table("../../../data/10x_PBMC_w_proteins/ab.txt.gz", index_col=0)
proteins_all = np.log2(proteins_all + 1)


# Looking at correlations between proteins and total protein per cell

# res = {}
# for p in proteins.columns:
#     pc = proteins[p]
#     pr = proteins.sum(axis=1) - pc
#     res[p] = np.corrcoef(pc, pr)[0, 1]
# res = pd.Series(res)
# plt.figure()
# plt.boxplot(res)
# plt.show()

# %% Mixture modeling

from sklearn.mixture import GaussianMixture
from scipy.stats import norm

# p = 'CD8a'
res = []
for p in proteins_all.columns:
    X = proteins_all[p].values.reshape((-1, 1))
    model = GaussianMixture(n_components=2)
    model.fit(X)

    xd = np.linspace(X.min(), X.max(), 1000)
    comp1 = norm.pdf(xd, loc=model.means_[0, 0], scale=model.covariances_[0, 0, 0]**.5)*model.weights_[0]
    comp2 = norm.pdf(xd, loc=model.means_[1, 0], scale=model.covariances_[1, 0, 0]**.5)*model.weights_[1]



    X2 = proteins_all.loc[proteins.index][p].values.reshape((-1, 1))
    pred = model.predict(X2)
    if model.means_[1, 0] < model.means_[0, 0]:
        pred = 1 - pred

    dm = abs(model.means_[0, 0] - model.means_[1, 0])

    p_active = pred.sum() / X2.size

    res.append([p, dm, p_active])

    # plt.figure()
    # plt.hist(proteins_all[p].values, 30, density=True)
    # plt.plot(xd, comp1)
    # plt.plot(xd, comp2)
    # plt.title(p)
    # plt.show()


res = pd.DataFrame(res, columns=['protein', 'dm', 'pactive']).set_index('protein')
res['pactive2'] = res['pactive']
res.loc[res.dm < 2, 'pactive2'] = 0
res = res.sort_values('pactive2')

valid = res.index[res.pactive2 > .01]


# %% Transform data

data_df = pd.concat(
    [v['Z'].rename(k) for k, v in data.items()],
    axis=1
)

# Drop proteins that aren't doing anything
# data_df = data_df.loc[
#     data_df.max(axis=1) > 10
# ]

# Just drop specific ones
# data_df = data_df.drop(['CD14', 'IgG2a', 'IgG2b', 'IgG1'], axis=0)

# Just keep specific proteins
protein_mean = proteins.mean(axis=0)
data_df = data_df.loc[
    #protein_mean.index[protein_mean > 10]
    valid
]

# Convert to delta-s from Threshold
data_df = data_df.subtract(
    data_df['Threshold'], axis=0
)

data_df = data_df.drop('Threshold', axis=1)
data_df.index.name = 'Protein'

data_df = data_df.reset_index() \
    .melt(id_vars=['Protein'], var_name='Method', value_name='DeltaZ')

# %% Plot results

order = ['NBDisp', 'HVG', 'PCA', 'Hotspot']

plt.figure()
sns.boxplot(
    data=data_df,
    x="Method",
    y="DeltaZ",
    whis="range",
    order=order,
    showcaps=False,
    color="#cccccc",
)
sns.stripplot(
    data=data_df, x="Method", y="DeltaZ", color="#000000", order=order, size=3
)
plt.grid(color='#cccccc', lw=.5, linestyle=(0, (5, 5)), axis='y')
plt.gca().set_axisbelow(True)
plt.ylabel('$Z - Z_{Threshold}$')
# plt.show()
plt.savefig('deltaZ.svg')

# %% Compare significance

from scipy.stats import ranksums


for mm1 in data_df.Method.unique():
    for mm2 in data_df.Method.unique():
        if mm1 == mm2: continue

        rh = data_df.loc[(data_df.Method == mm1)] \
            .set_index('Protein') \
            .loc[valid] \
            .DeltaZ

        rx = data_df.loc[(data_df.Method == mm2)] \
            .set_index('Protein') \
            .loc[valid] \
            .DeltaZ

        test_result = ranksums(rh.values, rx.values)
        print('{} vs {}'.format(mm1, mm2), ':', test_result.statistic, test_result.pvalue)
