"""
Let's compare the computation of pair-wise Z-scores for both
this replicate and the other dataset with the proteins
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import loompy
from hotspot.local_stats_pairs import create_centered_counts

# %% Load gene table
loom_file = "../../data/10x_PBMC/data.loom"
with loompy.connect(loom_file, mode='r') as ds:
    e = ds.ra['EnsID'][:]
    s = ds.ra['Symbol'][:]
    ens_map = {e[i]: s[i] for i in range(len(e))}
    sym_map = {s[i]: e[i] for i in range(len(e))}

# %% Load data

pairs1 = pd.read_table("../hotspot/hotspot_pairs_lc.txt.gz", index_col=0)
pairs2 = pd.read_table("../../CD4_w_protein/hotspot/hotspot_pairs_lc.txt.gz", index_col=0)

common = pairs1.index.intersection(pairs2.index)

print(pairs1.shape[0], pairs2.shape[0], len(common))


pairs1 = pairs1.loc[common, common]
pairs2 = pairs2.loc[common, common]

# %%

def convert_to_long(pairs):
    pairs.index.name = "EnsID"
    pairs = pairs.reset_index()
    pairs.index.name = None
    pairs.columns.name = None
    pairs = pairs.melt(id_vars="EnsID")
    pairs.columns = ["Gene1", "Gene2", "Z"]
    pairs = pairs.set_index(["Gene1", "Gene2"])
    return pairs

pairs1 = convert_to_long(pairs1)
pairs1.columns = ["Z1"]

pairs2 = convert_to_long(pairs2)
pairs2.columns = ["Z2"]

data = pairs1.join(pairs2)


# %% Plot the correspondence

plt.figure()
plt.plot(data.Z1, data.Z2, 'o', ms=2)
plt.show()



# %% What if we look at regular correlations?

pairs1_c = pd.read_table("../hotspot/regular_pairs_lc.txt.gz", index_col=0)
pairs2_c = pd.read_table("../../CD4_w_protein/hotspot/regular_pairs_lc.txt.gz", index_col=0)

common = pairs1_c.index.intersection(pairs2_c.index)

print(pairs1_c.shape[0], pairs2_c.shape[0], len(common))


pairs1_c = pairs1_c.loc[common, common]
pairs2_c = pairs2_c.loc[common, common]

# %%

pairs1_c = convert_to_long(pairs1_c)
pairs1_c.columns = ["Z1_c"]

pairs2_c = convert_to_long(pairs2_c)
pairs2_c.columns = ["Z2_c"]

data_c = pairs1_c.join(pairs2_c)
data = data.join(data_c)


# %% Plot the correspondence

plt.figure()
plt.plot(data.Z1_c, data.Z2_c, 'o', ms=2)
plt.show()

# %% Compare between?

plt.figure()
plt.plot(data.Z2_c, data.Z2, 'o', ms=2, alpha=0.1, mew=0)
plt.plot([-1, 1], [-1, 1], '--', lw=1)
plt.show()

# %% How can we model the error? MA Plot?

x = data.Z1
y = data.Z2

a = (x+y)/2
m = np.abs(x-y)

fig, axs = plt.subplots(1, 2, figsize=(9, 5), sharex=True, sharey=True)
plt.sca(axs[0])
plt.plot(a, m, 'o', ms=1)
print(m.mean())

# With actual correlations

x = data.Z1_c
y = data.Z2_c

a = (x+y)/2
m = np.abs(x-y)

plt.sca(axs[1])
plt.plot(a, m, 'o', ms=1)
plt.xlim(-1, 1)
plt.ylim(-.1, 1)
print(m.mean())
plt.show()


# %% How do clustermaps look?

pairs_hs = pd.read_table("../hotspot/hotspot_pairs_lc.txt.gz", index_col=0)
pairs_hs_z = pd.read_table("../hotspot/hotspot_pairs_z.txt.gz", index_col=0)
pairs_reg = pd.read_table("../hotspot/regular_pairs_lc.txt.gz", index_col=0)
ii = pairs_hs.index
pairs_reg = pairs_reg.loc[ii, ii]

def symbol_row_cols(df):
    df.index = df.index.map(ens_map)
    df.columns = df.columns.map(ens_map)
    return df

pairs_hs = symbol_row_cols(pairs_hs)
pairs_hs_z = symbol_row_cols(pairs_hs_z)
pairs_reg = symbol_row_cols(pairs_reg)

# %% 

from scipy.cluster.hierarchy import linkage, leaves_list
from scipy.spatial.distance import squareform

dd = 2 - pairs_hs.values
np.fill_diagonal(dd, 0)
Z1 = linkage(squareform(dd), method='average')

sns.clustermap(pairs_hs, vmin=-.2, vmax=.2, cmap="RdBu_r",
               col_linkage=Z1, row_linkage=Z1,
               yticklabels=True)
plt.show()

# %% Other thing

dd = 1 - pairs_reg.values
np.fill_diagonal(dd, 0)
# They think it's not symmetric?  Odd as max difference is ~4e-16.  Cludge to fix
dd = (dd + dd.T)/2
Z2 = linkage(squareform(dd), method='average')

sns.clustermap(pairs_reg, vmin=-.2, vmax=.2, cmap="RdBu_r",
               col_linkage=Z2, row_linkage=Z2,
               yticklabels=True)
plt.show()

# %% Z values
from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
results = pairs_hs_z
offset = pairs_hs_z.values.max()
dd = offset - pairs_hs_z.values
np.fill_diagonal(dd, 0)
ZZ = linkage(squareform(dd), method='average')

sns.clustermap(pairs_hs_z, vmin=-20, vmax=20, cmap="RdBu_r",
               col_linkage=ZZ, row_linkage=ZZ,
               yticklabels=True)
plt.show()

# %% How well do they correspond between the two?
ll = leaves_list(Z2)


plt.figure()
sns.heatmap(pairs_hs.iloc[ll, ll],
            vmin=-.2, vmax=.2, cmap="RdBu_r")
plt.show()

dendrogram(ZZ, labels=pairs_hs_z.index)
plt.hlines(offset, 0, ZZ.shape[0])
plt.show()

dendrogram(Z2, labels=pairs_reg.index, leaf_font_size=10)
plt.show()


# %% How do these correspond so well?  Why doesn't the noise matter?
with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

latent_file = "../scvi/hvg/latent.txt.gz"
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)
num_umi = pd.Series(num_umi, index=barcodes)
latent = pd.read_table(latent_file, index_col=0)
counts = counts.loc[:, latent.index]
num_umi = num_umi[latent.index]

counts = counts.loc[
    counts.sum(axis=1) > 3
]

# %%

x = counts.loc[sym_map['FOXP3']].values
y = counts.loc[sym_map['IL2RA']].values

np.corrcoef(np.vstack((x, y)))

plt.figure()
plt.plot(x + np.random.randn(len(x))*.5,
         y + np.random.randn(len(x))*.5,
         'o', ms=2)
plt.show()

n = ((x-x.mean())*(y-y.mean())).sum() / (len(x))
d = x.std() * y.std()

n/d


# Have to do this because data_slideseq makes it a numpy array

# Align to latent space

# need counts, latent, and num_umi

import hotspot
n_neighbors = 30
hs = hotspot.Hotspot(counts, latent, num_umi)

hs.create_knn_graph(
    weighted_graph=False, n_neighbors=n_neighbors, neighborhood_factor=3
)

from numba import njit


@njit
def neighborhood_impute(vals, neighbors, weights):

    N_CELLS = len(vals)
    out = np.zeros(N_CELLS)
    out_denom = np.zeros(N_CELLS)

    for i in range(N_CELLS):
        for k in range(neighbors.shape[1]):
            j = neighbors[i, k]
            wij = weights[i, k]

            out[i] += vals[j] * wij
            out_denom[i] += wij

            out[j] += vals[i] * wij
            out_denom[j] += wij

    for i in range(N_CELLS):
        out[i] = out[i] / out_denom[i]

    return out


x_impute = neighborhood_impute(x, hs.neighbors.values, hs.weights.values)
y_impute = neighborhood_impute(y, hs.neighbors.values, hs.weights.values)

plt.figure()
plt.plot(x_impute, y_impute, 'o', ms=2)
plt.show()

np.corrcoef(np.vstack((x_impute, y_impute)))

(np.corrcoef(np.vstack((x_impute, y)))[0, 1] + np.corrcoef(np.vstack((y_impute, x)))[0, 1])/2

# What about centered?

c_counts = create_centered_counts(
    counts.loc[[sym_map['FOXP3'], sym_map['IL2RA']]].values , 'danb', num_umi.values)

x_c = c_counts[0, :]
y_c = c_counts[1, :]

plt.figure()
plt.plot(x_c, y_c, 'o', ms=2)
plt.show()

np.corrcoef(np.vstack((x_c, y_c)))

x_c_impute = neighborhood_impute(x_c, hs.neighbors.values, hs.weights.values)
y_c_impute = neighborhood_impute(y_c, hs.neighbors.values, hs.weights.values)

plt.figure()
plt.plot(x_c_impute, y_c_impute, 'o', ms=2)
plt.show()

np.corrcoef(np.vstack((x_c_impute, y_c_impute)))

hs_results = hs.compute_hotspot(model='danb', jobs=20, centered=True)

hs_results.loc[sym_map['FOXP3']]

hs_genes = [sym_map['FOXP3'], sym_map['IL2RA'], sym_map['FOXP3']]
lc, lcz = hs.compute_modules(hs_genes, model='danb', centered=True, jobs=1)

# %% Ok, let's look at a really weird synthetic dataset

N_CELLS = 5000
latent = pd.DataFrame(np.linspace(0, 1, N_CELLS).reshape((-1, 1)))
num_umi = pd.Series(np.ones(N_CELLS)*1000)
counts = pd.DataFrame({
    'Gene1': latent.values.ravel()**2 + np.random.binomial(2, .5, N_CELLS),
    'Gene2': latent.values.ravel() + np.random.binomial(2, .5, N_CELLS),
    'Gene3': -1*latent.values.ravel()+1,
}).T

n_neighbors = 30
hs = hotspot.Hotspot(counts, latent, num_umi)

hs.create_knn_graph(
    weighted_graph=False, n_neighbors=n_neighbors, neighborhood_factor=3
)
hs_results = hs.compute_hotspot(model='danb', jobs=20, centered=True)

lc, lcz = hs.compute_modules(counts.index, model='danb', centered=True, jobs=1)

# This looks fine actually, meaning our limits aren't incorrect?

from scipy.stats import pearsonr
from scipy.stats import t as tdist
pearsonr(x_c, y_c)
pearsonr(x_c_impute, y_c_impute)


np.corrcoef(counts.loc['Gene1'].values, counts.loc['Gene2'].values)

g1 = neighborhood_impute(counts.loc['Gene1'].values, hs.neighbors.values, hs.weights.values)
g2 = neighborhood_impute(counts.loc['Gene2'].values, hs.neighbors.values, hs.weights.values)

np.corrcoef(g1, g2)

# Maybe the local correlation is bounded by the local autocorrelation?

lc_x = pd.read_table("../hotspot/hotspot_pairs_lc.txt.gz", index_col=0)
hsr = pd.read_table("../hotspot/hotspot_hvg.txt", index_col=0)

x = hsr.loc[lc_x.index].C.values
max_vals = (x.reshape((-1, 1)) + x.reshape((1, -1)))/2

plt.figure()
plt.plot(lc_x.values.ravel(), max_vals.ravel(), 'o', ms=2)
plt.plot([-1, 1], [-1, 1], '--', lw=1)
plt.plot([-1, 1], [1, -1], '--', lw=1)
plt.show()

lc_scaled = lc_x / max_vals
lc_scaled = symbol_row_cols(lc_scaled)

import seaborn as sns
sns.clustermap(lc_scaled, vmin=-1, vmax=1, method='average', cmap='RdBu_r')
plt.show()

from scipy.cluster.hierarchy import linkage
from scipy.spatial.distance import squareform

dd = lc_scaled.values
np.fill_diagonal(dd, 1.3)
dd = 1.3 - dd
Z = linkage(squareform(dd), method='average')
sns.clustermap(lc_scaled, vmin=-1, vmax=1, row_linkage=Z, col_linkage=Z, cmap='RdBu_r')
plt.show()

a = pd.read_table("../hotspot/hotspot_pairs_lc.txt.gz", index_col=0)
b = pd.read_table("../hotspot/hotspot_pairs_z.txt.gz", index_col=0)

plt.plot(
    lc_scaled.values.ravel(),
    b.values.ravel(),
    'o', ms=2
)
plt.show()


from hotspot.local_stats_pairs import local_cov_pair
from hotspot.local_stats import local_cov_weights
from hotspot.knn import compute_node_degree

c_counts = create_centered_counts(
    counts.loc[[sym_map['FOXP3'], sym_map['IL2RA']]].values , 'danb', num_umi.values)

x_c = c_counts[0, :]
y_c = c_counts[1, :]

g = local_cov_pair(x_c, y_c, hs.neighbors.values, hs.weights.values)
g = local_cov_pair(x_c+1, y_c, hs.neighbors.values, hs.weights.values)
local_cov_pair(x_c, x_c, hs.neighbors.values, hs.weights.values)

D = compute_node_degree(hs.neighbors.values, hs.weights.values)

def norm_vals(vals, D):
    mu = (vals*D).sum() / D.sum()
    vals_c = vals - mu

    std = ((vals_c**2*D).sum() / D.sum())**.5
    vals_n = vals_c / std
    return vals_n

x_cn = norm_vals(x_c, D)
y_cn = norm_vals(y_c, D)

g0 = local_cov_pair(x_cn, x_cn, hs.neighbors.values, hs.weights.values)
g0 = local_cov_weights(x_cn, hs.neighbors.values, hs.weights.values)

from tqdm import tqdm

gs = []
for i in tqdm(range(1000)):
    noise = np.random.randn(len(x_cn))*1e-1
    x2_cn = norm_vals(x_cn + noise, D)
    gs.append(
        local_cov_pair(x_cn, x2_cn, hs.neighbors.values, hs.weights.values)
    )
gs = np.array(gs)


@njit
def local_cov_pair_d(x, y, neighbors, weights):
    """Test statistic for local pair-wise autocorrelation"""
    out = 0

    for i in range(len(x)):
        for k in range(neighbors.shape[1]):

            j = neighbors[i, k]
            w_ij = weights[i, k]

            xi = x[i]
            yj = y[j]

            out += w_ij*xi*yj

    return out

g0 = local_cov_weights(x_cn, hs.neighbors.values, hs.weights.values)

from tqdm import tqdm

gs = []
lcs = []
for i in tqdm(range(1000)):
    noise = np.random.randn(len(x_cn))*1e-1
    x2_cn = norm_vals(x_cn + noise, D)
    lcs.append(
        local_cov_weights(x2_cn, hs.neighbors.values, hs.weights.values)
    )
    gs.append(
        local_cov_pair_d(x_cn, x2_cn, hs.neighbors.values, hs.weights.values)
    )
gs = np.array(gs)
lcs = np.array(lcs)

plt.figure()
lmax = np.array([max(g0, x) for x in lcs])
plt.plot(gs)
plt.plot(lmax)
plt.show()

plt.hist(gs - lmax)
plt.show()

(lmax - gs).min()

"""
This whole thing has been a giant mess.

Ok, so originally I tried to compare the LC values with regular correlations to
try to show they were better.

In fact, they were consistently of a lower magnitude.
