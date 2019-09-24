"""
Plotting p-value distributions for hotspot and spatialDE

Both the general distributions and the distribution of 'positive'
genes taken from the DropViz data
"""
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
plt.rcParams['svg.fonttype'] = 'none'

tab10 = plt.get_cmap('tab10').colors

hs_data = pd.read_table(
    "../Puck_180819_12/hotspot/hotspot_300.txt", index_col=0
)

hs_data_shuffled = pd.read_table(
    "../Puck_180819_12/hotspot/hotspot_shuffled.txt", index_col=0
)


plt.figure()
plt.hist(hs_data.Pval, 100, alpha=0.5, label='Real')
plt.hist(hs_data_shuffled.Pval, 100, alpha=0.5, label='Shuffled')
plt.legend()
plt.show()

# %% Load positives:
pos_data = pd.read_table(
    "../DropViz/edger_markers_1vAll.txt"
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
plt.show()

# %% What about for SDE?
sd_data = pd.read_table(
    "../Puck_180819_12/spatialDE/spatialDE_fixed.txt", index_col=0
)

sd_data_shuffled = pd.read_table(
    "../Puck_180819_12/spatialDE/spatialDE_shuffled.txt", index_col=0
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
plt.show()

# %% Need log-value PDFs instead?

is_marker = [
    x.lower() in markers_all for x in hs_data.Symbol
]

y_real, bins = np.histogram(
    np.log10(hs_data.Pval),
    bins=[float('-inf')] + list(range(-10, 1)),
    density=False
)
y_markers, bins = np.histogram(
    np.log10(hs_data.Pval[is_marker]),
    bins=[float('-inf')] + list(range(-10, 1)),
    density=False
)
y_shuffled, bins = np.histogram(
    np.log10(hs_data_shuffled.Pval),
    bins=[float('-inf')] + list(range(-10, 1)),
    density=False
)
y_theo = np.array([10**bins[i] - 10**bins[i-1] for i in range(1, len(bins))])
y_theo = y_theo * hs_data.shape[0]

x = [bins[i]-.5 for i in range(1, len(bins))]

plt.bar(x, y_real, alpha=0.6, label='All')
#plt.bar(x, y_markers, color=tab10[2], alpha=0.6, label='Markers')
plt.bar(x, y_shuffled, color=tab10[3], alpha=0.6, label='Shuffled')
plt.plot(x, y_theo, '--', color=tab10[1], label='Theoretical')
plt.yscale('log')
plt.legend()
plt.show()

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

    plt.plot(np.log10(expected), np.log10(sorted_pvals), 'o-', ms=3)
    plt.plot(np.log10(expected), np.log10(expected), '--', linewidth=1)


fig, axs = plt.subplots(1, 2, figsize=(9, 5))
plt.sca(axs[0])
log_qq_plot(hs_data_shuffled.Pval.values)
plt.xlabel('Expected')
plt.ylabel('Actual')
plt.title('Hotspot')

plt.sca(axs[1])
log_qq_plot(sd_data_shuffled.pval.values)
plt.xlabel('Expected')
plt.ylabel('Actual')
plt.title('spatialDE')
plt.suptitle('Q-Q plot for shuffled data')

plt.show()


# %%

import loompy
loom_file = "../../data/SlideSeq/Puck_180819_12/data.loom"
ds = loompy.connect(loom_file, mode="r")
counts = ds[:, :]
genes = ds.row_attrs["Symbol"]
num_umi = pd.Series(ds.col_attrs['NumUmi'])
ds.close()

gene_detects = (counts > 0).sum(axis=1)
gene_detects = pd.Series(gene_detects, index=genes)

gene_detects = gene_detects.loc[hs_data_shuffled.index]
plt.figure()
plt.plot(hs_data_shuffled.Z, gene_detects, 'o', ms=2)
plt.show()


n_genes = pd.Series((counts > 0).sum(axis=0))

import hotspot
with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

latent = pd.read_table("../Puck_180819_12/positions/positions_shuffled.txt", index_col=0)

# Have to do this because data_slideseq makes it a numpy array
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)
num_umi = pd.Series(num_umi, index=barcodes)

# Align to latent space
counts = counts.loc[:, latent.index]
num_umi = num_umi[latent.index]

import hotspot.sim_data
shuffled_exp = hotspot.sim_data.generate_permutation_null(
    counts.loc['mt-Rnr1'], 2000
)
shuffled_exp = pd.DataFrame(shuffled_exp, columns=counts.columns)

const_umi = pd.Series(np.ones_like(num_umi)*1000, index=latent.index)
hs = hotspot.Hotspot(shuffled_exp, latent, const_umi)

hs.create_knn_graph(
    weighted_graph=False, n_neighbors=300, neighborhood_factor=3
)

results = hs.compute_hotspot(model='bernoulli', jobs=20, centered=True)

from hotspot.bernoulli_model import find_gene_p

# %% Plot the fit of the bernoulli model for a single gene

gene = 'mt-Rnr1'
gene = 'Malat1'
gene = 'Car8'
gene = 'Gpr39'
gene = 'Cox4i1'  # This one fits great
gene = 'Gpm6b'  # This one does the worst
data = pd.DataFrame({
    'gene': counts.loc[gene] > 0,
    'umi': num_umi
})
data['bin'] = pd.qcut(data.umi, 30)

rr = data.groupby('bin')['gene'].mean()
i_mid = np.array([x.mid for x in rr.index])

gene_p = find_gene_p(data['umi'].values, data['gene'].sum())
detect_p_x = np.logspace(np.log10(rr.index[0].left),
                         np.log10(rr.index[-1].right)
                         , 300)
detect_p = 1-(1-gene_p)**detect_p_x

plt.figure()
plt.plot(i_mid, rr.values, 'o')
plt.plot(detect_p_x, detect_p)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('UMIs/Barcode')
plt.ylabel('P(detect)')
plt.title('Actual vs. Predicted P(detect) for {}'.format(gene))
plt.show()


# %% Evaluate the likelihood on a per-gene basis

from numba import jit

@jit(nopython=True)
def log_likelihood(gene_counts, umi_counts):

    tot = 0
    N = gene_counts.size

    for i in range(N):
        if gene_counts[i] > 0:
            tot += 1

    gene_p = find_gene_p(umi_counts, tot)
    detect_p = 1-(1-gene_p)**umi_counts

    ll = 0

    for i in range(gene_counts.size):
        if gene_counts[i] > 0:
            ll += np.log(detect_p[i])
        else:
            ll += np.log(1-detect_p[i])

    return ll


@jit(nopython=True)
def log_likelihood_const_umi(gene_counts):

    tot = 0
    N = gene_counts.size

    for i in range(N):
        if gene_counts[i] > 0:
            tot += 1

    detect_p = tot / N

    ll = 0

    for i in range(gene_counts.size):
        if gene_counts[i] > 0:
            ll += np.log(detect_p)
        else:
            ll += np.log(1-detect_p)

    return ll


@jit(nopython=True)
def log_likelihood_p(gene_counts, detect_p):

    ll = 0

    for i in range(gene_counts.size):
        if gene_counts[i] > 0:
            ll += np.log(detect_p[i])
        else:
            ll += np.log(1-detect_p[i])

    return ll


umi_counts = num_umi.values
gene = 'Fth1'

log_likelihood(counts.loc[gene].values, umi_counts)
log_likelihood_const_umi(counts.loc[gene].values)

from tqdm import tqdm
LLs = []
LLCs = []

cc_num = counts.values

for i in tqdm(range(cc_num.shape[0])):
    LLs.append(log_likelihood(cc_num[i], umi_counts))
    LLCs.append(log_likelihood_const_umi(cc_num[i]))

LLs = pd.Series(LLs, counts.index)
LLCs = pd.Series(LLCs, counts.index)


# %% Plot the change in LL

plt.figure()
plt.hist(LLs-LLCs, 200)
plt.xlabel('$\Delta LL$')
plt.ylabel('# of Genes')
plt.yscale('log')
plt.title('Change in Likelihood (UMI-correcting vs. non-correcting)')
plt.show()


# %% Now, let's plot the change in LL vs. the shuffled Z score

plt.figure()
plt.plot(
    hs_data_shuffled.Z[LLs.index],
    LLs-LLCs, 
    'o', ms=2
)
plt.show()


# %% Should we just be doing logistic regression?

from sklearn.linear_model import LogisticRegression

gene = 'Gpm6b'
gene_i = counts.index.tolist().index(gene)

umi_counts_vec = num_umi.values.reshape((-1, 1))
y = (cc_num[gene_i] > 0).astype('int')



data = pd.DataFrame({
    'gene': counts.loc[gene] > 0,
    'umi': num_umi
})
data['bin'] = pd.qcut(data.umi, 30)

rr = data.groupby('bin')['gene'].mean()
i_mid = np.array([x.mid for x in rr.index])

gene_p = find_gene_p(data['umi'].values, data['gene'].sum())
detect_p_x = np.logspace(np.log10(rr.index[0].left),
                         np.log10(rr.index[-1].right)
                         , 300)
detect_p = 1-(1-gene_p)**detect_p_x

import time
a = time.time()
model = LogisticRegression(penalty='none', solver='lbfgs', max_iter=1000)
model.fit(umi_counts_vec, y)
p_pred = model.predict_proba(detect_p_x.reshape((-1, 1)))[:, 1]
b = time.time()
print("{:.2f} Elapsed".format(b-a))


plt.figure()
plt.plot(i_mid, rr.values, 'o', label='Actual')
plt.plot(detect_p_x, detect_p, label='Theoretical')
plt.plot(detect_p_x, p_pred, label='LogisticRegression')
plt.xscale('log')
plt.yscale('log')
plt.xlabel('UMIs/Barcode')
plt.ylabel('P(detect)')
plt.title('Actual vs. Predicted P(detect) for {}'.format(gene))
plt.legend()
plt.show()

l1 = log_likelihood(counts.loc[gene].values, umi_counts)
l2 = log_likelihood_const_umi(counts.loc[gene].values)
detect_p = model.predict_proba(umi_counts_vec)[:, 1]
l3 = log_likelihood_p(counts.loc[gene].values, detect_p)
