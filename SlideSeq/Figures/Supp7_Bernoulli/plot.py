import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

plt.rcParams['svg.fonttype'] = 'none'


# %% Load files

N_BINS_BERNOULLI = [5, 10, 20, 30, 50, 75, 100]

results = {}
for n_bins in N_BINS_BERNOULLI:

    path = "../../Puck_180819_12/bernoulli_bins/{}_bins/hotspot.txt".format(n_bins)

    hs_result = pd.read_table(path, index_col=0)
    results[n_bins] = hs_result


# %% Compute overlap of top 1000 genes with nbins=30 result
def get_top_n(result):
    return result.sort_values("Z", ascending=False).index[0:1000]


REF_BINS = 30
ref = get_top_n(results[REF_BINS])

avg_rank = {}
for n_bin in results:
    top_genes = get_top_n(results[n_bin])
    avg_rank[n_bin] = results[n_bin]['Z'].rank(ascending=False).loc[ref].mean()

avg_rank = pd.Series(avg_rank)

# %% Make the plot

plt.figure()
ax = plt.gca()
plt.plot(avg_rank.index, avg_rank.values, 'o-')
plt.xscale('log')
plt.xticks(avg_rank.index)
ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
ax.xaxis.set_ticks([], minor=True)
plt.xlabel('# of Neighbors')
plt.ylabel('Average Rank of\nHS(bins={}) Genes'.format(REF_BINS))
plt.gca().set_axisbelow(True)
plt.grid(color='#CCCCCC', ls=(0, (5, 5)), lw=.5)
plt.subplots_adjust(left=0.15)
plt.savefig('AverageRank_bernoulli_nbin_sensitivity.svg', dpi=300)
# plt.show()
