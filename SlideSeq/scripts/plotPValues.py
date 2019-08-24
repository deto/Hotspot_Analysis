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
    "../Puck_180819_12/hotspot/hotspot.txt", index_col=0
)


plt.figure()
plt.hist(hs_data.Pval, 100)
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


# %% Plot both distributions

is_marker = [
    x.lower() in markers_all for x in hs_data.Symbol
]

plt.figure()
plt.hist(hs_data.Pval, 100, range=(0, 1), density=True, alpha=0.6, label='All')
plt.hist(hs_data.Pval[is_marker], 100, range=(0, 1), density=True, alpha=0.6, label='Markers')
plt.legend()
plt.xlabel('P-Value')
plt.ylabel('PDF')
plt.show()

# %% What about for SDE?
sd_data = pd.read_table(
    "../Puck_180819_12/spatialDE/spatialDE.txt", index_col=0
)

is_marker = [
    x.lower() in markers_all for x in sd_data.index
]

plt.figure()
plt.hist(sd_data.pval, 100, range=(0, 1), density=True, alpha=0.6, label='All')
plt.hist(sd_data.pval[is_marker], 100, range=(0, 1), density=True, alpha=0.6, label='Markers')
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
    density=True
)
y_markers, bins = np.histogram(
    np.log10(hs_data.Pval[is_marker]),
    bins=[float('-inf')] + list(range(-10, 1)),
    density=True
)
y_theo = [10**bins[i] - 10**bins[i-1] for i in range(1, len(bins))]

x = [bins[i]-.5 for i in range(1, len(bins))]

plt.bar(x, y_real, alpha=0.6)
plt.bar(x, y_markers, color=tab10[2], alpha=0.6)
plt.plot(x, y_theo, '--', color=tab10[1])
plt.yscale('log')
plt.show()

# This doesn't look good either
