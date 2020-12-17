import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import re

from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, LinearSegmentedColormap


plt.rcParams['svg.fonttype'] = 'none'

# %% Load HS Results for each puck

puck_analysis_dirs = {
    '9': '../../Puck_180819_9/',
    '10': '../../Puck_180819_10/',
    '11': '../../Puck_180819_11/',
    '12': '../../Puck_180819_12/',
}

name_map_12_r = json.load(open("../MainFigure/module_name_map.json"))
name_map_12 = {v: k for k, v in name_map_12_r.items()}


def load_hs_modules(analysis_dir):
    file_name = "hotspot/modules.txt"

    results = pd.read_table(
        os.path.join(analysis_dir, file_name),
        index_col=0
    )
    return results


def load_hs_modules_nb(analysis_dir):
    file_name = "neg_binom/hotspot/modules.txt"

    results = pd.read_table(
        os.path.join(analysis_dir, file_name),
        index_col=0
    )
    return results


hs_modules = {k: load_hs_modules(v) for k, v in puck_analysis_dirs.items()}
hs_modules_nb = {
    k: load_hs_modules_nb(v) for k, v in puck_analysis_dirs.items()
}

# Map the cluster ids for Puck 12 so they are consistent with other plots
name_map_12_int = {
    k: int(re.search('Module (.*)$', v).group(1))
    for k, v in name_map_12.items()
}

name_map_12_int[-1] = -1

hs_modules['12'].Cluster = hs_modules['12'].Cluster.map(name_map_12_int)
hs_modules_nb['12'].Cluster = hs_modules_nb['12'].Cluster.map(name_map_12_int)


# %%  plot confusion


def plot_pair(x1, ax):

    m1 = hs_modules[x1]
    m2 = hs_modules_nb[x1]
    m1 = m1.rename({'Cluster': 'm1'}, axis=1)
    m2 = m2.rename({'Cluster': 'm2'}, axis=1)

    all_dat = m1.join(m2, how='outer')
    all_dat = all_dat.fillna(-1)

    out = pd.DataFrame(
        0, index=all_dat['m1'].unique(),
        columns=all_dat['m2'].unique()
    )

    for i1 in out.index:
        for i2 in out.columns:
            g1 = all_dat.index[all_dat['m1'] == i1]
            g2 = all_dat.index[all_dat['m2'] == i2]

            intersection = set(g1) & set(g2)
            union = set(g1) | set(g2)

            jacc = len(intersection)/len(union)

            out.loc[i1, i2] = jacc

    out = out.T

    # re-order rows to match as much as possible

    assignments = {}

    out_s = out.drop(-1, axis=0).drop(-1, axis=1)

    while out_s.size > 0:

        iref = (out_s == out_s.values.max()).any(axis=1).idxmax()
        cref = (out_s == out_s.values.max()).any(axis=0).idxmax()

        assignments[iref] = cref

        out_s = out_s.drop(iref, axis=0).drop(cref, axis=1)

    i_order = np.sort(out.index).tolist()
    i_order = i_order[1:] + i_order[0:1]
    out = out.loc[i_order]

    c_order = []
    for x in out.index:
        if x in assignments:
            c_order.append(assignments[x])
    for y in out.columns:

        if y == -1:
            continue

        if y not in c_order:
            c_order.append(y)

    c_order.append(-1)
    out = out.loc[:, c_order]

    # %% Plot result

    def mod_to_name(i):
        if i == -1:
            return 'Unassigned'
        else:
            return 'Module {}'.format(int(i))

    plot_data = out.copy().T
    plot_data.index = [mod_to_name(x) for x in plot_data.index]
    plot_data.columns = [mod_to_name(x) for x in plot_data.columns]

    plt.sca(ax)
    sns.heatmap(plot_data, cmap="Blues", vmin=0, vmax=.75, ax=ax, cbar=False)

    plt.xticks(rotation=45)
    plt.title('Puck_180819_{}'.format(x1))
    plt.xlabel('Neg-Binom')
    plt.ylabel('Bernoulli')


# %%
x1s = ['9', '10', '11', '12']

fig, axs = plt.subplots(2, 2, figsize=(9, 8))
cbar_ax = fig.add_axes([.92, .11, .01, .1])

for x1, ax in zip(x1s, axs.ravel()):
    plot_pair(x1, ax)

norm = Normalize(vmin=0, vmax=0.75)
cmap = LinearSegmentedColormap.from_list(
    name='blues', colors=sns.color_palette("Blues")
)
sm = ScalarMappable(norm, cmap)

plt.colorbar(sm, cax=cbar_ax, label='Jaccard\nOverlap', ticks=[0, .25, .5, .75])

plt.subplots_adjust(
    left=0.15, bottom=0.15, wspace=0.7, hspace=0.7, right=0.86, top=0.95
)
# plt.show()
plt.savefig('ModuleAssignmentComparison.svg', dpi=300)
