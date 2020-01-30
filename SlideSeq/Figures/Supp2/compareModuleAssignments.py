import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

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


hs_modules = {k: load_hs_modules(v) for k, v in puck_analysis_dirs.items()}

# %%  plot confusion


def plot_pair(x1, x2, ax):

    m1 = hs_modules[x1]
    m2 = hs_modules[x2]
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
        if y == -1: continue
        if y not in c_order:
            c_order.append(y)

    c_order.append(-1)
    out = out.loc[:, c_order]

    # %% Plot result

    def mod_to_name(i):
        if i == -1:
            return 'Unassigned'
        else:
            return 'Module {}'.format(int(i+1))

    plot_data = out.copy().T
    plot_data.index = [mod_to_name(x) for x in plot_data.index]
    plot_data.columns = [name_map_12[x] for x in plot_data.columns]

    plt.sca(ax)
    sns.heatmap(plot_data, cmap="Blues", vmin=0, vmax=.75, ax=ax, cbar=False)

    plt.xticks(rotation=45)
    plt.xlabel('Puck_180819_{}'.format(x2))
    plt.ylabel('Puck_180819_{}'.format(x1))

# %%

x2 = '12'

x1s = ['9', '10', '11']

fig, axs = plt.subplots(3, 1, figsize=(6, 10))
cbar_ax = fig.add_axes([.85, .11, .01, .1])

for x1, ax in zip(x1s, axs.ravel()):
    plot_pair(x1, x2, ax)

from matplotlib.cm import ScalarMappable
from matplotlib.colors import Normalize, LinearSegmentedColormap

norm = Normalize(vmin=0, vmax=0.75)
cmap = LinearSegmentedColormap.from_list(name='blues', colors=sns.color_palette("Blues"))
sm = ScalarMappable(norm, cmap)

plt.colorbar(sm, cax=cbar_ax, label='Jaccard\nOverlap', ticks=[0, .25, .5, .75])

plt.subplots_adjust(left=0.3, bottom=0.1, hspace=0.7, right=0.80, top=0.98)
# plt.show()
plt.savefig('ModuleAssignmentComparison.svg')
