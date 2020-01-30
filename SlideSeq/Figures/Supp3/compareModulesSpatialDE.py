import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json

plt.rcParams['svg.fonttype'] = 'none'

# %% Load module genes

name_map_12_r = json.load(open("../MainFigure/module_name_map.json"))
name_map_12 = {v: k for k, v in name_map_12_r.items()}

modules_hs_df = pd.read_table(
    "../../Puck_180819_12/hotspot/modules.txt", index_col=0
)

modules_hs = {}

for g, c in modules_hs_df.Cluster.iteritems():
    if c == -1: continue
    if c not in modules_hs:
        modules_hs[c] = set()

    modules_hs[c].add(g)

modules_sd_df = pd.read_table(
    "../../Puck_180819_12/spatialDE/histology_windowed_10.txt", index_col=0
)


modules_sd = {}

for row in modules_sd_df.itertuples():
    if row.pattern not in modules_sd:
        modules_sd[row.pattern] = set()

    modules_sd[row.pattern].add(row.g)

# %% Compute overlap

def jaccard(setA, setB):
    setA = {x.lower() for x in setA}
    setB = {x.lower() for x in setB}

    intersection = len(setA & setB)
    union = len(setA | setB)

    return intersection/union


hs_keys = sorted(modules_hs.keys())
sd_keys = sorted(modules_sd.keys())

overlaps = pd.DataFrame(
    0, index=hs_keys, columns=sd_keys
)

for hs, hs_genes in modules_hs.items():
    for sd, sd_genes in modules_sd.items():

        overlaps.loc[hs, sd] = jaccard(hs_genes, sd_genes)

hs_order = [3, 4, 2, 6, 1, 5]
sde_order = [0, 9, 1, 4, 5, 2, 3, 6, 8, 7]

overlaps = overlaps.loc[hs_order, sde_order]

# %% Plot

import re

def make_overlaps_plot():
    plot_data = overlaps.copy()
    plot_data.index = [
        re.match("Module (\d+)", name_map_12[i]).group(1)
        for i in plot_data.index]
    plot_data.columns = ["{}".format(i+1) for i in plot_data.columns]

    sns.heatmap(
        plot_data,
        cmap="Blues",
        vmin=0,
        vmax=0.5,
        cbar_kws=dict(shrink=0.5, ticks=[0, 0.25, 0.5], label="Jaccard\nOverlap"),
    )
    plt.xlabel("SpatialDE Module")
    plt.ylabel("Hotspot Module")

plt.figure()
make_overlaps_plot()
# plt.show()
plt.savefig('Module_Matching_Heatmap.svg')

# %% Plot modules for spatial DE

patterns_sd_df = pd.read_table(
    "../../Puck_180819_12/spatialDE/patterns_windowed_10.txt", index_col=0
)

patterns_hs_df = pd.read_table(
    "../../Puck_180819_12/hotspot/module_scores.txt.gz", index_col=0
)

positions = pd.read_table(
    "../../Puck_180819_12/positions/positions.txt", index_col=0
)

# Z-normalize module scores for plotting

patterns_sd_df = patterns_sd_df \
    .subtract(patterns_sd_df.mean(axis=0), axis=1) \
    .divide(patterns_sd_df.std(axis=0), axis=1)

patterns_hs_df = patterns_hs_df \
    .subtract(patterns_hs_df.mean(axis=0), axis=1) \
    .divide(patterns_hs_df.std(axis=0), axis=1)

data_sd = positions.join(patterns_sd_df).fillna(0)
data_hs = positions.join(patterns_hs_df).fillna(0)

# %%

import matplotlib.colors
import matplotlib.cm

norm_sd = matplotlib.colors.Normalize(vmin=-1, vmax=1, clip=True)
cmap_sd = plt.get_cmap('viridis')
sm = matplotlib.cm.ScalarMappable(norm=norm_sd, cmap=cmap_sd)

# Do them in six columns

pairs = [
    (3, 0),
    (4, 9),
    (2, 1),
    (6, 4),
    (1, 5),
    (5, 3)
]

fig = plt.figure(figsize=(12, 8))
gs0 = fig.add_gridspec(2, 1, left=0.05, right=0.95, hspace=.35,
                       height_ratios=[1, 1.5])
gs01 = gs0[1].subgridspec(1, 3, wspace=0.4)
gs01L = gs01[0].subgridspec(2, 2, hspace=0.45)
gs01M = gs01[1].subgridspec(2, 2, hspace=0.45)
gs01R = gs01[2].subgridspec(2, 2, hspace=0.45)

cax = fig.add_axes(
    [.95, .1, .005, .1]
)

# plot bottom left

kwargs = dict(
    s=.3, rasterized=True, edgecolor='none'
)


def despine():
    ax = plt.gca()
    for x in ax.spines.values():
        x.set_visible(False)


p = 0

for gs in [gs01L, gs01M, gs01R]:
    for r in range(2):

        pair = pairs[p]
        ax = fig.add_subplot(gs[r, 0])
        ax1 = fig.add_subplot(gs[r, 1])

        plt.sca(ax)
        c = sm.to_rgba(data_hs[str(pair[0])])
        plt.scatter(
            data_hs.Comp1, data_hs.Comp2, c=c,
            **kwargs,
        )
        plt.xticks([])
        plt.yticks([])
        despine()
        hs_mod = re.match("Module (\d+)", name_map_12[pair[0]]).group(1)
        plt.title('Hotspot-{}'.format(hs_mod))

        plt.sca(ax1)
        vals = data_sd[str(pair[1])]
        c = sm.to_rgba(vals)
        c[vals == 0] = [.7, .7, .7, 1.]
        plt.scatter(
            data_hs.Comp1, data_hs.Comp2, c=c,
            **kwargs,
        )

        plt.xticks([])
        plt.yticks([])
        despine()
        plt.title('SpatialDE-{}'.format(pair[1]+1))

        p = p+1

cbar = plt.colorbar(sm, cax=cax, ticks=[-1, 0, 1])
plt.sca(cax)
plt.ylabel('Standardized\nModule Scores', size=9)
cax.tick_params(labelsize=8)
# plt.show()
plt.savefig('Module_Matching_Positions.svg', dpi=300)
