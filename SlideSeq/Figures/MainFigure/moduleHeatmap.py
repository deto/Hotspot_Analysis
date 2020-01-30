import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import json
from scipy.cluster.hierarchy import linkage, leaves_list
plt.rcParams['svg.fonttype'] = 'none'

# %% Load data

results_hs = pd.read_table(
    "../../Puck_180819_12/hotspot/hotspot.txt",
    index_col=0
)

clusters = pd.read_table(
    "../../Puck_180819_12/hotspot/modules.txt", index_col=0
).Cluster

Z = pd.read_table(
    "../../Puck_180819_12/hotspot/linkage.txt", header=None
).values

results_z = pd.read_table(
    "../../Puck_180819_12/hotspot/hotspot_pairs_z.txt.gz", index_col=0
)

# %%


ii = leaves_list(Z)

results_z = results_z.iloc[ii, ii]

# %% Define a module: name map and save for other plotting scripts

name_map = {}
next_mod = 1
for i in results_z.index:
    mod = clusters[i]
    if mod in name_map:
        continue
    if mod == -1:
        continue
    else:
        name_map[int(mod)] = "Module {}".format(next_mod)
        next_mod = next_mod + 1
name_map[-1] = "Module 0"

name_map_r = {v: k for k, v in name_map.items()}
json.dump(name_map_r, open("module_name_map.json", "w"))

# %%

colors = sns.color_palette("muted")
colormap = {i: colors[i % 10] if i != -1 else (1.0, 1.0, 1.0) for i in clusters.unique()}
colormap_r = {v: k for k, v in colormap.items()}
row_colors1 = pd.Series(
    [colormap[i] for i in clusters],
    index=clusters.index,
)

row_colors = pd.DataFrame({"Module": row_colors1})
row_colors = row_colors.loc[results_z.index]

cmap = sns.diverging_palette(230, 20, sep=20, s=85, l=50, as_cmap=True)

cm = sns.clustermap(
    results_z,
    row_cluster=False,
    col_cluster=False,
    vmin=-8,
    vmax=8,
    cmap=cmap,
    xticklabels=False,
    yticklabels=False,
    row_colors=row_colors,
    rasterized=True,
    cbar_kws=dict(
        label="Z-Score",
        ticks=[-8, 0, 8],
    ),
)

fig = plt.gcf()
fig.patch.set_visible(False)
plt.sca(cm.ax_heatmap)
plt.ylabel("")
plt.xlabel("")

# Assuming row_colors is ordered, label the groups
plt.sca(cm.ax_row_colors)
for i in row_colors.Module.unique():
    if i == (1.0, 1.0, 1.0):
        continue
    mid = 0
    count = 0
    cl = colormap_r[i]
    for j, x in enumerate(row_colors.Module.values):
        if x == i:
            mid += j
            count += 1

    mid = mid/count

    plt.text(-.4, mid, name_map[cl],
             horizontalalignment='right',
             verticalalignment='center')

plt.savefig('Module_Heatmap.svg', dpi=300)
# plt.show()

