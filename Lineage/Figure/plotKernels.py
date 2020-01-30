import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

from scipy.cluster.hierarchy import linkage, leaves_list

plt.rcParams['svg.fonttype'] = 'none'

# %% Load kernel data and modules

module_file = "../Embryo3/hotspot/modules_lineage_tree.txt"
modules = pd.read_table(module_file, index_col=0)
modules.columns = ['Module']

hs_results_file = "../Embryo3/hotspot/hotspot_lineage_tree.txt"
hs_results = pd.read_table(hs_results_file, index_col=0)
hs_results_mod = hs_results.join(modules)

kernels = pd.read_excel("../../data/Lineage/GSE122187_CellStateKernels.xls", skiprows=1).set_index('Gene')
kernels = kernels.loc[
    :, kernels.columns.sort_values()
]

kernels_z = kernels \
    .subtract(kernels.mean(axis=1), axis=0) \
    .divide(kernels.std(axis=1), axis=0)

kernel_name_map = {
    0: "trophoblast stem cells",
    1: "neural ectoderm anterior",
    2: "primitive streak late",
    3: "anterior primitive streak",
    4: "primitive/definitive endoderm",
    5: "allantois",
    6: "secondary heart field/splanchnic lateral plate",
    7: "gut endoderm",
    8: "ectoderm early 1",
    9: "primitive blood early",
    10: "preplacodal ectoderm",
    11: "neural ectoderm posterior",
    12: "posterior lateral plate mesoderm",
    13: "hematopoietic/endothelial progenitors",
    14: "parietal endoderm",
    15: "amnion mesoderm early",
    16: "surface ectoderm",
    17: "epiblast",
    18: "somites",
    19: "ectoderm early 2",
    20: "splanchnic lateral plate/anterior paraxial mesoderm",
    21: "primitive heart tube",
    22: "primitive blood late",
    23: "notochord",
    24: "fore/midbrain",
    25: "distal extraembryonic ectoderm",
    26: "neuromesodermal progenitor early",
    27: "primordial germ cells",
    28: "differentiated trophoblasts",
    29: "visceral endoderm early",
    30: "presomitic mesoderm",
    31: "neuromesodermal progenitor late",
    32: "angioblasts",
    33: "neural crest",
    34: "pharyngeal arch mesoderm",
    35: "similar to neural crest",
    36: "primitive blood progenitors",
    37: "primitive streak early",
    38: "node",
    39: "future spinal cord",
    40: "visceral endoderm late",
    41: "amnion mesoderm late",
}

kernel_name_map = {k: v.capitalize() for k, v in kernel_name_map.items()}

kernels_z.columns = [kernel_name_map[x] for x in kernels_z.columns]

# %% Order module genes

layers = {
    "Epiblast": [17],
    "PGC": [27],
    "Embryonic Mesoderm": [2, 3, 12, 13, 18, 20, 23, 26, 30, 31, 34, 37, 38, 6, 21, 32, 9, 22, 36],
    "Extraembryonic Mesoderm": [15, 41, 5],
    "Endoderm": [14, 7, 4, 29, 40],
    "Embryonic Ectoderm": [8, 19, 1, 11, 24, 39, 33, 35, 10, 16],
    "Extraembryonic Ectoderm": [0, 25, 28],
}

all_mod_genes = hs_results_mod.loc[hs_results_mod.Module >= 0].Symbol

k_sub = kernels_z.loc[pd.Index(all_mod_genes) & kernels.index]
ii_auto = leaves_list(linkage(k_sub.T, metric='cosine'))

ii_manual = []
for v in layers.values():
    ii_manual.extend(v)

ii = ii_manual

# %% Plot them
mods = hs_results_mod.Module.unique()
mods = [x for x in mods if x != -1 and not np.isnan(x)]
mods = sorted(mods)

sizes = [(hs_results_mod.Module == x).sum() for x in mods]
sizes = [np.log10(x) for x in sizes]

fig, axs = plt.subplots(len(mods), 1, gridspec_kw=dict(height_ratios=sizes),
                        figsize=(6, 8.5))

for ax, mod in zip(axs.ravel(), mods):

    plt.sca(ax)
    mod_genes = hs_results_mod.Symbol[hs_results_mod.Module == mod]
    mod_genes = pd.Index(mod_genes) & kernels.index

    sns.heatmap(kernels_z.loc[mod_genes].iloc[:, ii], vmin=-2, vmax=2,
                cmap="RdBu_r", ax=ax, xticklabels=True, rasterized=True)
    plt.ylabel(int(mod))

    plt.yticks([])

    if mod < max(mods):
        plt.xticks([])
    else:
        ax.tick_params(labelsize=7)

plt.tight_layout()
# plt.show()
plt.savefig('KernelHeatmap.svg', dpi=300)
