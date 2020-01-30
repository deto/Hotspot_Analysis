import loompy
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors
import seaborn as sns

from scipy.spatial.distance import cdist


# %%

loom_file = "../../data/Lineage/Embryo3/data.loom"
umap_file = "../Embryo3/umap/umap.txt"
kernels_file = "../../data/Lineage/GSE122187_CellStateKernels.xls"

# %% Load data

with loompy.connect(loom_file, 'r') as ds:
    barcodes = ds.ca['Barcode'][:]
    counts = ds[:, :]
    gene_info = ds.ra['EnsID', 'Symbol']
    num_umi = ds.ca['NumUmi'][:]

# Have to do this because data_slideseq makes it a numpy array
gene_info = pd.DataFrame(
    gene_info, columns=['EnsID', 'Symbol']).set_index('EnsID')
counts = pd.DataFrame(counts, index=gene_info.index, columns=barcodes)
num_umi = pd.Series(num_umi, index=barcodes)

umap = pd.read_table(umap_file, index_col=0)
umap = umap.loc[counts.columns]

kernels = pd.read_excel(kernels_file, skiprows=1).set_index('Gene')
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

# %% Create assignments

scaled = counts.divide(num_umi, axis=1)*10000

common = kernels.index & gene_info.Symbol
common_ens = gene_info.index[gene_info.Symbol.isin(kernels.index)]

sc_sub = scaled.loc[common_ens & scaled.index] \
    .join(gene_info) \
    .groupby('Symbol') \
    .sum()

X = np.log2(sc_sub+1)
Y = kernels

common = X.index & Y.index

X = X.loc[common]
Y = Y.loc[common]

dd = cdist(X.T.values, Y.T.values, metric='cosine')
dd = pd.DataFrame(dd, index=X.columns, columns=Y.columns)
dd.columns = [kernel_name_map[i] for i in dd.columns]

assignment = dd.idxmin(axis=1)
assignment = assignment.loc[counts.columns]

# %% Plot Assignments

colors = sns.color_palette("deep")
angio = colors[3]
rest = '#DDDDDD'
c = [angio if x == "Angioblasts" else rest for x in assignment]

cmap2 = matplotlib.colors.LinearSegmentedColormap.from_list(
    'grays', ['#DDDDDD', '#000000']
)

# Angio markers are: Pecam1, Cdh5, Tek, Lyve1

fig, axs = plt.subplots(1, 2, figsize=(8, 4))

cbar_ax = fig.add_axes(
    [.9, .1, .005, .2]
)

plt.sca(axs[0])
plt.scatter(
    umap.iloc[:, 0],
    umap.iloc[:, 1],
    s=2, c=c,
    rasterized=True
)
plt.title('Angioblasts')
for sp in plt.gca().spines.values():
    sp.set_visible(False)
plt.xticks([])
plt.yticks([])

marker_gene = 'Pecam1'
marker_ens = gene_info.index[gene_info.Symbol == marker_gene][0]

plt.sca(axs[1])

sc = plt.scatter(
    umap.iloc[:, 0],
    umap.iloc[:, 1],
    s=2,
    c=np.log2(scaled.loc[marker_ens]+1),
    cmap=cmap2, vmin=0, vmax=2,
    rasterized=True
)
plt.colorbar(sc, cax=cbar_ax, label='$log_2(CP10K+1)$')
plt.title(marker_gene)
for sp in plt.gca().spines.values():
    sp.set_visible(False)
plt.xticks([])
plt.yticks([])

# plt.show()
plt.savefig('Angioblast_UMAPS.svg', dpi=300)

# %% Plot Z vs. Z for Angioblast genes?

# kk = kernels.idxmax(axis=1)
# angio_genes = kk.index[kk == 32]
# 
# # Alternately:
# angio_genes = ['Pecam1', 'Cdh5', 'Tek', 'Lyve1', 'Gja5', 'Ephb4']
# 
# hs_res_tree = pd.read_table(
#     "../Embryo3/hotspot/hotspot_lineage_tree.txt",
#     index_col=0
# )
# 
# hs_res_tx = pd.read_table(
#     "../Embryo3/hotspot/hotspot_pca.txt",
#     index_col=0
# )
# 
# hs_res_tx = hs_res_tx.loc[hs_res_tree.index]
# 
# 
# ens_map = {s: ens for s, ens in zip(gene_info.Symbol, gene_info.index)}
# 
# angio_genes_ens = [ens_map[x] for x in angio_genes]
# 
# 
# plt.figure(figsize=(5, 5))
# 
# plt.plot(hs_res_tree.Z, hs_res_tx.Z, 'o', ms=2)
# plt.plot(hs_res_tree.Z[angio_genes_ens], hs_res_tx.Z[angio_genes_ens], 'o', ms=2)
# 
# plt.show()
