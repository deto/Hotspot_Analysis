from scipy.cluster.hierarchy import linkage, leaves_list, fcluster, cut_tree
from scipy.spatial.distance import squareform

import numpy as np
import pandas as pd
import __main__ as main
if hasattr(main, '__file__'):
    import matplotlib
    matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable
import seaborn as sns
from tqdm import tqdm
import hotspot.modules
from statsmodels.stats.multitest import multipletests

plt.rcParams["svg.fonttype"] = "none"

results_file = snakemake.input["results_z"]

MIN_CLUSTER_GENES = snakemake.params["min_cluster_genes"]  # 50
FDR_THRESHOLD = snakemake.params["fdr_threshold"]


cluster_heatmap_file = snakemake.output["cluster_heatmap"]
cluster_output_file = snakemake.output["cluster_output"]
linkage_output_file = snakemake.output["linkage_output"]

reg = pd.read_table(results_file, index_col=0)

# %% Compute Linkage and Ordering

# Make 'symmetric'
# Matrix is already symmetric but differences in lower-order bits
# are causing scipy to complain here
# Same with diagonal

cdist = reg.values.max()-reg.values
np.fill_diagonal(cdist, 0)
Z = linkage(squareform(cdist), method='average')

cut_height = reg.values.max() - np.percentile(reg.values.ravel(), 90)
modules = cut_tree(Z, height=cut_height)

modules = pd.Series(modules.ravel(), reg.index)

# Remove modules with too few genes
module_counts = modules.value_counts()
bad_modules = module_counts.index[module_counts < MIN_CLUSTER_GENES]

for bm in bad_modules:
    modules.loc[modules == bm] = -1

# Rename modules to be ascending
module_counts = modules.value_counts()
module_counts = module_counts.loc[module_counts.index != -1]

module_map = {x: int(i+1) for i, x in enumerate(module_counts.index)}
module_map[-1] = -1
modules = modules.map(module_map)

modules.rename('Cluster').to_frame().to_csv(cluster_output_file, sep="\t")
linkage_out = pd.DataFrame(Z).to_csv(
    linkage_output_file, header=False, index=False, sep="\t")


# %% Plot the clusters

from scipy.cluster.hierarchy import leaves_list
ii = leaves_list(Z)

colors = list(plt.get_cmap("tab10").colors)
row_colors1 = pd.Series(
    [colors[i % 10] if i != -1 else "#ffffff" for i in modules],
    index=reg.index,
)

row_colors = pd.DataFrame({"Cluster": row_colors1})

cm = sns.clustermap(
    reg.iloc[ii, ii],
    row_cluster=False,
    col_cluster=False,
    vmin=-15,
    vmax=15,
    cmap="RdBu_r",
    xticklabels=False,
    yticklabels=False,
    row_colors=row_colors.iloc[ii],
    rasterized=True,
)

fig = plt.gcf()
fig.patch.set_visible(False)
plt.sca(cm.ax_heatmap)
plt.ylabel("")
plt.xlabel("")
plt.savefig(cluster_heatmap_file, dpi=300)
# plt.show()
