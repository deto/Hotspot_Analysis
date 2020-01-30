from scipy.cluster.hierarchy import linkage, leaves_list, fcluster
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
latent_file = snakemake.input["latent"] # Only used for N cells

MIN_CLUSTER_GENES = snakemake.params["min_cluster_genes"]  # 50
CORE_ONLY = snakemake.params["core_only"]  # 50
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

cdist = 1-reg.values
cdist = cdist/2 + cdist.T/2  
np.fill_diagonal(cdist, 0)

Z = linkage(squareform(cdist), method='average')


# Find critical T value
N_CELLS = pd.read_table(latent_file, index_col=0).shape[0]
from scipy.stats import t

regX = reg.values.copy()
regX = regX/2 + regX.T/2
np.fill_diagonal(regX, 0)
allR = squareform(regX)
allR = np.sort(allR)
allT = allR * np.sqrt((N_CELLS-2) / (1-allR**2))
allP = t.sf(allT, df=N_CELLS-2)
allP_c = multipletests(allP, method='fdr_bh')[1]
ii = np.nonzero(allP_c < FDR_THRESHOLD)[0][0]
Z_THRESHOLD = allR[ii]

print(Z_THRESHOLD)

import hotspot.modules
modules = hotspot.modules.assign_modules(
    Z, offset=1, Z_THRESHOLD=Z_THRESHOLD, MIN_THRESHOLD=MIN_CLUSTER_GENES,
    leaf_labels=reg.index)


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
