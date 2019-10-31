import numpy as np
import pandas as pd
import hotspot.modules
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.colors import Normalize
from matplotlib.cm import ScalarMappable


# %% Load data

z_scores = pd.read_table(
    "../../CD4_w_protein/hotspot/hotspot_pairs_z.txt.gz",
    index_col=0
)

hs_results = pd.read_table(
    "../../CD4_w_protein/hotspot/hotspot_hvg.txt",
    index_col=0
)

ens_map = {x: y for x, y in zip(hs_results.index, hs_results.Symbol)}

# %% Compute Modules

modules, Z = hotspot.modules.compute_modules(
    z_scores, min_gene_threshold=10, core_only=False
)

# %% Modules to gene sets

gene_sets = {}
for i in modules.unique():
    if i == -1: continue
    genes = modules[modules == i].index
    genes = {ens_map[x].upper() for x in genes}
    gene_sets[i] = genes

all_genes = {ens_map[x].upper() for x in hs_results.index}

# %% Load GO sets

from gene_enrich import load_gene_set_gmt

go_sets = load_gene_set_gmt("/data/yosef2/users/david.detomaso/Signatures/GO/GO_biological_process.gmt")
go_sets2 = load_gene_set_gmt("/data/yosef2/users/david.detomaso/Signatures/Enrichr/GO_Biological_Process_2015.txt")

# The Enrichr one is better!

# %% Overlap analysis
from  gene_enrich.analyses import gene_set_enrichment


results = gene_set_enrichment(
    gene_sets[6], all_genes, go_sets2, correct_false_neg=True
)

# %% Why are 0 and 6 different?

g0 = modules.index[modules == 0]
g6 = modules.index[modules == 6]

z0 = z_scores.loc[g0, g0].values.ravel()
z6 = z_scores.loc[g6, g6].values.ravel()
z06 = z_scores.loc[g0, g6].values.ravel()


plt.figure()
kwargs = dict( 
    range=(-100, 100),
    alpha=0.5,
    bins=30,
)
plt.hist(z0, **kwargs, label='z0')
plt.hist(z6, **kwargs, label='z6')
plt.hist(z06, **kwargs, label='z06')
plt.legend()
plt.show()


zz06 = z_scores.loc[g0, g6]
